# Import packages

using Revise

include("../../HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs

using ALPHAio, ADerrors
import ADerrors: err

using Statistics
using Distributions

using BDIO
using JLD2

using Plots
using PyPlot
using PyCall
using Colors

mpl_fold = pyimport("mpl_fold_axis")

using ProgressBars
using Suppressor

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "..", "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "..", "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "..", "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_bPert   = joinpath(julia_script_directory, "..", "..", "..", "PertSD")
path_FVCcont = joinpath(julia_script_directory, "..", "..", "..", "FVCcont")

path_plot    = joinpath(julia_script_directory)

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

# Charge factor

charge_factor = Dict(
    "g33" => 1., "g88" => 1/3., "g88conn" => 1/3., "gSS" => 1/9., "gCCconn" => 4/9., "gCCdisc" => 4/9., "gC8disc" => 2/(3*sqrt(3)), "∆ls_amu" => 1/3., "∆ls_amuconn" => 1., "∆lc_b" => 4/9.,
    "g33s" => "", "g88s" => "(1/3)", "g88conns" => "(1/3)", "gSSs" => "(1/9)", "gCCconns" => "(4/9)", "gCCdiscs" => "(4/9)", "gC8discs" => "(2/3\\sqrt{3})", "∆ls_amus" => "(1/3)", "∆ls_amuconns" => "", "∆lc_bs" => "(4/9)",
    "g3333" => 1.,  "g3388" => 2/3., "g33CC" => 8/9., "g8888" => 1/9., "g88CC" => 8/27., "gCCCC" => 16/81., 
    "g3333s" => "",  "g3388s" => "(2/3)", "g33CCs" => "(8/9)", "g8888s" => "(1/9)", "g88CCs" => "(8/27)", "gCCCCs" => "(16/81)",
)

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

FITdata = true

# Set plot parameters

diag = "NLOc"  #  LO  NLOa  NLOb  NLOc  NLOa&b
wind = "NW"  #  NW  SD  SDsub  ID  LD  ILD
COMP = ["g3333","g3388","g33CC"]  #  g33  g88  gSS  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆ls_amuconn  ∆lc_b

Q = 5.0  # virtuality for SDsub

BLIND     = false

STD_DERIV = false
tl_IMPR   = false
VREF      = false
RESC      = false

SAVE      = true
OVERSAVE  = false

SHIFT_t0  = true

path_bdio = path_bdio_dict["local"]
# path_bdio = joinpath(julia_script_directory,"..","..","..","HVP lepton mass","ObsBDIO")

# Data reading and definitions

# fig = figure(figsize=(8,10))
fig = figure(figsize=(10,10))

funcOrder_conv = [
    "a2","a2loga","a3","a4",
    "a2phi2","a2phi4","a3phi2","phi2","phi2sqr","phi2log","phi2inv","logphi2","phi4","phi4sqr","phi4log","phi4inv","logphi4",
    "a2y","a2z","a3y","y","ysqr","ylog","yinv","logy","z","zsqr","zlog","zinv","logz"
]
order_map = Dict(s => i for (i, s) in enumerate(funcOrder_conv))
function func_str(func_vec::Vector{Function};Order::Bool=true)::String
    str_vec = string.(func_vec)
    if Order
        str_vec = sort(str_vec, by = x -> order_map[x])
    end
    return paste_str(str_vec)
end

for comp in COMP
    if wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]
        @info("Reading obsBDIO data:\n - wind: $wind\n - diag: $diag\n - comp: $comp\n - Q: $Q")
    else
        @info("Reading obsBDIO data:\n - wind: $wind\n - diag: $diag\n - comp: $comp")
    end

    println("- Reading FINAL result...")

    RES, INFO = BDIOread_MAtot(path_bdio,diag,wind,comp,resc=RESC,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q); uwerr(RES)
    SYST = INFO["syst"]

    FITCUT   = sort(collect(keys(INFO["FITCUTtoMODEL"])), by = x -> Dict(s => i for (i, s) in enumerate(["None","beta","mass","beta&mass","beta_ext"]))[x])
    MultFunc = INFO["MultFunc"]
    IMPR_SET = comp in ["gCCdisc","gC8disc"] ? [""] : INFO["IMPR_SET"]
    mykeys   = INFO["Keys"]

    println("- Reading info...")
    # MA = Dict()
    info = Dict()
    for impr_set in IMPR_SET
        println("   - Impr. set $impr_set")
        # MA[impr_set] = Dict()
        info[impr_set] = Dict()
        for FitCut in FITCUT
            println("      - Fit Cut $FitCut")

            # MA[impr_set][FitCut], info[impr_set][FitCut] = BDIOread_MA(path_bdio,diag,wind,comp,INFO["FITCUTtoMODEL"][FitCut][2],FitCut,impr_set,resc=RESC,SimpleBase=INFO["FITCUTtoMODEL"][FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)
            model_str = func_str(INFO["FITCUTtoMODEL"][FitCut][2])
            tl_str = (tl_IMPR && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
            SIMstr = INFO["FITCUTtoMODEL"][FitCut][1] ? "SIMPLE" : ""
            MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
            SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
            IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
            VREFstr = VREF ? "_Vref" : ""
            RESstr = RESC ? "_resc" : ""
            DERstr = STD_DERIV ? "_std" : "" 
            BLINstr = BLIND ? "Blind" : ""
            pJDL2 = joinpath(path_bdio,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCut)]","MA","$(BLINstr)MAinfo$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jld2")
            info[impr_set][FitCut] = load(pJDL2,"MAinfo")
        end
        println("      - Fit Cut Average")
        # MA[impr_set]["average"], info[impr_set]["average"] = BDIOread_MAtot(path_bdio,diag,wind,comp,read="impr",resc=RESC,impr_set=impr_set,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)
        tl_str = (tl_IMPR && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
        SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
        VREFstr = VREF ? "_Vref" : ""
        RESstr = RESC ? "_resc" : ""
        DERstr = STD_DERIV ? "_std" : "" 
        BLINstr = BLIND ? "Blind" : ""
        IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : "_set"
        pJDL2 = joinpath(path_bdio,"Fit&MA",diag,wind,comp*tl_str,"MAtot","$(BLINstr)MAinfo$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jld2")
        info[impr_set]["average"] = load(pJDL2,"MAinfo")
    end

    if FITdata
        println("- Reading FITs data...")

        xdata = Dict(); ydata = Dict()
        fitres = Dict()
        res =  Dict(); param = Dict()
        modelinfo = Dict()
        for impr_set in IMPR_SET
            println("   - Impr. set $impr_set")
            
            fitres[impr_set] = Dict()
            res[impr_set]    =  Dict()
            param[impr_set]  = Dict()
            ydata[impr_set]  = Dict()
            for FitCut in FITCUT
                println("      - Fit Cut $FitCut")

                println("         - Reading X & Y data...")
                xdata[FitCut], ydata[impr_set][FitCut] = BDIOread_XYdata(path_bdio,diag,wind,comp,INFO["FITCUTtoMODEL"][FitCut][2],FitCut,impr_set,resc=RESC,SimpleBase=INFO["FITCUTtoMODEL"][FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

                println("         - Reading FitRes...")

                fitres[impr_set][FitCut] = JDL2read_FitRes(path_bdio,diag,wind,comp,INFO["FITCUTtoMODEL"][FitCut][2],FitCut,impr_set,resc=RESC,SimpleBase=INFO["FITCUTtoMODEL"][FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

                println("         - Reading results and parameters...")

                res[impr_set][FitCut], param[impr_set][FitCut] = BDIOread_res(path_bdio,diag,wind,comp,INFO["FITCUTtoMODEL"][FitCut][2],FitCut,impr_set,resc=RESC,SimpleBase=INFO["FITCUTtoMODEL"][FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,Q=Q,BLIND=BLIND,param=true)

                println("         - Reading model information...")

                modelinfo[FitCut] = JDL2read_ModelInfo(path_bdio,diag,wind,comp,INFO["FITCUTtoMODEL"][FitCut][2],FitCut,impr_set,resc=RESC,SimpleBase=INFO["FITCUTtoMODEL"][FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)
            end
        end

        xydata = Dict()

        model_len = Dict()
        a2RESCAL  = Dict()
        nens      = Dict()
        ensInfo   = Dict()

        f_tot_isov = Dict()
        n_par_tot_isov = Dict()
        label_tot_isov = Dict()

        for FitCut in FITCUT
            xydata[FitCut] = Dict{String,Any}("xdata" => xdata[FitCut])
            xydata[FitCut]["ydata"] = Dict(); [xydata[FitCut]["ydata"][impr_set] = ydata[impr_set][FitCut] for impr_set in IMPR_SET]

            model_len[FitCut] = modelinfo[FitCut]["length"]
            a2RESCAL[FitCut]  = modelinfo[FitCut]["a2Rescaling"]
            nens[FitCut]      = modelinfo[FitCut]["nens"]
            ensInfo[FitCut]   = EnsInfo.(modelinfo[FitCut]["ensList"])

            f_tot_isov[FitCut], n_par_tot_isov[FitCut], label_tot_isov[FitCut] = call_models(INFO["FITCUTtoMODEL"][FitCut][2],ensInfo[FitCut],4,SimpleBase=INFO["FITCUTtoMODEL"][FitCut][1],a2resc=a2RESCAL[FitCut],MultFunc=MultFunc,fPiresc=RESC,na_max=modelinfo[FitCut]["na_max"],nmPi_max=modelinfo[FitCut]["nmPi_max"],nmK_max=modelinfo[FitCut]["nmK_max"])
        end
    end


    @info("Continiuum extrapolation plot")

    # IMPR_SET = readIMPR_SET  #  readIMPR_SET  ["1"]  ["2"]  ["1","2"]  ["1old","2"]  ["1","1old","2"]

    npoints = 50
    Wdump   = 1e5

    println("   - Computing 'x' and 'y' plot points")

    xarr  = Dict()
    xproj = Dict()
    for FitCut in FITCUT
        xarr[FitCut] = [Float64.(range(1e-5, maximum(value.(xydata[FitCut]["xdata"][:,1]))+0.005, length=npoints)) fill(value(phi2_ph), npoints) fill(value(phi4_ph), npoints)]; yarr = Dict()
    end
    yproj = Dict(); yproj_syst = Dict()
    yarr  = Dict()
    weig  = Dict()
    for (j,impr_set) in enumerate(IMPR_SET)

        weig[impr_set]  = Dict()
        yarr[impr_set]  = Dict()
        yproj[impr_set] = Dict(); yproj_syst[impr_set] = Dict()
        for (k,key) in enumerate(mykeys)
            yarr[impr_set][key] = Dict()
            weig[impr_set][key] = Dict()

            # valid_indices = Dict()
            # valid_indices = filter(i -> vcat([label_tot_isov[FitCut] for FitCut in FITCUT]...)[i][1] ∉ ["baseResc","baseSimpResc"], eachindex(info[impr_set]["average"]["weight"][key]))

            # warg_average = sortperm(info[impr_set]["average"]["weight"][key],rev=true)
            # warg_average = warg_average[[warg in valid_indices for warg in warg_average]]

            Wcut = maximum(info[impr_set]["average"]["weight"][key])/Wdump
            
            i = 0
            for FitCut in FITCUT
                # yarr[impr_set][key][FitCut] = Vector{Vector{uwreal}}()
                yarr[impr_set][key][FitCut] = Vector{Vector{Float64}}()
                weig[impr_set][key][FitCut] = Vector{Float64}()

                println("- Impr. set = $impr_set; Comp = $key; Fit cut = $FitCut")
                println("   - [Computing bands...]")

                valid_indices = filter(i -> label_tot_isov[FitCut][i][1] ∉ ["baseResc","baseSimpResc"], eachindex(info[impr_set][FitCut]["weight"][key]))

                # warg[impr_set][key][FitCut] = sortperm(info[impr_set][FitCut]["weight"][key][valid_indices[FitCut]]; rev=true) .|> i -> valid_indices[FitCut][i]
                # warg[impr_set][key][FitCut] = warg[impr_set][key][FitCut][info[impr_set][FitCut]["weight"][key][warg[impr_set][key][FitCut]] .> info[impr_set][FitCut]["weight"][key][warg[impr_set][key][FitCut][1]]/100]
                
                len = length(info[impr_set][FitCut]["weight"][key])
                W = info[impr_set]["average"]["weight"][key][(i+1):(i+len)][valid_indices]
                for (i,ind) in ProgressBar(enumerate(valid_indices))
                    if W[i] > Wcut
                        # my = f_tot_isov[FitCut][ind](xarr[FitCut], param[impr_set][FitCut][key][ind]); uwerr.(my)
                        my = f_tot_isov[FitCut][ind](xarr[FitCut], value.(param[impr_set][FitCut][key][ind]))
                        push!(yarr[impr_set][key][FitCut], my)
                        push!(weig[impr_set][key][FitCut], W[i])
                    end
                end
                i += len
            end
        end
    end

    # Shift to new t0
    if SHIFT_t0
        der_sqrtt0 = mchist(RES, "sqrtt0 [fm]")[1] / artificial_err

        RES += der_sqrtt0 * (0.1440 - 0.1442); uwerr(RES)
        for impr_set in IMPR_SET
            for key in mykeys
                for FitCut in FITCUT
                    for i=1:length(yarr[impr_set][key][FitCut])
                        yarr[impr_set][key][FitCut][i] .+= der_sqrtt0 * (0.1440 - 0.1442)
                    end
                end
            end
        end
    end

    wPen = 1.0

    ls_beta_ext = (0, (5, 2))      # long dashed
    color_dict = Dict(
        ""     => Dict("cc" => "gold"),
        "1"    => Dict("ll" => "blue",  "lc" => "red"  ),
        "2"    => Dict("ll" => "green", "lc" => "brown")
    ); color_dict["1old"] = color_dict["1"]
    ls_dict = Dict(
        "None"      => "solid",
        "mass"      => "dashed",
        "beta"      => "dotted",
        "beta&mass" => "dashdot",
        "beta_ext" => ls_beta_ext
    )
    DictFITCUTtoSTR = Dict(
        "None" =>  latexstring("\\mathrm{None}"),
        "beta" => latexstring("\\beta>3.34"), 
        "mass" => latexstring("m_\\pi<400"),
        "beta&mass" => latexstring("\\beta>3.34 \\& m_\\pi<400"),
        "beta_ext" => latexstring("\\beta>3.4"),
        )

    myFactor = charge_factor[comp] * (diag == "LO" ? 1 : (diag == "NLOc" ? 1000 : 10))

    for (j,impr_set) in enumerate(IMPR_SET)
        for (k,key) in enumerate(mykeys)
            for (f,FitCut) in enumerate(FITCUT)
                for i=1:length(yarr[impr_set][key][FitCut])
                    # fill_between(xarr[FitCut][:,1], value.(yarr[impr_set][key][FitCut][i]).-err.(yarr[impr_set][key][FitCut][i]), value.(yarr[impr_set][key][FitCut][i]).+err.(yarr[impr_set][key][FitCut][i]), alpha=(!a2RESCAL[FitCut] ? 1 : 2)*info[impr_set][FitCut]["weight"][key][warg[impr_set][key][FitCut][i]]*wPen, color=color_dict[impr_set]["$(key[end-1:end])"])
                    # PyPlot.plot(xarr[FitCut][:,1], myFactor*value.(yarr[impr_set][key][FitCut][i]), alpha=(weig[impr_set][key][FitCut][i])^(2/3)*wPen, linestyle = ls_dict[FitCut], color=color_dict[impr_set]["$(key[end-1:end])"])
                    PyPlot.plot(xarr[FitCut][:,1], myFactor*yarr[impr_set][key][FitCut][i], alpha=(weig[impr_set][key][FitCut][i])^(2/3)*wPen, linestyle = ls_dict[FitCut], color=color_dict[impr_set]["$(key[end-1:end])"])
                end
            end
        end
    end
    errorbar([0.0],myFactor*value.(RES),myFactor*err.(RES),fmt="o",mfc="none",color="black", ms=5, capsize=3)
    errorbar([0.0],myFactor*value.(RES),myFactor*sqrt.(err.(RES).^2 .+ (SYST)^2),fmt="o",color="black", ms=5, capsize=3)
    # digits = comp in ["gCCdisc","gC8disc"] ? 7 : 
    # res_str = SHOWRES ? "  [$(round(value(amu[1]),digits=digits))($(round(err(amu[1]),digits=digits)))($(round(syst,digits=digits)))]" : ""
    # PyPlot.title("Projection to continuum extrapolation$res_str")
end
axvline(0.0, color="black", lw=0.2, alpha=0.8)
for beta in b_values
    axvline(value(1/(8 * t0sym(beta))) ,ls="dotted", color="black", lw=0.2, alpha=0.7)
end
xlabel(L"$a^2/8t_0$")
diag_str = diag == "NLOa&b" ? "NLOa\\&b" : diag
comp_str = diag == "NLOc" ? "d,e-f,g" : "d,e"
wind_str = wind == "SDsub" ? "^{\\mathrm{SD}}_{\\mathrm{sub}}" : "^{\\mathrm{$wind}}"
V_str    = VREF ? "(V_{\\mathrm{ref}})" : ""
fact_str = diag == "LO" ? "\\times10^{10}" : (diag == "NLOc" ? "\\times10^{13}" : "\\times10^{11}")
BLIND_str = BLIND ? "\\mathrm{BLIND}\\times" : ""
# if wind == "NW"
#     ylabel(latexstring("$(BLIND_str)a_{\\mu}^{\\mathrm{hvp}}[\\mathrm{$(diag_str)}]$(V_str)$fact_str"))
# else
#     ylabel(latexstring("$(BLIND_str)\\left(a_{\\mu}^{\\mathrm{hvp}}[\\mathrm{$(diag_str)}]$(V_str)\\right)$wind_str$fact_str"))
# end
ylabel(latexstring("$(BLIND_str)a_{\\mu}^{\\mathrm{hvp,\\,nlo(c)}}$(V_str)$fact_str"))
mpl = pyimport("matplotlib.lines")  # Import the `lines` module from Matplotlib
Line2D = mpl.Line2D  # Get the Line2D class
handles = []
ls_beta_ext = (0, (5, 2))      # long dashed
color_dict = Dict(
    ""     => Dict("cc" => "gold"),
    "1"    => Dict("ll" => "blue",  "lc" => "red", "ll-ll" => "blue", "lc-lc" => "red"),
    "2"    => Dict("ll" => "green", "lc" => "brown", "ll-ll" => "green", "lc-lc" => "brown")
); color_dict["1old"] = color_dict["1"]
ls_dict = Dict(
    "None"      => "solid",
    "mass"      => "dashed",
    "beta"      => "dotted",
    "beta&mass" => "dashdot",
    "beta_ext" => ls_beta_ext
)
DictFITCUTtoSTR = Dict(
    "None" =>  latexstring("\\mathrm{None}"),
    "beta" => latexstring("\\beta>3.34"), 
    "mass" => latexstring("m_\\pi<400"),
    "beta&mass" => latexstring("\\beta>3.34 \\& m_\\pi<400"),
    "beta_ext" => latexstring("\\beta>3.4"),
    )
for impr_set in ["1","2"]
    for (k,key) in enumerate(diag == "NLOc" ? ["ll-ll","lc-lc"] : ["ll","lc"])
        push!(handles,Line2D([], [], color=color_dict[impr_set][key], linestyle="-", label="discr. $(key); set $impr_set"))
    end
end
for FitCut in ["None","beta","mass","beta&mass"]
    push!(handles,Line2D([], [], color="gray", linestyle=ls_dict[FitCut], label=DictFITCUTtoSTR[FitCut]))
end
# legend(handles=handles, loc="upper left", bbox_to_anchor=(0.2, 0.35), ncol=2)
legend(handles=handles, loc="upper left", bbox_to_anchor=(0.25, 0.35), ncol=2)
ax = gca()
mpl_fold.fold_axis(ax, [(120, 230, 0.05)], axis="y", which="both")
ax.set_yticks([0, 25, 50, 75, 100, 250, 275])
ax.set_yticklabels(["0", "25", "50", "75", "100", "250", "275"])
tight_layout()
xlim(-0.001,0.06)
ylim(0.0,280)
display(gcf())
if SAVE
    RESCstr  = !RESC ? "" : "_resc"
    SHIFTstr = !SHIFT_t0 ? "" : "_shift"
    p = create_path(path_plot,["NLOc$(RESCstr)$(SHIFTstr).pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p)
end
close()
