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
    "g33s" => "", "g88s" => "(1/3)", "g88conns" => "(1/3)", "gSSs" => "(1/9)", "gCCconns" => "(4/9)", "gCCdiscs" => "(4/9)", "gC8discs" => "(2/3\\sqrt{3})", "∆ls_amus" => "\\frac{1}{3}\\,", "∆ls_amuconns" => "", "∆lc_bs" => "(4/9)",
    "g3333" => 1.,  "g3388" => 2/3., "g33CC" => 8/9., "g8888" => 1/9., "g88CC" => 8/27., "gCCCC" => 16/81., 
    "g3333s" => "",  "g3388s" => "(2/3)", "g33CCs" => "(8/9)", "g8888s" => "(1/9)", "g88CCs" => "(8/27)", "gCCCCs" => "(16/81)",
)

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

# Set plot parameters

diag = "NLOa&b"  #  LO  NLOa  NLOb  NLOc  NLOa&b
wind = "SDsub"  #  NW  SD  SDsub  ID  LD  ILD
comp = "g33"  #  g33  g88  gSS  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆ls_amuconn  ∆lc_b

Q = 5.0  # virtuality for SDsub

BLIND = false

STD_DERIV = false
tl_IMPR   = true
VREF      = true
RESC      = false


path_bdio = path_bdio_dict["local"]
# path_bdio = joinpath(julia_script_directory,"..","..","..","HVP lepton mass","ObsBDIO")

# Data reading and definitions

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

@info("All data ready")


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Pion chiral-continiuum best extrapolation")

ShowGHOST = false

SAVE     = true
OVERSAVE = false

ARGPLOT = 1  #  set to 1 for best plot

IMPR_SET = IMPR_SET  #  readIMPR_SET  ["1"]  ["2"]  ["1","2"]  ["1old","2"]  ["1","1old","2"]

Pi_ph = !RESC ? phi2_ph : y_ph
K_ph  = !RESC ? phi4_ph : z_ph

myFactor = charge_factor[comp] * (diag == "LO" ? 1 : 10)

impr_set = "2"
key = mykeys[2]

for (i,impr_set) in enumerate([impr_set])
    for (j,key) in enumerate([key])
        label = []
        color = ["orange","red","purple","blue","green","brown"]
        fmt   = ["o","^","s","d","h","v"]

        argSort = sortperm(eachindex(info[impr_set]["average"]["weight"][key]), 
                    by=i -> vcat([label_tot_isov[FitCut] for FitCut in FITCUT]...)[i][1] in ["baseResc", "baseSimpResc"] ? -Inf : info[impr_set]["average"]["weight"][key][i], 
                    rev=true)
        argPlot = argSort[ARGPLOT]
        discr   = key[end-1:end]

        offset = 0
        newArg   = nothing
        myFitCut = nothing
        for FitCut in FITCUT
            len = length(label_tot_isov[FitCut])
            if offset < argPlot <= offset + len
                newArg = argPlot - offset
                myFitCut = FitCut
            end
            offset += len
        end
        argPlot  = newArg

        if myFitCut in ["beta","beta&mass"]
            color = color[2:end]
            fmt   = fmt[2:end]
        elseif myFitCut == "beta_ext"
            color = color[3:end]
            fmt   = fmt[3:end]
        end

        fit      = fitres[impr_set][myFitCut][key][argPlot]
        model    = label_tot_isov[myFitCut][argPlot]
        myres    = res[impr_set][myFitCut][key][argPlot]; uwerr(myres)

        # myparam  = param[impr_set][myFitCut][key][argPlot]; uwerr.(myparam)

        # fit to obtain params
        myfit, fitresid = fit_routine(f_tot_isov[myFitCut][argPlot], value.(xydata[myFitCut]["xdata"]), xydata[myFitCut]["ydata"][impr_set][key], n_par_tot_isov[myFitCut][argPlot], pval=true, info=false, lineprint=false, fitRes=true)
        myparam = myfit.param

        xproj_ch = deepcopy(xydata[myFitCut]["xdata"])
        xproj_ch[:,3] = fill(K_ph, length(xydata[myFitCut]["xdata"][:,1]))::Vector{uwreal}

        Deltay = f_tot_isov[myFitCut][argPlot](xproj_ch,myparam) - f_tot_isov[myFitCut][argPlot](xydata[myFitCut]["xdata"],myparam)
        ydata = myFactor * (xydata[myFitCut]["ydata"][impr_set][key] + Deltay); uwerr.(ydata)

        println("Chosen model (diag=$diag, set=$impr_set, discr=$discr, cut=$myFitCut): $(label_tot_isov[myFitCut][argPlot]) (n: $argPlot)")
        println("- chi2/dof     = $(fit.chi2/fit.dof)")
        println("- chi2/chi2exp = $(fit.chi2/fit.chi2exp)")
        println("- pval = $(myfit.pval)")

        # fig = figure(figsize=(10,7.5))
        fig = figure(figsize=(8,6))

        mpl = pyimport("matplotlib.lines")  # Import the `lines` module from Matplotlib
        Line2D = mpl.Line2D  # Get the Line2D class
        handles = []

        for (k,b) in  enumerate(sort(unique(getfield.(ensInfo[myFitCut], :beta))))
            # push!(label,L"$\beta=$"*"$b")
            push!(handles,Line2D([], [], 
                color=color[k], 
                marker=fmt[k],
                markerfacecolor="none",
                markeredgecolor=color[k],
                markersize=8, 
                linestyle="--", 
                label=L"$\beta=$"*"$b"))
            n_ = findall(x->x.beta == b, ensInfo[myFitCut])
            a2_aux = mean(value.(xydata[myFitCut]["xdata"][:,1][n_]))
            errorbar(value.(xydata[myFitCut]["xdata"][n_,2]), value.(ydata[n_]), err.(ydata[n_]), fmt=fmt[k], capsize=2, color=color[k], mfc="none")
            if ShowGHOST
                uwerr.(xydata[myFitCut]["ydata"][impr_set][key])
                errorbar(value.(xydata[myFitCut]["xdata"][n_,2]), myFactor*value.(xydata[myFitCut]["ydata"][impr_set][key][n_]), myFactor*err.(xydata[myFitCut]["ydata"][impr_set][key][n_]), fmt=fmt[k], capsize=2, color=color[k], mfc="none", alpha=0.1)
            end
            xxx = !RESC ? [fill(a2_aux,100) Float64.(range(0.039, 0.8, length=100)) fill(phi4_ph.mean, 100)] : [fill(a2_aux,100) Float64.(range(0.02, 0.4, length=100)) fill(z_ph.mean, 100)]
            yyy = myFactor*f_tot_isov[myFitCut][argPlot](xxx, myparam); uwerr.(yyy)
            # yyy = myFactor*f_tot_isov[myFitCut][argPlot](xxx, value.(myparam))
            PyPlot.plot(xxx[:,2],value.(yyy),ls="--",color=color[k],lw=0.5)
            fill_between(xxx[:,2], value.(yyy)-err.(yyy), value.(yyy)+err.(yyy), alpha=0.2, color=color[k])
            # PyPlot.plot(xxx[:,2],yyy,ls="--",color=color[k],lw=0.5)
        end
        xxx_ph = !RESC ? [fill(0.0,100) Float64.(range(0.039, 0.8, length=100)) fill(phi4_ph.mean, 100)] : [fill(0.0,100) Float64.(range(0.02, 0.4, length=100)) fill(z_ph.mean, 100)]
        yyy_ph = myFactor*f_tot_isov[myFitCut][argPlot](xxx_ph, myparam); uwerr.(yyy_ph)

        errorbar(Pi_ph.mean, myFactor*myres.mean, myFactor*myres.err, fmt="^", capsize=2, color="black",label="ph.")
        push!(label,"ph.")
        PyPlot.plot(xxx_ph[:,2],value.(yyy_ph),ls="--",color="gray",lw=0.5)
        fill_between(xxx_ph[:,2], value.(yyy_ph)-err.(yyy_ph), value.(yyy_ph)+err.(yyy_ph), alpha=0.2, color="gray")

        axvline(x=Pi_ph.mean, ls="dashed", color="black", lw=0.2, alpha=0.7) 
        if !RESC
            xlabel(L"$\Phi_2$")
        else
            xlabel(L"$y$")
        end
        diag_str = diag == "LO" ? "\\mathrm{lo}" : (diag == "NLOa&b" ? "\\mathrm{nlo(a\\&b)}" : "\\mathrm{nlo($(diag[end]))}")
        comp_str = comp[1] == '∆' ? (comp == "∆lc_b" ? "\\Delta_{lc}(\\tilde{b}^{$diag_str})" : comp == "∆ls_amu" ? ("\\Delta_{ls}\\left(a_{\\mu}^{$diag_str}\\right)") : ("\\Delta^{\\mathrm{conn.}}_{ls}\\left(a_{\\mu}^{$diag_str}\\right)")) : (diag != "NLOc" ? "\\mathrm{$(comp[2]),$(comp[3])}" : "\\mathrm{$(comp[2]),$(comp[3])-$(comp[4]),$(comp[5])}")
        if comp in ["gCCdisc","gC8disc"]; comp_str *= "\\mathrm{(disc)}"; end
        main_str  = comp[1] == '∆' ? comp_str : "a_{\\mu}^{$comp_str,\\,$diag_str}"
        V_str    = VREF ? "(V_{\\mathrm{ref}})" : ""
        wind_str = wind == "SDsub" ? (comp != "∆ls_amu" ? "^{\\mathrm{SD}}_{\\mathrm{sub}}" : "^{\\mathrm{SD}}") : "^{\\mathrm{$wind}}"
        Q_str    = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "($Q\\ \\mathrm{GeV})" : ""
        fact_str = diag == "LO" ? "\\times10^{10}" : "\\times10^{11}"
        BLIND_str = BLIND ? "\\mathrm{BLIND}\\times" : ""
        if wind == "NW" || comp == "∆lc_b"
            ylabel(latexstring("$(BLIND_str)$(charge_factor[comp*"s"]) $(main_str)$(V_str)$(Q_str)$(fact_str)"))
        else
            if comp != "∆ls_amu"
                ylabel(latexstring("$(BLIND_str)$(charge_factor[comp*"s"])\\left($(main_str)$(V_str)\\right)$(wind_str)$(Q_str)$(fact_str)"))
            else
                ylabel(latexstring("$(BLIND_str)$(charge_factor[comp*"s"])$(main_str)$(V_str)$(wind_str)$(Q_str)$(fact_str)"))
            end
        end
        if myFitCut in ["mass","beta&mass"]
            !RESC ? (wind == "LD" ? xlim(0.06,0.62) : xlim(0.04,0.62)) : xlim(0.02,0.31)
        else
            !RESC ? (wind == "LD" ? xlim(0.06,0.8) : xlim(0.04,0.8)) : xlim(0.02,0.4)
        end
        legend(handles=handles, loc="best", ncol=1)
        # legend(label,loc="best")
        tight_layout()
        display(gcf())
        if SAVE
            RESCstr = !RESC ? "" : "_resc"
            p = create_path(path_plot,["extr_$(diag)_$(wind)_$(key)_$(impr_set)$(RESCstr).pdf"],OVERWRITE=OVERSAVE)
            PyPlot.savefig(p)
        end
        close()
    end
end

##


using Plots
using Colors