# Import packages

using Revise

include("../HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs

using ALPHAio, ADerrors
import ADerrors: err

using Revise

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

# include uwreal constants

include("../HVPtool/uwConst.jl")

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_bPert   = joinpath(julia_script_directory, "..", "..", "PertSD")
path_FVCcont = joinpath(julia_script_directory, "..", "..", "FVCcont")

path_plot    = joinpath(julia_script_directory,"..","..","Slides & Plots","Plots")

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

# Charge factor

factor = Dict(
    "g33" => 1., "g88" => 1/3., "gCCconn" => 4/9., "∆ls_amu" => 1/3., "∆lc_b" => 4/9.,
    "g33s" => "", "g88s" => "(1/3)", "gCCconns" => "(4/9)", "∆ls_amus" => "(1/3)", "∆lc_bs" => "(4/9)"
)

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

FITdata = true

# Set plot parameters

diag = "NLOa&b"  #  LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
wind = "ID"  #  NW  SD  SDsub  ID  LD  ILD
comp = "g33"  #  g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b

Q = 5.0  # virtuality for SDsub

BLIND = false

STD_DERIV  = false
tl_IMPR    = false
VREF       = true
RESC       = false

path_bdio = path_bdio_dict["local"]


# Data reading and definitions

if wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]
    @info("Reading obsBDIO data:\n - wind: $wind\n - diag: $diag\n - comp: $comp\n - Q: $Q")
else
    @info("Reading obsBDIO data:\n - wind: $wind\n - diag: $diag\n - comp: $comp")
end

println("- Reading FINAL result...")

RES, INFO = BDIOread_MAtot(path_bdio,diag,wind,comp,resc=RESC,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q); uwerr(RES)

FITCUT   = sort(collect(keys(INFO["FITCUTtoMODEL"])), by = x -> Dict(s => i for (i, s) in enumerate(["None","beta","mass","beta&mass"]))[x])
MultFunc = INFO["MultFunc"]
IMPR_SET = INFO["IMPR_SET"]
mykeys   = INFO["Keys"]

println("- Reading PARTIAL result...")
MA = Dict(); info = Dict()
for impr_set in IMPR_SET
    println("   - Impr. set $impr_set")
    MA[impr_set] = Dict(); info[impr_set] = Dict()
    for FitCut in FITCUT
        println("      - Fit Cut $FitCut")

        MA[impr_set][FitCut], info[impr_set][FitCut] = BDIOread_MA(path_bdio,diag,wind,comp,INFO["FITCUTtoMODEL"][FitCut][2],FitCut,impr_set,resc=RESC,SimpleBase=INFO["FITCUTtoMODEL"][FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)
    end
    println("      - Fit Cut Average")
    MA[impr_set]["average"], info[impr_set]["average"] = BDIOread_MAtot(path_bdio,diag,wind,comp,read="impr",resc=RESC,impr_set=impr_set,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)
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
        res[impr_set] =  Dict(); param[impr_set] = Dict()
        for FitCut in FITCUT
            println("      - Fit Cut $FitCut")

            println("         - Reading X & Y data..")
            ydata[FitCut] = Dict()
            xdata[FitCut], ydata[FitCut][impr_set] = BDIOread_XYdata(path_bdio,diag,wind,comp,INFO["FITCUTtoMODEL"][FitCut][2],FitCut,impr_set,resc=RESC,SimpleBase=INFO["FITCUTtoMODEL"][FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

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
        xydata[FitCut] = Dict("xdata" => xdata[FitCut], "ydata" => ydata[FitCut])

        model_len[FitCut] = modelinfo[FitCut]["length"]
        a2RESCAL[FitCut]  = modelinfo[FitCut]["a2Rescaling"]
        nens[FitCut]      = modelinfo[FitCut]["nens"]
        ensInfo[FitCut]   = EnsInfo.(modelinfo[FitCut]["ensList"])

        f_tot_isov[FitCut], n_par_tot_isov[FitCut], label_tot_isov[FitCut] = call_models(INFO["FITCUTtoMODEL"][FitCut][2],ensInfo[FitCut],4,SimpleBase=INFO["FITCUTtoMODEL"][FitCut][1],a2resc=a2RESCAL[FitCut],MultFunc=MultFunc,fPiresc=RESC,na_max=modelinfo[FitCut]["na_max"],nmPi_max=modelinfo[FitCut]["nmPi_max"],nmK_max=modelinfo[FitCut]["nmK_max"])
    end
end

@info("All data ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##


@info("Continiuum 'magic' plot")

# IMPR_SET = readIMPR_SET  #  readIMPR_SET  ["1"]  ["2"]  ["1","2"]  ["1old","2"]  ["1","1old","2"]

npoints = 30

println("   - Computing 'x' and 'y' plot points")

xarr  = Dict()
xproj = Dict()
for FitCut in FITCUT
    xarr[FitCut] = [Float64.(range(1e-5, 1.5*maximum(value.(xydata[FitCut]["xdata"][:,1])), length=npoints)) fill(value(phi2_ph), npoints) fill(value(phi4_ph), npoints)]; yarr = Dict()
end
yproj = Dict(); yproj_syst = Dict()
yarr = Dict()
warg = Dict()
for (j,impr_set) in enumerate(IMPR_SET)

    warg[impr_set]  = Dict()
    yarr[impr_set]  = Dict()
    yproj[impr_set] = Dict(); yproj_syst[impr_set] = Dict()
    for (k,key) in enumerate(mykeys)
        yarr[impr_set][key] = Dict()
        warg[impr_set][key] = Dict()

        valid_indices = Dict()
        for FitCut in FITCUT
            yarr[impr_set][key][FitCut] = Vector{Vector{uwreal}}()
            println("- Impr. set = $impr_set; Comp = $key; Fit cut = $FitCut")
            println("   - [Computing bands...]")

            valid_indices[FitCut] = filter(i -> label_tot_isov[FitCut][i][1] ∉ ["baseResc","baseSimpResc"], eachindex(info[impr_set][FitCut]["weight"][key]))

            warg[impr_set][key][FitCut] = sortperm(info[impr_set][FitCut]["weight"][key][valid_indices[FitCut]]; rev=true) .|> i -> valid_indices[FitCut][i]
            warg[impr_set][key][FitCut] = warg[impr_set][key][FitCut][info[impr_set][FitCut]["weight"][key][warg[impr_set][key][FitCut]] .> info[impr_set][FitCut]["weight"][key][warg[impr_set][key][FitCut][1]]/10]

            for i in ProgressBar(warg[impr_set][key][FitCut])
                my = f_tot_isov[FitCut][i](xarr[FitCut], param[impr_set][FitCut][key][i]); uwerr.(my)
                push!(yarr[impr_set][key][FitCut], my)
            end
        end
    end
end

#-- Stop for only plot compilation

SAVE     = false
OVERSAVE = false

# IMPR_SET = readIMPR_SET

# mykeys = [""]

wPen = 0.1

# SHOWRES = false

# color_list = comp ∉ ["cc conn","cc disc","c8 disc"] ? ["blue","red","green","brown"] : ["blue","green"]
color_dict = Dict(
    ""     => Dict("cc" => "gold"),
    "1"    => Dict("ll" => "blue",  "lc" => "red"  ),
    "2"    => Dict("ll" => "green", "lc" => "brown")
); color_dict["1old"] = color_dict["1"]
ls_dict = Dict(
    "None"      => "solid",
    "mass"      => "dashed",
    "beta"      => "dotted",
    "beta&mass" => "dashdot"
)

fig = figure(figsize=(10,7.5))
for (j,impr_set) in enumerate(IMPR_SET)
    for (k,key) in enumerate(mykeys)
        for (f,FitCut) in enumerate(FITCUT)
            for i=1:length(yarr[impr_set][key][FitCut])
                # fill_between(xarr[FitCut][:,1], value.(yarr[impr_set][key][FitCut][i]).-err.(yarr[impr_set][key][FitCut][i]), value.(yarr[impr_set][key][FitCut][i]).+err.(yarr[impr_set][key][FitCut][i]), alpha=(!a2RESCAL[FitCut] ? 1 : 2)*info[impr_set][FitCut]["weight"][key][warg[impr_set][key][FitCut][i]]*wPen, color=color_dict[impr_set]["$(key[end-1:end])"])
                PyPlot.plot(xarr[FitCut][:,1], value.(yarr[impr_set][key][FitCut][i]), alpha=3*(!a2RESCAL[FitCut] ? 1 : 2)*info[impr_set][FitCut]["weight"][key][warg[impr_set][key][FitCut][i]]*wPen, linestyle = ls_dict[FitCut], color=color_dict[impr_set]["$(key[end-1:end])"])
            end
        end
    end
end
SYST = INFO["syst"]
errorbar([0.0],value.(RES),err.(RES),fmt="o",mfc="none",color="black", ms=5, capsize=3)
errorbar([0.0],value.(RES),sqrt.(err.(RES).^2 .+ (SYST)^2),fmt="o",color="black", ms=5, capsize=3)
digits = comp in ["gCCdisc","gC8disc"] ? 7 : 
# res_str = SHOWRES ? "  [$(round(value(amu[1]),digits=digits))($(round(err(amu[1]),digits=digits)))($(round(syst,digits=digits)))]" : ""
# PyPlot.title("Projection to continuum extrapolation$res_str")
axvline(0.0, color="black", lw=0.2, alpha=0.8)
for beta in b_values
    axvline(value(1/(8 * t0sym(beta))) ,ls="dotted", color="black", lw=0.2, alpha=0.7)
end
xlabel(L"$a^2/8t_0$")
if diag != "LO" # we multiply the y axis by a factor 10
    formatter(x, pos) = string(round(10 * x, digits=2))  # Round to 2 decimal places
    ax = gca()
    ax.yaxis.set_major_formatter(PyPlot.matplotlib.ticker.FuncFormatter(formatter))
end
# xlim(right=0.065)
# ylim(top=-32/10)
# ylim(top=-40/10)
# ylim(top=-60/10)
# ylim(bottom=27/10)
diag_str = diag == "LO" ? "\\rm{LO}" : (diag == "NLOa&b" ? "\\rm{NLOa\\&b}" : "\\rm{NLO$(diag[end])}")
comp_str = comp[1] == '∆' ? (comp[end] != 'b' ? "\\Delta_{ls}(a_{\\mu})" : "\\Delta_{lc}(b)") : "\\mathrm{$(comp[2]),$(comp[3])}"
Q_str    = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "($Q\\ \\rm{GeV})" : ""
fact_str = diag == "LO" ? "\\times10^{10}" : "\\times10^{11}"
if wind == "NW"
    ylabel(latexstring("a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]$Q_str$fact_str"))
else
    ylabel(latexstring("\\left(a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]\\right)^{\\rm{$wind}}$Q_str$fact_str"))
end
mpl = pyimport("matplotlib.lines")  # Import the `lines` module from Matplotlib
Line2D = mpl.Line2D  # Get the Line2D class
handles = []
for impr_set in IMPR_SET
    for key in mykeys
        push!(handles,Line2D([], [], color=color_dict[impr_set]["$(key[end-1:end])"], linestyle="-", label="discr. $(key[end-1:end]); set $impr_set"))
    end
end
legend(handles=handles, loc="upper center")
tight_layout()
display(gcf())
if SAVE
    RESCstr = !RESC ? "" : "_resc"
    p = create_path(path_plot,[diag,wind,"$(diag)_$(comp)_ContExtr.$(RESCstr).pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p)
end
close()


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Ensemble cut convergence plot")

SAVE     = true
OVERSAVE = false

sqrtt0_ph_TAR = nothing  # 0.1443  sqrtt0_ph_CLS  nothing

path_bdio = path_bdio_dict["local"]


DictFITCUTtoSTR = Dict(
    "None" =>  latexstring("\\rm{None}"),
    "beta" => latexstring("\\beta>3.34"), 
    "mass" => latexstring("m_\\pi<400"),
    "beta&mass" => latexstring("\\beta>3.34 \\& m_\\pi<400"))

myFactor = factor[comp] * (diag == "LO" ? 1 : 10)

myRES = myFactor * RES; uwerr(myRES)
myERR = myFactor * sqrt(INFO["syst"]^2 + err(RES)^2)

fig, (ax1, ax2) = subplots(1, 2,
    gridspec_kw = Dict("width_ratios" => [6, 1], "wspace" => 0),
    figsize = (12, 6)
)

Expec_RES = nothing  # nothing  [362.0,3.7,2.7]  [378.7,3.7,3.1]

y = 0
yTicksPos = []
yTicks    = []
groupLabel = []
groupPos   = []
for (i,impr_set) in enumerate(IMPR_SET)
    for (k,key) in enumerate(mykeys)
        argw0 = 1
        for FitCut in FITCUT
            y += 1

            res  = myFactor * MA[impr_set][FitCut]["res"][key]; uwerr(res)
            syst = myFactor * info[impr_set][FitCut]["syst"][key]

            if !isnothing(sqrtt0_ph_TAR)
                der = mchist(res, "sqrtt0 [fm]")[1] / artificial_err
                res_t0shift = res + value(sqrtt0_ph_TAR - sqrtt0_ph) * der; uwerr(res_t0shift)
            else
                res_t0shift = res
            end

            ax1.errorbar(value(res), xerr=err(res), y, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2, alpha=0.4)
            ax1.errorbar(value(res), xerr=sqrt(err(res)^2+syst^2), y, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2, alpha=0.4)
            if !isnothing(sqrtt0_ph_TAR)
                ax1.errorbar(value(res_t0shift), xerr=err(res_t0shift), y+0.1, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2, alpha=0.4)
                ax1.errorbar(value(res_t0shift), xerr=sqrt(err(res_t0shift)^2+syst^2), y+0.1, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2, alpha=0.4)
            end
            argwf = argw0 + length(info[impr_set][FitCut]["weight"][key]) - 1
            ax2.barh(y, sum(info[impr_set]["average"]["weight"][key][argw0:argwf]), height = 0.7, color = "orange", alpha = 0.9)
            argw0 = argwf + 1

            push!(yTicksPos, y)
            push!(yTicks, DictFITCUTtoSTR[FitCut])

            if key[end-1:end] == "ll"
                curr = "VV"
            elseif key[end-1:end] == "lc"
                curr = "VVc"
            end
            str = curr * " - set " * impr_set
            if str ∉ groupLabel
                push!(groupLabel, str)
                p0 = (i-1)*length(mykeys)*(length(FITCUT)+2) + (k-1)*(length(FITCUT)+2)
                pf = p0 + length(FITCUT)+1
                push!(groupPos, (p0+pf)/2)
            end
        end
        y += 1

        res  = myFactor * MA[impr_set]["average"]["res"][key]; uwerr(res)
        syst = myFactor * info[impr_set]["average"]["syst"][key]

        if !isnothing(sqrtt0_ph_TAR)
            der = mchist(res, "sqrtt0 [fm]")[1] / artificial_err
            res_t0shift = res + value(sqrtt0_ph_TAR - sqrtt0_ph) * der; uwerr(res_t0shift)
        else
            res_t0shift = res
        end

        ax1.errorbar(value(res), xerr=err(res), y, 0.0, fmt="o", mfc="none", color="green", ms=10, capsize=2, alpha=0.4)
        ax1.errorbar(value(res), xerr=sqrt(err(res)^2+syst^2), y, 0.0, fmt="o", mfc="none", color="green", ms=10, capsize=2, alpha=0.4)
        if !isnothing(sqrtt0_ph_TAR)
            ax1.errorbar(value(res_t0shift), xerr=err(res_t0shift), y+0.1, 0.0, fmt="o", mfc="none", color="green", ms=10, capsize=2, alpha=0.4)
            ax1.errorbar(value(res_t0shift), xerr=sqrt(err(res_t0shift)^2+syst^2), y+0.1, 0.0, fmt="o", mfc="none", color="green", ms=10, capsize=2, alpha=0.4)
        end
        push!(yTicksPos, y)
        push!(yTicks, latexstring("\\rm{Average}"))
        y += 1
    end
end
if !isnothing(Expec_RES)
    y += 1
    ax1.errorbar(Expec_RES[1], xerr=Expec_RES[2], y, 0.0, fmt="s", color="gray", ms=10, capsize=2, alpha=0.4)
    ax1.errorbar(Expec_RES[1], xerr=sqrt(Expec_RES[2]^2+Expec_RES[3]^2), y, 0.0, fmt="s", color="gray", ms=10, capsize=2, alpha=0.4)
end
ax1.axvspan(value(myRES)-err(myRES), value(myRES)+err(myRES), color="limegreen", alpha=0.3)
ax1.axvspan(value(myRES)-myERR, value(myRES)+myERR, color="limegreen", alpha=0.3, label="SD paper result")
diag_str = diag == "LO" ? "\\rm{LO}" : (diag == "NLOa&b" ? "\\rm{NLO}_{\\rm{a}\\&\\rm{b}}" : "\\rm{NLO}_{\\rm{$(diag[end])}}")
comp_str = comp[1] == '∆' ? (comp[end] != 'b' ? "\\Delta_{ls}(a_{\\mu})" : "\\Delta_{lc}(b)") : "\\mathrm{$(comp[2]),$(comp[3])}"
V_str    = VREF ? "(V_{\\rm{ref}})" : ""
Q_str    = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "($Q\\ GeV)" : ""
fact_str = diag == "LO" ? "\\times10^{10}" : "\\times10^{11}"
if wind == "NW"
    ax1.set_xlabel(latexstring("$(factor[comp*"s"])a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]$Q_str$fact_str"))
else
    ax1.set_xlabel(latexstring("$(factor[comp*"s"])\\left(a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]$(V_str)\\right)^{\\rm{$wind}}$Q_str$fact_str"))
end
ax1.set_yticks(yTicksPos, yTicks, rotation = 30, fontsize=15)
labels = ax1.get_yticklabels()
for label in labels
    if label.get_text() == "\$\\rm{Average}\$"
        label.set_color("green")  # Change the color
    end
end

ax2.set_xlim(0, 1)  # Assuming weights range from 0 to 1
ax2.set_ylim(ax1.set_ylim())
ax2.set_xlabel("Weight")
ax2.set_yticks([])  # Remove redundant y-ticks on the second plot

# Create an additional axis for brackets
axE = ax1.twiny()  # Create a twin x-axis
axE.set_xlim(-1, 0)  # Adjust limits for space
axE.set_xticks([])  # Remove x-ticks
axE.spines["top"].set_visible(false)  # Hide top border
axE.spines["bottom"].set_visible(false)  # Hide bottom border
axE.spines["right"].set_visible(false)  # Hide right border
axE.spines["left"].set_visible(false)  # Hide left border
FITCUT
for i in eachindex(groupLabel)
    xText = "beta&mass" in FITCUT ? -1.24 :  -1.14
    axE.text(xText, groupPos[i], groupLabel[i], ha="right", va="center", rotation = 90, fontsize=15)
end

tight_layout()
display(gcf())
if SAVE
    RESCstr = !RESC ? "" : "_resc"
    p = create_path(path_plot,[diag,wind,"$(diag)_$(comp)_StabilityPlot$(RESCstr).pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p)
end
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##


@info("Model average plot")

nBEST = 50

Norm2DOF = false

color_dict = Dict(
    ""     => Dict("cc" => "cyan"),
    "1"    => Dict("ll" => "blue",  "lc" => "red"  ),
    "2"    => Dict("ll" => "green", "lc" => "brown")
); color_dict["1old"] = color_dict["1"]

alpha_dict = Dict(
    "None" => 0.8,
    "beta" => 0.6,
    "mass" => 0.4,
    "beta&mass" => 0.2
)


myFactor = factor[comp] # * (diag == "LO" ? 1 : 10)

argBEST = sortperm(INFO["weight"],rev=true)[1:nBEST]
sort!(argBEST)

PVAL = sum(isnothing.(vcat([vcat([vcat([getfield.(fitres[impr_set][FitCut][key],:pval) for FitCut in FITCUT]...) for key in mykeys]...) for impr_set in IMPR_SET]...))) < 1

w_ = Dict(); p_ = Dict(); r_ = Dict(); c_ = Dict(); m_ = Dict(); x_ = Dict();
i = 1
for impr_set in IMPR_SET
    w_[impr_set] = Dict(); p_[impr_set] = Dict(); r_[impr_set] = Dict(); c_[impr_set] = Dict(); m_[impr_set] = Dict(); x_[impr_set] = Dict();
    for key in mykeys
        w_[impr_set][key] = Dict(); p_[impr_set][key] = Dict(); r_[impr_set][key] = Dict(); c_[impr_set][key] = Dict(); m_[impr_set][key] = Dict(); x_[impr_set][key] = Dict();
        for FitCut in FITCUT
            j = i + length(res[impr_set][FitCut][key]) - 1

            arg  = argBEST[(i .<= argBEST) .& (argBEST .<= j)]
            arg_ = arg .- (i - 1)

            x_[impr_set][key][FitCut] = i:j
            w_[impr_set][key][FitCut] = INFO["weight"][arg]
            p_[impr_set][key][FitCut] = getfield.(fitres[impr_set][FitCut][key],:pval)[arg_]
            r_[impr_set][key][FitCut] = res[impr_set][FitCut][key][arg_]
            c_[impr_set][key][FitCut] = !Norm2DOF ? (getfield.(fitres[impr_set][FitCut][key],:chi2)./getfield.(fitres[impr_set][FitCut][key],:chi2exp))[arg_] : (getfield.(fitres[impr_set][FitCut][key],:chi2)./getfield.(fitres[impr_set][FitCut][key],:dof))[arg_]
            m_[impr_set][key][FitCut] = modelinfo[FitCut]["label_tot_isov"][arg_]
            
            i = j + 1
        end
    end
end

v = myFactor*RES.mean
e = myFactor*RES.err; esyst = myFactor*sqrt(RES.err^2+INFO["syst"]^2)

fig = figure(figsize=(15,10))

subplots_adjust(hspace=0.1)
subplot(411)
ax1 = gca()
x = 1:nBEST

PyPlot.title("$nBEST best fits")
fill_between(x, (v-e) * 1e-10, (v+e) * 1e-10, color="limegreen", alpha=0.2)
fill_between(x, (v-esyst) * 1e-10, (v+esyst) * 1e-10, color="limegreen", alpha=0.2)

i = 1
for impr_set in IMPR_SET
    for key in mykeys
        for FitCut in FITCUT
            r = myFactor .* r_[impr_set][key][FitCut] * 1e-10; uwerr.(r)
            j = i + length(r) - 1
            errorbar(collect(i:j), value.(r), err.(r), fmt="o", mfc="none", color=color_dict[impr_set]["$(key[end-1:end])"], alpha=alpha_dict[FitCut], ms=10, capsize=2)
            i = j + 1
        end
    end
end
PyPlot.xticks([], [])
# ax1.yaxis.offsetText.set_visible(false)   # hide automatic offset
# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))  # force scientific
# ax1.yaxis.get_offset_text().set_text(L"\times10^{-11}")
# ax1.yaxis.offsetText.set_visible(true)
# setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels
diag_str = diag == "LO" ? "\\rm{LO}" : (diag == "NLOa&b" ? "\\rm{NLO}_{\\rm{a}\\&\\rm{b}}" : "\\rm{NLO}_{\\rm{$(diag[end])}}")
comp_str = comp[1] == '∆' ? (comp[end] != 'b' ? "\\Delta_{ls}(a_{\\mu})" : "\\Delta_{lc}(b)") : "\\mathrm{$(comp[2]),$(comp[3])}"
Q_str    = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "($Q)" : ""
# fact_str = diag == "LO" ? "\\times10^{10}" : "\\times10^{11}"
if wind == "NW"
    ylabel(latexstring("a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]$Q_str"))
else
    ylabel(latexstring("\\left(a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]\\right)^{\\rm{$wind}}$Q_str"))
end
# ylim([v-6*e,v+6*e])

subplot(412)
ax2=gca()

i = 1
for impr_set in IMPR_SET
    for key in mykeys
        for FitCut in FITCUT
            w = w_[impr_set][key][FitCut]
            j = i + length(w) - 1
            PyPlot.bar(collect(i:j), w, color=color_dict[impr_set]["$(key[end-1:end])"], alpha=alpha_dict[FitCut])
            i = j + 1
        end
    end
end
PyPlot.xticks([], [])
# setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels
ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))

subplot(413)
ax3=gca()

i = 1
for impr_set in IMPR_SET
    for key in mykeys
        for FitCut in FITCUT
            c = c_[impr_set][key][FitCut]
            j = i + length(c) - 1
            PyPlot.bar(collect(i:j), c, color=color_dict[impr_set]["$(key[end-1:end])"], alpha=alpha_dict[FitCut])
            i = j + 1
        end
    end
end
PyPlot.xticks([], [])
# setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels
!Norm2DOF ? ylabel(L"$\chi^2/\chi^2_{\mathrm{exp}}$") : ylabel(L"$\chi^2/\mathrm{d.o.f.}$")
# ylim([0.0,2.0])

if PVAL
    subplot(414)

    i = 1
    for impr_set in IMPR_SET
        for key in mykeys
            for FitCut in FITCUT
                p = p_[impr_set][key][FitCut]
                j = i + length(p) - 1
                PyPlot.bar(collect(i:j), p, color=color_dict[impr_set]["$(key[end-1:end])"], alpha=alpha_dict[FitCut])
                i = j + 1
            end
        end
    end
    ylabel(L"$\rm{p-values}$")
    ylim([0.0,1.0])
end

m = vcat([ vcat([vcat([m_[impr_set][key][FitCut] for FitCut in FITCUT]...) for key in mykeys]...) for impr_set in IMPR_SET]...)
m_str = []; x_str = []
i = 0
while i < length(m)
    i += 1
    str = m[i][2]
    for s in m[i][3:end]
        str *= ","*s
    end
    push!(m_str, str)
    if m[i][2:end] == m[i+1][2:end]
        push!(x_str, i+0.5)
        i += 1
    else
        push!(x_str, i)
    end
end
PyPlot.xticks(x_str, m_str, rotation=70, ha="right")
tight_layout()
display(fig)
close("all")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Hit plot")

SAVE     = true
OVERSAVE = false

path_bdio = path_bdio_dict["local"]


color_dict = Dict(
    ""     => Dict("cc" => "gold"),
    "1"    => Dict("ll" => "blue",  "lc" => "red"  ),
    "2"    => Dict("ll" => "green", "lc" => "brown")
); color_dict["1old"] = color_dict["1"]
fmt_dict = Dict(
    "None"      => "o",
    "mass"      => "s",
    "beta"      => "d",
    "beta&mass" => "x"
)

myFactor = factor[comp] * (diag == "LO" ? 1 : 10)

RESTOT = []; WEIGHT = []
for (i,impr_set) in enumerate(IMPR_SET)
    for (k,key) in enumerate(mykeys)
        for FitCut in FITCUT
            restot = myFactor .* MA[impr_set][FitCut]["res_tot"][key]; uwerr.(restot)
            weight = info[impr_set][FitCut]["weight"][key]

            # to make the plot go faster we only plot the best 1000 prediction
            # bitvec = weight .> maximum(weight)/100000

            # restot = restot[bitvec]; weight = weight[bitvec]

            color = color_dict[impr_set][key[end-1:end]]
            fmt   = fmt_dict[FitCut]

            errorbar(value.(restot), xerr=err.(restot), weight, 0.0, fmt=fmt, color=color, ms=5, capsize=1, alpha=0.2, label="impr. set $impr_set; discr. $(key[end-1:end])")
            push!(RESTOT,restot)
            push!(WEIGHT,weight)
        end
    end
end

x_lim = xlim()
y_lim = ylim()

RESTOT = vcat(RESTOT...)
WEIGHT = vcat(WEIGHT...)

valerr = [value.(RESTOT) err.(RESTOT)]
dist = [Normal(val,err) for (val,err) in [collect(row) for row in eachrow(valerr)]]

x = collect(range(x_lim[1],x_lim[2],length=1000))
Npdf = pdf.(dist, Ref(x))

y = zeros(length(x))
for (i,npdf) in enumerate(Npdf)
    y += npdf.*WEIGHT[i]
end
y *= maximum(WEIGHT)/maximum(y)

fill_between(x, y, 0.0, color="orange", alpha=0.2)
PyPlot.plot(x, y, color="orange", ls="--", alpha=0.5)

fill_betweenx([0.0,1.0], myFactor*(RES.mean-RES.err), myFactor*(RES.mean+RES.err), color="gray", alpha=0.2)
fill_betweenx([0.0,1.0], myFactor*(RES.mean-sqrt(RES.err^2+INFO["syst"]^2)), myFactor*(RES.mean+sqrt(RES.err^2+INFO["syst"]^2)), color="gray", alpha=0.2)

diag_str = diag == "LO" ? "\\rm{LO}" : (diag == "NLOa&b" ? "\\rm{NLO}_{\\rm{a}\\&\\rm{b}}" : "\\rm{NLO}_{\\rm{$(diag[end])}}")
comp_str = comp[1] == '∆' ? (comp[end] != 'b' ? "\\Delta_{ls}(a_{\\mu})" : "\\Delta_{lc}(b)") : "\\mathrm{$(comp[2]),$(comp[3])}"
V_str    = VREF ? "(V_{\\rm{ref}})" : ""
Q_str    = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "($Q\\ GeV)" : ""
fact_str = diag == "LO" ? "\\times10^{10}" : "\\times10^{11}"
if wind == "NW"
    xlabel(latexstring("$(factor[comp*"s"])a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]$Q_str$fact_str"))
else
    xlabel(latexstring("$(factor[comp*"s"])\\left(a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]$(V_str)\\right)^{\\rm{$wind}}$Q_str$fact_str"))
end
# legend()
ylim(y_lim)
ylabel(latexstring("w"))
display(gcf())
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##


@info("SD (aµ + b) stability plot")

diag = "NLOb"
comp = "g33"  #  g33  gCCconn

StdDer = false
tlImpr = true

QRES  = "best"  # "average"  "best"
QLIST = Qlist  #  [5.0]  [5.0,8.0]  Qlist

path_bdio = path_bdio_dict["local"]


factor = diag == "LO" ? 1 : 10

if comp == "gCCconn"
    factor *= 4/9 
end

b33Pert = TXTread_bQ(path_bPert,diag); uwerr.(b33Pert)

amu0, info0 = BDIOread_MAtot(path_bdio,diag,"SD",comp,StdDer=StdDer,tlImpr=tlImpr)

RES  = [amu0]
SYST = [info0["syst"]]
for Q in QLIST
    amu, amuInfo = BDIOread_MAtot(path_bdio,diag,"SDsub",comp,StdDer=StdDer,tlImpr=tlImpr,Q=Q)
    if comp == "g33"
        bq = value(b33Pert[Q .== Qlist][1])
        syst = sqrt(amuInfo["syst"]^2 + Window("SD")(0)^2 * err(b33Pert[Q .== Qlist][1])^2)
    elseif comp == "gCCconn"
        ∆lc, ∆lcInfo = BDIOread_MAtot(path_bdio,diag,"SDsub","∆lc_b",StdDer=StdDer,Q=Q)
        bq = 2*value(b33Pert[Q .== Qlist][1]) + ∆lc
        syst = sqrt(amuInfo["syst"]^2 + Window("SD")(0)^2 * (4*err(b33Pert[Q .== Qlist][1])^2 + ∆lcInfo["syst"]^2))
    end
    push!(RES ,amu + Window("SD")(0) * bq)
    push!(SYST,syst)
end

RES  .*= factor; uwerr.(RES)
SYST .*= factor

ERR = sqrt.(err.(RES).^2 .+ SYST.^2)

fig = figure(figsize=(10,8))
# PyPlot.title("Convergence plot [diag. $diag]")
errorbar(union([0.0],1.0./QLIST), value.(RES), err.(RES), fmt="o", color="black", ms=5, capsize=2)
errorbar(union([0.0],1.0./QLIST), value.(RES), ERR, fmt="o", color="black", ms=5, capsize=2)

xmin, xmax = xlim()
if QRES == "average"
    weight = vcat([1 / ERR[1]^2,1 ./ ERR[2:end].^2 ./ length(ERR[2:end])]...); weight ./= sum(weight)  # 1 ./ (ERR.^2 .* sum(1 ./ ERR.^2))
    VAL = sum(weight .* value.(RES))
    UNC = sqrt(sum(weight .* ERR)^2 + abs(sum(weight .* value.(RES).^2) - sum(weight .* value.(RES))^2))
elseif QRES == "best"
    VAL = value(RES[argmin(ERR)])
    UNC = ERR[argmin(ERR)]
else
    VAL = value(RES[findfirst(QRES .== QLIST) + 1])
    UNC = ERR[findfirst(QRES .== QLIST) + 1]
    # fill_between([0,10],(value(RES[findfirst(QRES .== QLIST) + 1])-ERR[findfirst(QRES .== QLIST) + 1]),(value(RES[findfirst(QRES .== QLIST) + 1])+ERR[findfirst(QRES .== QLIST) + 1]), color="gray", alpha=0.4)
end
fill_between([0,10],VAL-UNC,VAL+UNC, color="gray", alpha=0.4)

xlabel(latexstring("1/Q\\ [\\rm{GeV}^{-1}]"))
diag_str = diag == "NLOa&b" ? "\\rm{NLOa}\\&\\rm{b}" : diag
comp_str = comp[1] == '∆' ? (comp[end] != 'b' ? "\\Delta_{ls}(a_{\\mu})" : "\\Delta_{lc}(b)") : "\\mathrm{$comp}"
fact_str = diag == "LO" ? "10" : "11"
if comp == "g33"
    ylabel(latexstring("\\left[\\left(a_{\\mu}^{(\\rm{3,3})}[\\rm{$diag}]\\right)^{\\rm{SD}}_{\\rm{sub}}(Q^2)+\\omega_{\\rm{SD}}(0)b^{\\rm{(3,3)}}[\\rm{$diag}](Q^2)\\right]\\times 10^{$(fact_str)}"))
elseif comp == "gCCconn"
    ylabel(latexstring("\\frac{4}{9}\\left[\\left(a_{\\mu}^{(\\rm{C,C})}[\\rm{$diag}]\\right))^{\\rm{SD}}_{\\rm{sub}}(Q^2)+\\omega_{\\rm{SD}}(0)b^{\\rm{(C,C)}}[\\rm{$diag}](Q^2)\\right]\\times 10^{$(fact_str)}"))
end
xlim(xmin,xmax)

ax2 = gca()[:twiny]()
ax2[:set_xlim](xmin, xmax)
ax2[:set_xticks](union([0.0],1.0./QLIST))
ax2[:set_xticklabels](union([latexstring("\\infty")], QLIST))
ax2[:set_xlabel](latexstring("Q\\ [\\rm{GeV}]"))

ax2[:xaxis][:set_label_coords](0.5, 1.12)  # Adjust the second value to move the label upwards

tight_layout()
display(gcf())
close()


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##


@info("Scale and Vref comparison")

diag = "NLOa"  # LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
wind = "ID"  # NW  SD  ID  LD  ILD
comp = "g33"  # g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b

Q = 5.0  # virtuality for SDsub

BLIND = false

STD_DERIV  = false
tl_IMPR    = false

path_bdio = path_bdio_dict["local"]


# PyPlot.title("Scale and Vref crosscheck")

x = 0.0
for RESC in [false,true]
    for VREF in [false,true]
        x+=1.
        AMU, INFO = BDIOread_MAtot(path_bdio,diag,wind,comp,resc=RESC,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q); uwerr(AMU)
        SYSTerr   = INFO["syst"]

        scale_ph = !RESC ? sqrtt0_ph_Regensburg :  fPi_ph_PDGFLAG
        SCALEerr = get_t0err([AMU],scale_ph,resc=RESC)[1]

        if VREF
            FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,wind)
            AMU += FVC_ChPT
        end
        uwerr(AMU)

        errorbar(x, AMU.mean, AMU.err, fmt="o", color="black", ms=5, capsize=2)
        errorbar(x, AMU.mean, sqrt(AMU.err^2+SYSTerr^2), fmt="o", color="black", ms=5, capsize=2)
        errorbar(x, AMU.mean, sqrt(AMU.err^2+SYSTerr^2+SCALEerr^2), fmt="o", color="black", ms=5, capsize=2)
    end
end
PyPlot.xticks([1,2,3,4], ["t0","t0 (Vref)","fPi","fPi (Vref)"], rotation = 30, fontsize=15)
display(gcf())
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

RES = 
[
    [-35.59972,0.09114,0.06019,0.04326],
    # [-35.58123,0.09663,0.05868,0.04395],
    [-35.65578,0.09550,0.07430,0.04395],
    # [-35.63630,0.10050,0.07760,0.04326],
    [-35.67449,0.09408,0.06824,0.04098],
    # [-35.65642,0.09909,0.07088,0.04156],
    [-35.74577,0.09734,0.08227,0.03355],
    # [-35.72728,0.09734,0.08227,0.03355],
    [-35.68598,0.10070,0.05890,0.04098],
    [-35.69013,0.13117,0.15917,0.00968],
]
myRES = [-35.65578,0.09550,0.07430,0.04395]

# AleRES = [-35.73,0.16]
AleRES = [-35.697,0.16]


for (y,res) in enumerate(RES)
    fmt = "o"
    mfc = "none"
    if y == 2
        color = "blue"
    elseif y == 6
        color = "gray"
        fmt = "s"
    else
        color = "black"
    end
    errorbar(res[1], xerr=res[2], y, 0.0, fmt=fmt, mfc=mfc, color=color, ms=10, capsize=2, alpha=0.4)
    errorbar(res[1], xerr=sqrt(res[2]^2+res[3]^2), y, 0.0, fmt=fmt, mfc=mfc, color=color, ms=10, capsize=2, alpha=0.4)
    errorbar(res[1], xerr=sqrt(res[2]^2+res[3]^2+res[4]^2), y, 0.0, fmt=fmt, mfc=mfc, color=color, ms=10, capsize=2, alpha=0.4)
end
errorbar(AleRES[1], xerr=AleRES[2], 7, 0.0, fmt="x", color="purple", ms=10, capsize=2, alpha=0.4)
fill_betweenx([0,8], myRES[1]-sqrt(myRES[2]^2+myRES[3]^2+myRES[4]^2), myRES[1]+sqrt(myRES[2]^2+myRES[3]^2+myRES[4]^2), color="gray", alpha=0.4)
xlabel(latexstring("\\left(a_{\\mu}^{3,3}[\\rm{NLOa\\&b}]\\right)^{\\rm{ID}}\\times10^{11}"))
# PyPlot.yticks([1,2,3,4,5,6,7,8,9,10], ["Vref, t0","Vinf, t0","(w. a2phi4) Vref, t0","(w. a2phi4) Vinf, t0","(w. D201) Vref, t0","(w. D201) Vinf, t0","(w. D201, a2phi4) Vref, t0","(w. D201, a2phi4) Vinf, t0","Vref, fpi","(w. a2phi4) Alessandro C."], rotation = 30, fontsize=15)
PyPlot.yticks([1,2,3,4,5,6,7], ["(no a2phi4) t0","t0","(w. D201, no a2phi4) t0","(w. D201) t0","(no b=3.34) t0","(no a2phi4) fpi","(no b=3.34) Alessandro C."], rotation = 30, fontsize=15)
# ylim()
tight_layout()
display(gcf())
close()