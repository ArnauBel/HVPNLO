using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using Plots
using PyPlot
using Colors

using BDIO
using DataStructures
using Statistics

using ProgressBars

# Export correlator, fv corr. and other usefull functions from "data_management.jl"
# Export TMR functions from "KernelTMRNLO.jl"
# Export HVP and FV from "amuNLO.jl"
# Export read_BDIO from "IO_BDIO.jl"

include("data_management.jl")
export get_Z3, get_Z8, get_Z08, corr33, corr88_conn, corr08_conn, corrR, ∆GHP, ensCheck

include("TMRKernel_NLO.jl")
export Tildef4aInner, Tildef4a, Tildef4bInner, Tildef4b, Tildef4cInner, Tildef4c

include("amuNLO.jl")
export amuHVPNLO, amu∆G

include("isovModel.jl")

include("tools.jl")
export get_w_from_fitcat, model_average

include("const.jl")

# Plot specifications

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

##=============> Path definition and ensamble choice <=============##

julia_script_directory = @__DIR__

path_HVP  = joinpath(julia_script_directory, "..", "LatticeData", "HVP_data")
path_rw   = joinpath(julia_script_directory, "..", "LatticeData", "rwf_deflated")
path_ms   = joinpath(julia_script_directory, "..", "LatticeData", "ms_t0_dat")
path_fvc  = joinpath(julia_script_directory, "..", "LatticeData", "JKMPI")

path_bdio = joinpath(julia_script_directory, "..", "obsBDIO")

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302"]
myensList, badensList = ensCheck(ensList, path_HVP, path_rw, path_ms, path_fvc, showbad=true)

ensInfo = EnsInfo.(myensList)

isempty(badensList) ? @info("Enough information has been found for all ensembles") : @info("Not enough information has been found concerning ensables $(join(badensList, ", "))\n")

##=============> Set of diagrams to compute <=============##

#setDiagrams = ["a","b"]
setDiagrams = ["c"]

##=============> BDIO; HVP extraction <=============##

# Initialize the Dictionaries
t0 = Dict{String, uwreal}()
HVP_33ll = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_33ll[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_33lc = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_33lc[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_88ll = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_88ll[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_88lc = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_88lc[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_08ll = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_08ll[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_08lc = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_08lc[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP∆G    = Dict{String, Dict{String, Dict{String, Vector{uwreal}}}}()
[HVP∆G[ensid] = Dict{String, Dict}("a" => Dict{String, Vector{uwreal}}(), "b" => Dict{String, Vector{uwreal}}(), "c" => Dict{String, Vector{uwreal}}()) for ensid in getfield.(ensInfo,:id)]

HVP_Rll  = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_Rll[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_Rlc  = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_Rlc[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]

# Read BDIO
for ensid in getfield.(ensInfo,:id)
    p = joinpath(path_bdio, ensid, string(ensid,"a","_amuObs_Set1.bdio"))
    t0[ensid] = read_BDIO(p, "HVP", "t0")[1]

    for diagram in setDiagrams
        for IMPR_SET in ["1","2"]
            p = joinpath(path_bdio, ensid, string(ensid,diagram,"_amuObs_Set$(IMPR_SET).bdio"))
            HVP_33ll[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "33_ll")[1]
            HVP_33lc[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "33_lc")[1]
            HVP_88ll[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "88_ll")[1]
            HVP_88lc[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "88_lc")[1]
            HVP_08ll[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "08_ll")[1]
            HVP_08lc[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "08_lc")[1]
            HVP∆G[ensid][diagram][IMPR_SET]    = read_BDIO(p, "HVP", "FVC_HP")

            HVP_Rll[ensid][diagram][IMPR_SET]  = HVP_33ll[ensid][diagram][IMPR_SET] + HVP_88ll[ensid][diagram][IMPR_SET] + HVP_08ll[ensid][diagram][IMPR_SET]
            HVP_Rlc[ensid][diagram][IMPR_SET]  = HVP_33lc[ensid][diagram][IMPR_SET] + HVP_88lc[ensid][diagram][IMPR_SET] + HVP_08lc[ensid][diagram][IMPR_SET]
        end
    end
end

## <<<<==================================>>>> #
## <<<<==========>>>> FITS <<<<==========>>>> #
## <<<<==================================>>>> # 

##=============> Data point creation for all ens in ensInfo <=============##

phi2_ph = 8*(t0_ph*mpi_ph/hc)^2; uwerr(phi2_ph)
phi4_ph = 8*t0_ph^2*((mK_ph/hc)^2 + 0.5*(mpi_ph/hc)^2); uwerr(phi4_ph)     # for TrM=cst ensembles

ydata_ll   = Dict{String, Dict}("a" => Dict{String, Vector{uwreal}}(), "b" => Dict{String, Vector{uwreal}}(), "c" => Dict{String, Vector{uwreal}}())
ydata_lc   = Dict{String, Dict}("a" => Dict{String, Vector{uwreal}}(), "b" => Dict{String, Vector{uwreal}}(), "c" => Dict{String, Vector{uwreal}}())
fit_ll     = Dict{String, Dict}("a" => Dict{String, Vector{FitRes}}(), "b" => Dict{String, Vector{FitRes}}(), "c" => Dict{String, Vector{FitRes}}())
fit_lc     = Dict{String, Dict}("a" => Dict{String, Vector{FitRes}}(), "b" => Dict{String, Vector{FitRes}}(), "c" => Dict{String, Vector{FitRes}}())
res_val_ll = Dict{String, Dict}("a" => Dict{String, Vector{uwreal}}(), "b" => Dict{String, Vector{uwreal}}(), "c" => Dict{String, Vector{uwreal}}())
res_val_lc = Dict{String, Dict}("a" => Dict{String, Vector{uwreal}}(), "b" => Dict{String, Vector{uwreal}}(), "c" => Dict{String, Vector{uwreal}}())

xdata = hcat([[1 / (8*t0[ens.id]), 8*t0[ens.id]*m_ens[ens.id]["m_pi"]^2, 8*t0[ens.id]*(m_ens[ens.id]["m_K"]^2+0.5*m_ens[ens.id]["m_pi"]^2)] for ens in ensInfo]...)
xdata = [xdata[1,:] xdata[2,:] xdata[3,:]]

for diagram in setDiagrams
    for IMPR_SET in ["1","2"]
        ydata_ll[diagram][IMPR_SET] = [HVP_Rll[ens.id][diagram][IMPR_SET] + HVP∆G[ens.id][diagram][IMPR_SET][end] for ens in ensInfo]
        ydata_lc[diagram][IMPR_SET] = [HVP_Rlc[ens.id][diagram][IMPR_SET] + HVP∆G[ens.id][diagram][IMPR_SET][end] for ens in ensInfo]

        fit_ll[diagram][IMPR_SET]=[]; fit_lc[diagram][IMPR_SET]=[]; res_val_ll[diagram][IMPR_SET]=[]; res_val_lc[diagram][IMPR_SET]=[]; 
        for i in collect(1:length(f_tot_isov))
            push!(fit_ll[diagram][IMPR_SET], fit_routine(f_tot_isov[i], value.(xdata), ydata_ll[diagram][IMPR_SET], n_par_tot_isov[i]))
            push!(fit_lc[diagram][IMPR_SET], fit_routine(f_tot_isov[i], value.(xdata), ydata_lc[diagram][IMPR_SET], n_par_tot_isov[i]))
            push!(res_val_ll[diagram][IMPR_SET], fit_ll[diagram][IMPR_SET][i].param[1])
            push!(res_val_lc[diagram][IMPR_SET], fit_lc[diagram][IMPR_SET][i].param[1])
        end
        uwerr.(res_val_ll[diagram][IMPR_SET])
        uwerr.(res_val_lc[diagram][IMPR_SET])
    end
end

## <<<<===================================>>>> #
## <<<<==========>>>> PLOTS <<<<==========>>>> #
## <<<<===================================>>>> # 

##=============> Continuum projection Plot <=============##

# Projected points

xdata_proj = [xdata[:,1] fill(phi2_ph, length(xdata[:,1])) fill(phi4_ph, length(xdata[:,1]))]; uwerr.(xdata_proj)
ydata_proj_ll   = Dict{String, Dict}("a" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()), "b" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()), "c" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()))
ydata_proj_lc   = Dict{String, Dict}("a" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()), "b" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()), "c" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()))

@info("Projection of the data points")
for diagram in setDiagrams
    @info("Starting diagram $diagram")
    for IMPR_SET in ["1","2"]
        @info("Starting set $IMPR_SET")
        for i in ProgressBar(1:length(f_tot_isov))
            myll = f_tot_isov[i](xdata_proj, fit_ll[diagram][IMPR_SET][i].param); uwerr.(myll)
            mylc = f_tot_isov[i](xdata_proj, fit_lc[diagram][IMPR_SET][i].param); uwerr.(mylc)
            push!(ydata_proj_ll[diagram][IMPR_SET], myll)
            push!(ydata_proj_lc[diagram][IMPR_SET], mylc)
        end
    end
end

# Compute the points for the continuum projection

xarr = [Float64.(range(0.0, 1.4*maximum(value.(xdata[:,1])), length=100)) fill(phi2_ph, 100) fill(phi4_ph, 100)]
yarr_ll   = Dict{String, Dict}("a" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()), "b" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()), "c" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()))
yarr_lc   = Dict{String, Dict}("a" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()), "b" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()), "c" => Dict{String, Vector{Vector{uwreal}}}("1" => Vector{Vector{uwreal}}(), "2" => Vector{Vector{uwreal}}()))

@info("Computing points for continuum projection plot")
for diagram in setDiagrams
    @info("Starting diagram $diagram")
    for IMPR_SET in ["1","2"]
        @info("Starting set $IMPR_SET")
        for i in ProgressBar(1:length(f_tot_isov))
            myll = f_tot_isov[i](xarr, fit_ll[diagram][IMPR_SET][i].param); uwerr.(myll)
            mylc = f_tot_isov[i](xarr, fit_lc[diagram][IMPR_SET][i].param); uwerr.(mylc)
            push!(yarr_ll[diagram][IMPR_SET], myll)
            push!(yarr_lc[diagram][IMPR_SET], mylc)
        end
    end
end
 
##== Plot 1

for diagram in setDiagrams
    argDom1 = argmin(abs.(getfield.(fit_ll[diagram]["1"], :chi2)./getfield.(fit_ll[diagram]["1"], :chi2exp).-1))
    argDom2 = argmin(abs.(getfield.(fit_ll[diagram]["2"], :chi2)./getfield.(fit_ll[diagram]["2"], :chi2exp).-1))
    for i in collect(1:length(f_tot_isov))

        Set1_ll = yarr_ll[diagram]["1"][i] 
        Set1_lc = yarr_lc[diagram]["1"][i] 
        Set2_ll = yarr_ll[diagram]["2"][i] 
        Set2_lc = yarr_lc[diagram]["2"][i]

        alphaSet1_ll = (0.25*min(1,fit_ll[diagram]["1"][i].chi2exp/fit_ll[diagram]["1"][i].chi2))^2
        alphaSet1_lc = (0.25*min(1,fit_lc[diagram]["1"][i].chi2exp/fit_lc[diagram]["1"][i].chi2))^2
        alphaSet2_ll = (0.25*min(1,fit_ll[diagram]["2"][i].chi2exp/fit_ll[diagram]["2"][i].chi2))^2
        alphaSet2_lc = (0.25*min(1,fit_lc[diagram]["2"][i].chi2exp/fit_lc[diagram]["2"][i].chi2))^2

        if i != argDom1
            fill_between(xarr[:,1], value.(Set1_ll)-ADerrors.err.(Set1_ll), value.(Set1_ll)+ADerrors.err.(Set1_ll), alpha=alphaSet1_ll, color="blue")
            fill_between(xarr[:,1], value.(Set1_lc)-ADerrors.err.(Set1_lc), value.(Set1_lc)+ADerrors.err.(Set1_lc), alpha=alphaSet1_lc, color="purple")
        else
            fill_between(xarr[:,1], value.(Set1_ll)-ADerrors.err.(Set1_ll), value.(Set1_ll)+ADerrors.err.(Set1_ll), alpha=alphaSet1_ll, color="blue", label = "Set 1 local-local")
            fill_between(xarr[:,1], value.(Set1_lc)-ADerrors.err.(Set1_lc), value.(Set1_lc)+ADerrors.err.(Set1_lc), alpha=alphaSet1_lc, color="purple", label = "Set 1 local-conserved")
        end
        if i!= argDom2
            fill_between(xarr[:,1], value.(Set2_ll)-ADerrors.err.(Set2_ll), value.(Set2_ll)+ADerrors.err.(Set2_ll), alpha=alphaSet2_ll, color="red")
            fill_between(xarr[:,1], value.(Set2_lc)-ADerrors.err.(Set2_lc), value.(Set2_lc)+ADerrors.err.(Set2_lc), alpha=alphaSet2_lc, color="orange")
        else
            fill_between(xarr[:,1], value.(Set1_ll)-ADerrors.err.(Set1_ll), value.(Set1_ll)+ADerrors.err.(Set1_ll), alpha=alphaSet1_ll, color="blue", label = "Set 2 local-local")
            fill_between(xarr[:,1], value.(Set1_lc)-ADerrors.err.(Set1_lc), value.(Set1_lc)+ADerrors.err.(Set1_lc), alpha=alphaSet1_lc, color="purple", label = "Set 2 local-conserved")
        end
    end
    PyPlot.title("Projection to continuum extrapolation")
    axvline(ls="dashed", color="black", lw=0.2, alpha=0.7) 
    xlabel(L"$a^2/8t_0$")
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{NLO}_{$diagram}]"))
    diagram == "a" ? legend(["Set 1 local-local","Set 1 local-conserved","Set 2 local-local","Set 2 local-conserved"], loc="upper center") : legend(["Set 1 local-local","Set 1 local-conserved","Set 2 local-local","Set 2 local-conserved"], loc="lower center")
    tight_layout()
    display(gcf())
    close()
end

##== Plot 2

for diagram in setDiagrams
    argDom1ll = argmin(abs.(getfield.(fit_ll[diagram]["1"], :chi2)./getfield.(fit_ll[diagram]["1"], :chi2exp).-1))
    argDom1lc = argmin(abs.(getfield.(fit_lc[diagram]["1"], :chi2)./getfield.(fit_lc[diagram]["1"], :chi2exp).-1))
    argDom2ll = argmin(abs.(getfield.(fit_ll[diagram]["2"], :chi2)./getfield.(fit_ll[diagram]["2"], :chi2exp).-1))
    argDom2lc = argmin(abs.(getfield.(fit_lc[diagram]["2"], :chi2)./getfield.(fit_lc[diagram]["2"], :chi2exp).-1))

    Set1_ll = yarr_ll[diagram]["1"][argDom1ll] 
    Set1_lc = yarr_lc[diagram]["1"][argDom1lc] 
    Set2_ll = yarr_ll[diagram]["2"][argDom2ll] 
    Set2_lc = yarr_lc[diagram]["2"][argDom2lc]

    errorbar(value.(xdata_proj[:,1]), xerr=ADerrors.err.(xdata_proj[:,1]), value.(ydata_proj_ll[diagram]["1"][argDom1ll]), yerr=ADerrors.err.(ydata_proj_ll[diagram]["1"][argDom1ll]), fmt="o", capsize=2, color="blue", mfc="none", label = "Set 1 local-local")
    errorbar(value.(xdata_proj[:,1]), xerr=ADerrors.err.(xdata_proj[:,1]), value.(ydata_proj_lc[diagram]["1"][argDom1lc]), yerr=ADerrors.err.(ydata_proj_lc[diagram]["1"][argDom1lc]), fmt="s", capsize=2, color="purple", mfc="none", label = "Set 1 local-conserved")
    errorbar(value.(xdata_proj[:,1]), xerr=ADerrors.err.(xdata_proj[:,1]), value.(ydata_proj_ll[diagram]["2"][argDom2ll]), yerr=ADerrors.err.(ydata_proj_ll[diagram]["2"][argDom2ll]), fmt="d", capsize=2, color="red", mfc="none", label = "Set 2 local-local")
    errorbar(value.(xdata_proj[:,1]), xerr=ADerrors.err.(xdata_proj[:,1]), value.(ydata_proj_lc[diagram]["2"][argDom2lc]), yerr=ADerrors.err.(ydata_proj_lc[diagram]["2"][argDom2lc]), fmt="^", capsize=2, color="orange", mfc="none", label = "Set 2 local-conserved")

    fill_between(xarr[:,1], value.(Set1_ll)-ADerrors.err.(Set1_ll), value.(Set1_ll)+ADerrors.err.(Set1_ll), alpha=0.5, color="blue")
    fill_between(xarr[:,1], value.(Set1_lc)-ADerrors.err.(Set1_lc), value.(Set1_lc)+ADerrors.err.(Set1_lc), alpha=0.5, color="purple")
    fill_between(xarr[:,1], value.(Set2_ll)-ADerrors.err.(Set2_ll), value.(Set2_ll)+ADerrors.err.(Set2_ll), alpha=0.5, color="red")
    fill_between(xarr[:,1], value.(Set2_lc)-ADerrors.err.(Set2_lc), value.(Set2_lc)+ADerrors.err.(Set2_lc), alpha=0.5, color="orange")

    PyPlot.title("Projection to continuum extrapolation")
    axvline(ls="dashed", color="black", lw=0.2, alpha=0.7) 
    xlabel(L"$a^2/8t_0$")
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{NLO}_{$diagram}]"))
    diagram == "a" ? legend(["Set 1 local-local","Set 1 local-conserved","Set 2 local-local","Set 2 local-conserved"], loc="upper center") : legend(["Set 1 local-local","Set 1 local-conserved","Set 2 local-local","Set 2 local-conserved"], loc="best")
    tight_layout()
    display(gcf())
    close()
end

## <<<<===========================================>>>> #
## <<<<==========>>>> MODEL AVERAGE <<<<==========>>>> #
## <<<<===========================================>>>> # 

##=============> MA for local-local, local-conserved and sets of impr. 1 and 2 <=============##

fitcat_ll = Dict{String, Dict}("a" => Dict{String, FitCat}(), "b" => Dict{String, FitCat}(), "c" => Dict{String, FitCat}())
fitcat_lc = Dict{String, Dict}("a" => Dict{String, FitCat}(), "b" => Dict{String, FitCat}(), "c" => Dict{String, FitCat}())

for diagram in setDiagrams
    for IMPR_SET in ["1","2"]
        str_ll = "all_data_$(diagram)_set$(IMPR_SET)_ll"
        str_lc = "all_data_$(diagram)_set$(IMPR_SET)_lc"

        fitcat_ll[diagram][IMPR_SET] = FitCat(xdata, ydata_ll[diagram][IMPR_SET], str_ll) 
        fitcat_lc[diagram][IMPR_SET] = FitCat(xdata, ydata_lc[diagram][IMPR_SET], str_lc)
        for i in collect(1:length(fit_ll[diagram][IMPR_SET]))
            push!(fitcat_ll[diagram][IMPR_SET].fit, fit_ll[diagram][IMPR_SET][i])
            push!(fitcat_lc[diagram][IMPR_SET].fit, fit_lc[diagram][IMPR_SET][i]) 
        end
    end
end

##---------------------------------------------------------------------------------------------

#weight_ll = Dict{String, Dict}("a" => Dict{String, Vector{Float64}}(), "b" => Dict{String, Vector{Float64}}(), "c" => Dict{String, Vector{Float64}}())
#weight_lc = Dict{String, Dict}("a" => Dict{String, Vector{Float64}}(), "b" => Dict{String, Vector{Float64}}(), "c" => Dict{String, Vector{Float64}}())
#
#for diagram in setDiagrams
#    for IMPR_SET in ["1","2"]
#        weight_ll[diagram][IMPR_SET] = get_w_from_fitcat([fitcat_ll[diagram][IMPR_SET]], norm=false)
#        weight_lc[diagram][IMPR_SET] = get_w_from_fitcat([fitcat_lc[diagram][IMPR_SET]], norm=false)
#    end
#end
#
#resa1_ll, systa1_ll = model_average(res_val_ll["a"]["1"], weight_ll["a"]["1"]); uwerr(resa1_ll)
#resa1_lc, systa1_lc = model_average(res_val_lc["a"]["1"], weight_lc["a"]["1"]); uwerr(resa1_lc)
#resa2_ll, systa2_ll = model_average(res_val_ll["a"]["2"], weight_ll["a"]["2"]); uwerr(resa2_ll)
#resa2_lc, systa2_lc = model_average(res_val_lc["a"]["2"], weight_lc["a"]["2"]); uwerr(resa2_lc)
#
#resb1_ll, systb1_ll = model_average(res_val_ll["b"]["1"], weight_ll["b"]["1"]); uwerr(resb1_ll)
#resb1_lc, systb1_lc = model_average(res_val_lc["b"]["1"], weight_lc["b"]["1"]); uwerr(resb1_lc)
#resb2_ll, systb2_ll = model_average(res_val_ll["b"]["2"], weight_ll["b"]["2"]); uwerr(resb2_ll)
#resb2_lc, systb2_lc = model_average(res_val_lc["b"]["2"], weight_lc["b"]["2"]); uwerr(resb2_lc)
#
#resa1_ll, systa1_ll
#resa1_lc, systa1_lc
#resa2_ll, systa2_ll
#resa2_lc, systa2_lc
#
#resb1_ll, systb1_ll
#resb1_lc, systb1_lc
#resb2_ll, systb2_ll
#resb2_lc, systb2_lc

##---------------------------------------------------------------------------------------------

##=============> General MA for all data points at once <=============##

result_tot = Dict{String, Vector{uwreal}}("a" => Vector{uwreal}(), "b" => Vector{uwreal}(), "c" => Vector{uwreal}())
weight_tot = Dict{String, Vector{Float64}}("a" => Vector{Float64}(), "b" => Vector{Float64}(), "c" => Vector{Float64}())

for diagram in setDiagrams
    #S1ll = [subvec[1] for subvec in getfield.(fitcat_ll[diagram]["1"].fit, :param)]
    #S1lc = [subvec[1] for subvec in getfield.(fitcat_lc[diagram]["1"].fit, :param)]
    #S2ll = [subvec[1] for subvec in getfield.(fitcat_ll[diagram]["2"].fit, :param)]
    #S2lc = [subvec[1] for subvec in getfield.(fitcat_lc[diagram]["2"].fit, :param)]
    
    #result_tot[diagram] = vcat(S1ll,S1lc,S2ll,S2lc)

    result_tot[diagram] = vcat(res_val_ll[diagram]["1"],res_val_lc[diagram]["1"],res_val_ll[diagram]["2"],res_val_lc[diagram]["2"])

    weight_tot[diagram] = get_w_from_fitcat([fitcat_ll[diagram]["1"],fitcat_lc[diagram]["1"],fitcat_ll[diagram]["2"],fitcat_lc[diagram]["2"]], norm=false)
end

#resa, systa = model_average(result_tot["a"], weight_tot["a"]); uwerr(resa)
#resb, systb = model_average(result_tot["b"], weight_tot["b"]); uwerr(resb)
resc, systc = model_average(result_tot["c"], weight_tot["c"]); uwerr(resc)

#print("Diagram a: $(value(resa)) ± $(ADerrors.err(resa)) ± $systa = $(value(resa)) ±  $(sqrt(ADerrors.err(resa)^2 + systa^2))\n")
#print("Diagram b: $(value(resb)) ± $(ADerrors.err(resb)) ± $systb = $(value(resb)) ±  $(sqrt(ADerrors.err(resb)^2 + systb^2))\n")
print("Diagram c: $(value(resc)) ± $(ADerrors.err(resc)) ± $systc = $(value(resc)) ±  $(sqrt(ADerrors.err(resc)^2 + systc^2))\n")


## <<<<=============================================>>>> #
## <<<<==========>>>> SU(3) sym. ens. <<<<==========>>>> #
## <<<<=============================================>>>> # 

ensListSU3 = ["H101", "B450", "N202", "N300"]
ensInfoSU3 = EnsInfo.(ensListSU3)

ydataSU3_ll   = Dict{String, Dict}("a" => Dict{String, Vector{uwreal}}(), "b" => Dict{String, Vector{uwreal}}(), "c" => Dict{String, Vector{uwreal}}())
ydataSU3_lc   = Dict{String, Dict}("a" => Dict{String, Vector{uwreal}}(), "b" => Dict{String, Vector{uwreal}}(), "c" => Dict{String, Vector{uwreal}}())
fitSU3_ll     = Dict{String, Dict}("a" => Dict{String, Vector{FitRes}}(), "b" => Dict{String, Vector{FitRes}}(), "c" => Dict{String, Vector{FitRes}}())
fitSU3_lc     = Dict{String, Dict}("a" => Dict{String, Vector{FitRes}}(), "b" => Dict{String, Vector{FitRes}}(), "c" => Dict{String, Vector{FitRes}}())
resSU3_val_ll = Dict{String, Dict}("a" => Dict{String, Vector{uwreal}}(), "b" => Dict{String, Vector{uwreal}}(), "c" => Dict{String, Vector{uwreal}}())
resSU3_val_lc = Dict{String, Dict}("a" => Dict{String, Vector{uwreal}}(), "b" => Dict{String, Vector{uwreal}}(), "c" => Dict{String, Vector{uwreal}}())

resultSU3_tot = Dict{String, Vector{uwreal}}("a" => Vector{uwreal}(), "b" => Vector{uwreal}(), "c" => Vector{uwreal}())
weightSU3_tot = Dict{String, Vector{Float64}}("a" => Vector{Float64}(), "b" => Vector{Float64}(), "c" => Vector{Float64}())

fitcatSU3_ll  = Dict{String, Dict}("a" => Dict{String, FitCat}(), "b" => Dict{String, FitCat}(), "c" => Dict{String, FitCat}())
fitcatSU3_lc  = Dict{String, Dict}("a" => Dict{String, FitCat}(), "b" => Dict{String, FitCat}(), "c" => Dict{String, FitCat}())

xdataSU3 = hcat([[1 / (8*t0[ens.id]), 8*t0[ens.id]*m_ens[ens.id]["m_pi"]^2, 8*t0[ens.id]*(m_ens[ens.id]["m_K"]^2+0.5*m_ens[ens.id]["m_pi"]^2)] for ens in ensInfoSU3]...)
xdataSU3 = [xdataSU3[1,:] xdataSU3[2,:] xdataSU3[3,:]]

for diagram in setDiagrams
    for IMPR_SET in ["1","2"]
        #=====> Data point creation for all ens in ensInfoSU3 <=====#
        ydataSU3_ll[diagram][IMPR_SET] = [HVP_Rll[ens.id][diagram][IMPR_SET] + HVP∆G[ens.id][diagram][IMPR_SET][end] for ens in ensInfoSU3]
        ydataSU3_lc[diagram][IMPR_SET] = [HVP_Rlc[ens.id][diagram][IMPR_SET] + HVP∆G[ens.id][diagram][IMPR_SET][end] for ens in ensInfoSU3]

        fitSU3_ll[diagram][IMPR_SET]=[]; fitSU3_lc[diagram][IMPR_SET]=[]; resSU3_val_ll[diagram][IMPR_SET]=[]; resSU3_val_lc[diagram][IMPR_SET]=[]; 
        for i in collect(1:length(f_tot_isov))
            push!(fitSU3_ll[diagram][IMPR_SET], fit_routine(f_tot_isov[i], value.(xdataSU3), ydataSU3_ll[diagram][IMPR_SET], n_par_tot_isov[i]))
            push!(fitSU3_lc[diagram][IMPR_SET], fit_routine(f_tot_isov[i], value.(xdataSU3), ydataSU3_lc[diagram][IMPR_SET], n_par_tot_isov[i]))
            push!(resSU3_val_ll[diagram][IMPR_SET], fitSU3_ll[diagram][IMPR_SET][i].param[1])
            push!(resSU3_val_lc[diagram][IMPR_SET], fitSU3_lc[diagram][IMPR_SET][i].param[1])
        end
        uwerr.(resSU3_val_ll[diagram][IMPR_SET])
        uwerr.(resSU3_val_lc[diagram][IMPR_SET])

        #=====> General MA for all SU(3) sym. points <=====#
        str_ll = "SU3_$(diagram)_set$(IMPR_SET)_ll"
        str_lc = "SU3_$(diagram)_set$(IMPR_SET)_lc"

        fitcatSU3_ll[diagram][IMPR_SET] = FitCat(xdataSU3, ydataSU3_ll[diagram][IMPR_SET], str_ll) 
        fitcatSU3_lc[diagram][IMPR_SET] = FitCat(xdataSU3, ydataSU3_lc[diagram][IMPR_SET], str_lc)
        for i in collect(1:length(fitSU3_ll[diagram][IMPR_SET]))
            push!(fitcatSU3_ll[diagram][IMPR_SET].fit, fitSU3_ll[diagram][IMPR_SET][i])
            push!(fitcatSU3_lc[diagram][IMPR_SET].fit, fitSU3_lc[diagram][IMPR_SET][i]) 
        end

    end
    resultSU3_tot[diagram] = vcat(resSU3_val_ll[diagram]["1"],resSU3_val_lc[diagram]["1"],resSU3_val_ll[diagram]["2"],resSU3_val_lc[diagram]["2"])

    weightSU3_tot[diagram] = get_w_from_fitcat([fitcatSU3_ll[diagram]["1"],fitcatSU3_lc[diagram]["1"],fitcatSU3_ll[diagram]["2"],fitcatSU3_lc[diagram]["2"]], norm=false)
end

resaSU3, systaSU3 = model_average(resultSU3_tot["a"], weightSU3_tot["a"]); uwerr(resaSU3)
resbSU3, systbSU3 = model_average(resultSU3_tot["b"], weightSU3_tot["b"]); uwerr(resbSU3)
rescSU3, systcSU3 = model_average(resultSU3_tot["c"], weightSU3_tot["c"]); uwerr(rescSU3)

print("Diagram a: $(value(resaSU3)) ± $(ADerrors.err(resaSU3)) ± $systaSU3 = $(value(resaSU3)) ±  $(sqrt(ADerrors.err(resaSU3)^2 + systaSU3^2))\n")
print("Diagram b: $(value(resbSU3)) ± $(ADerrors.err(resbSU3)) ± $systbSU3 = $(value(resbSU3)) ±  $(sqrt(ADerrors.err(resbSU3)^2 + systbSU3^2))\n")
print("Diagram c: $(value(rescSU3)) ± $(ADerrors.err(rescSU3)) ± $systcSU3 = $(value(rescSU3)) ±  $(sqrt(ADerrors.err(rescSU3)^2 + systcSU3^2))\n")




## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> 4-fit method [LO] <==========================##

COMP = "33" # 33, 88, cc conn, cc disc, c8 disc

nens = COMP in ["cc disc","c8 disc"] ? length(ensInfo)-length(ensSU3) : length(ensInfo)   # number of data ensembles
mpi_cut = "<360"  #  all  <360  <300

if mpi_cut != "all"
    cut = parse(Int64,mpi_cut[2:4])
end


baseDAdirect = joinpath(path_bdio,"DA","4-fit")
!ispath(baseDAdirect) ? mkdir(baseDAdirect) : nothing

DAdirect = joinpath(baseDAdirect,"mpi[$mpi_cut]")
!ispath(DAdirect) ? mkdir(DAdirect) : nothing


mykeys = DictComptoKey[COMP]

phi2_ph = 8*(t0_ph*(mpi_ph*1e-3)/hbarc)^2
phi4_ph = 8*t0_ph^2*(((mK_ph*1e-3)/hbarc)^2 + 0.5*((mpi_ph*1e-3)/hbarc)^2)

@info(" Fitting for LO and component $COMP")

xdata = []
ydata = Dict{String, Dict}("1" => Dict{String,  Vector{uwreal}}(), "2" => Dict{String, Vector{uwreal}}())

i = 0; j = 0
for ens in ensInfo
    println("- Reading data ensemble: $(ens.id)")
    if (COMP ∉ ["cc disc","c8 disc"]) || (ens.kappa_l != ens.kappa_s)

        println("   - Reading t0...")

        t0 = uwreal(0.0)
        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_t0"),"r")
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            t0 = ALPHAdobs_read_next(fb, size=sz, keys=ks)["t0"][1]
        end

        factor = hbarc * sqrt(t0)/t0_ph  
        mpi = m_ens[ens.id]["m_pi"]*factor * 1e3

        if mpi_cut == "all" || mpi < cut
            i += 1

            println("   - Reading HVP...")

            HVP = Dict{String, Dict}("1" => Dict{String, uwreal}(), "2" => Dict{String, uwreal}())

            for impr_set in ["1","2"]
                fb = BDIO_open(joinpath(path_bdio,"HVP&FVC",ens.id,"HVP","$(ens.id)_HVPLO_set$(impr_set)"),"r")
                val = Dict{String, Dict{String, uwreal}}()
                while ALPHAdobs_next_p(fb)
                    d = ALPHAdobs_read_parameters(fb)
                    nobs = d["nobs"]
                    dims = d["dimensions"]
                    ks = collect(d["keys"])
                    val["HVP"] =  ALPHAdobs_read_next(fb, keys=ks)
                end
                BDIO_close!(fb)
                info = load(joinpath(path_bdio,"HVP&FVC",ens.id,"HVP","$(ens.id)_HVPLOinfo_set$(impr_set).jld2"), "HVPinfo")
                HVPdict = merge(val,info)
                # apply syst.!
                if COMP != "cc conn"
                    uwreal_syst = Dict{String, uwreal}(); [uwreal_syst[key] = uwreal([0.0, HVPdict["HVPsyst"][key]], "syst from BM/cut-off") for key in mykeys]
                else
                    uwreal_syst = Dict{String, uwreal}(); uwreal_syst["gcc_ll_conn"] = uwreal_syst["gcc_lc_conn"] = uwreal([0.0,0.0], "syst from BM/cut-off")
                end
                [HVP[impr_set][key] = HVPdict["HVP"][key] + uwreal_syst[key] for key in mykeys]
            end

            println("   - Reading FVC...")

            FVC = Dict{String, Any}()

            fb_LO = BDIO_open(joinpath(path_bdio,"HVP&FVC",ens.id,"FVC","$(ens.id)_FVC_LO"),"r")
            val_LO = Dict{String, uwreal}()
            while ALPHAdobs_next_p(fb_LO)
                d = ALPHAdobs_read_parameters(fb_LO)
                sz = tuple(d["size"]...)
                ks = collect(d["keys"])
                val_LO = ALPHAdobs_read_next(fb_LO, size=sz, keys=ks)
            end
            BDIO_close!(fb_LO)
            [FVC[FVCtype] = val_LO[FVCtype][end] for FVCtype in ["FVCPi","FVCK"]]

            if ens.kappa_l == ens.kappa_s
                if COMP in ["33","88"]; myFVC = 1.5*FVC["FVCPi"]; else myFVC = 0.0; end
            else
                if COMP == "33"
                    myFVC = FVC["FVCPi"] + FVC["FVCK"]
                elseif COMP == "88"
                    myFVC = (2/9)*FVC["FVCK"] # (2/3).*FVC["FVCK"]
                else
                    myFVC = 0.0
                end
            end

            println("   - Creating 'x' data points...")

            push!(xdata, [1 / (8*t0), 8*t0*m_ens[ens.id]["m_pi"]^2, 8*t0*(m_ens[ens.id]["m_K"]^2+0.5*m_ens[ens.id]["m_pi"]^2)])

            println("   - Creating 'y' data points...")

            for impr_set in ["1","2"]
                for key in mykeys
                    i == 1 ? ydata[impr_set][key] = Vector{uwreal}() : nothing
                    push!(ydata[impr_set][key], HVP[impr_set][key] + myFVC)
                end
            end
        else
            println("   - The pion mass cut ($mpi_cut) has excluded this ensemble")
            j += 1
        end
    else
        println("   - The contributions from $COMP are 0 for this ensemble !!")
    end
end

nens -= j

@info(" Fitting points for LO (data ens = $nens)")

# isov_model1(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+  p[4] .* (log.(x[:,2]) .- log.(value.(phi2_ph)))
# isov_model2(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+  p[4] .* (x[:,2].^2 .- value.(phi2_ph).^2)
# isov_model3(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+  p[4] .* (x[:,2].^(-1) .- value.(phi2_ph).^(-1))
# isov_model4(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+  p[4] .* (x[:,2] .* log.(x[:,2]) .- value.(phi2_ph) .* log.(value.(phi2_ph)))

# f_tot_isov     = [isov_model1,isov_model2,isov_model3,isov_model4]
# label_tot_isov = ["logphi2","phi2sqr","phi2inv","phi2log"]
# n_par_tot_isov = [4,4,4,4]

isov_model1(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+ p[4] .* x[:,1].^(3/2) .+  p[5] .* (log.(x[:,2]) .- log.(value.(phi2_ph)))
isov_model2(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+ p[4] .* x[:,1].^(3/2) .+  p[5] .* (x[:,2].^2 .- value.(phi2_ph).^2)
isov_model3(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+ p[4] .* x[:,1].^(3/2) .+  p[5] .* (x[:,2].^(-1) .- value.(phi2_ph).^(-1))
isov_model4(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+ p[4] .* x[:,1].^(3/2) .+  p[5] .* (x[:,2] .* log.(x[:,2]) .- value.(phi2_ph) .* log.(value.(phi2_ph)))

f_tot_isov     = [isov_model1,isov_model2,isov_model3,isov_model4]
label_tot_isov = ["logphi2","phi2sqr","phi2inv","phi2log"]
n_par_tot_isov = [5,5,5,5]

xdata = hcat(xdata...)
xdata = [xdata[1,:] xdata[2,:] xdata[3,:]]

fit = Dict{String, Dict}(); par = Dict{String,Dict}()

fit = Dict{String, Dict}("1" => Dict{String, Vector{FitRes}}(), "2" => Dict{String, Vector{FitRes}}())
par = Dict{String, Dict}("1" => Dict{String, Vector{Vector{uwreal}}}(), "2" => Dict{String, Vector{Vector{uwreal}}}())
for impr_set in ["1","2"]
    println("      - Starting set "*impr_set)
    for key in mykeys
        println("         - Fitting for comp. $key ...")
        fit[impr_set][key] = Vector{FitRes}(); par[impr_set][key] = Vector{uwreal}()
        for i in ProgressBar(collect(1:length(f_tot_isov)))
            @suppress begin
                myfit = fit_routine(f_tot_isov[i], value.(xdata), ydata[impr_set][key], n_par_tot_isov[i], pval=true)
                push!(fit[impr_set][key], myfit)
                push!(par[impr_set][key], myfit.param)
            end
        end
    end
end


println("- Printing BDIO & JDL2...")

!ispath(joinpath(DAdirect,"Fit")) ? mkdir(joinpath(DAdirect,"Fit")) : nothing

diagpath = joinpath(DAdirect,"Fit","LO")
!ispath(diagpath) ? mkdir(diagpath) : nothing

pcomp = joinpath(diagpath,COMP)
!ispath(pcomp) ? mkdir(pcomp) : nothing

pFitRes = joinpath(pcomp,"FitRes.jld2")
save(pFitRes,"FitRes",fit)

io = IOBuffer()
write(io, "parameters")

fb = ALPHAdobs_create(joinpath(pcomp,"param"), io) 

for i in collect(1:length(par["1"][mykeys[1]]))
    for impr_set in ["1","2"]
        parDict = Dict{String,Array{uwreal}}()
        for key in mykeys
            parDict["diagLO_$(key)_set$(impr_set):[$i]"] = par[impr_set][key][i]
        end
        extra = Dict{String, Any}("Diag" => "LO", "Set" => impr_set)
        ALPHAdobs_write(fb, parDict, extra=extra)
    end
end
ALPHAdobs_close(fb)

io = IOBuffer()
write(io, "xydata")

fb = ALPHAdobs_create(joinpath(pcomp,"xydata"), io)

xDict = Dict{String,Array{uwreal}}("xdata" => xdata)
ALPHAdobs_write(fb, xDict)
for impr_set in ["1","2"]
    yDict = Dict{String,Array{uwreal}}()
    for key in mykeys
        yDict["$(key)_set$(impr_set)"] = ydata[impr_set][key]
    end
    extra = Dict{String, Any}("Set" => impr_set)
    ALPHAdobs_write(fb, yDict, extra=extra)
end
ALPHAdobs_close(fb)


println("- Printing Model information...")

infoDict = Dict{String,Any}(
    "length" => length(f_tot_isov),
    "nens" => nens,
    "n_par_tot_isov" => n_par_tot_isov,
    "label_tot_isov" => label_tot_isov,
)

pinfo = joinpath(pcomp,"ModelInfo.jld2")
save(pinfo,"info",infoDict)