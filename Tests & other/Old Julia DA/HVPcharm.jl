using ADerrors
using HVPobs

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

include("data_management.jl")              #already includes "const.jl"
export get_Z3, get_Z8, get_Z08, corr33, corr88_conn, corr08_conn, corrR, ∆GHP, ensCheck

include("TMRKernel_NLO.jl")
export Tildef4aInner, Tildef4a, Tildef4bInner, Tildef4b, Tildef4cInner, Tildef4c

include("amuNLO.jl")
export amuHVPNLO, amu∆G

include("IO_BDIO.jl")
export read_BDIO

# include("const.jl")

# Plot specifications

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

##=============> Path definition, ensamble choice and set of  <=============##

julia_script_directory = @__DIR__

path_HVP  = joinpath(julia_script_directory, "..", "LatticeData", "HVP_data")
path_rw   = joinpath(julia_script_directory, "..", "LatticeData", "rwf_deflated")
path_ms   = joinpath(julia_script_directory, "..", "LatticeData", "ms_t0_dat")
path_fvc  = joinpath(julia_script_directory, "..", "LatticeData", "JKMPI")

path_coef = joinpath(julia_script_directory, "..", "Coefficients")

path_bdio = joinpath(julia_script_directory, "..", "obsBDIO")

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302"]
myensList, badensList = ensCheck(ensList, path_HVP, path_rw, path_ms, path_fvc, showbad=true)
#subensList1 = ["H101", "B450", "N202", "N300"]
#subensList2 = ["H102", "N101", "C101", "S400"]
#subensList3 = ["N203", "N200", "D200", "N302"]
#myensList, badensList = ensCheck(subensList3, path_HVP, path_rw, path_ms, path_fvc, showbad=true)

ensInfo = EnsInfo.(myensList)

isempty(badensList) ? @info("Enough information has been found for all ensembles") : @info("Not enough information has been found concerning ensables $(join(badensList, ", "))\n")

## <<<<=============================================>>>> #
## <<<<==========>>>> DATA EXTRACTION <<<<==========>>>> #
## <<<<=============================================>>>> #

IMPR      = true
IMPR_SET  = "1"         # either "1" or "2"
RENORM    = true
STD_DERIV = false

@info("Impr. set: $IMPR_SET")

##=============> Reading, renormalization and impr. <=============##

corrCharm = Dict{String, Dict{String, Corr}}()

t0 = Dict{String, uwreal}()

TMRa = Dict{String, Vector{uwreal}}()
TMRb = Dict{String, Vector{uwreal}}()
# TMRc = Dict{String, Matrix{uwreal}}()

@info("==> Reading charm data; renor. and impr. <==")
for (i,ens) in enumerate(ensInfo)

    @info("Starting ensemble: $(ens.id)     ($i/$(length(ensInfo)))")
    @info("Reading corr. ...")
    Gll, Glc = corr_charm(path_HVP, ens, path_rw = path_rw, frw_bcwd = true, impr = IMPR, impr_set = IMPR_SET, cons = true, std = STD_DERIV)    
    corrCharm[ens.id] = Dict("ll" => Gll, "lc" => Glc)

    p = joinpath(path_bdio, ens.id, string(ens.id,"a","_amuObs_Set1.bdio"))
    t0[ens.id] = read_BDIO(p, "HVP", "t0")[1]

    @info("Computing TMR Kernels ...")

    sym_points = Int64(length(Gll.obs)/2+1)

    factor = hbarc * sqrt(t0[ens.id])/t0_ph  
    TMRa[ens.id] = factor^2 .* Tildef4a((massmu/factor) .* collect(0:sym_points-1),path_coef)
    TMRb[ens.id] = factor^2 .* Tildef4b((massmu/factor) .* collect(0:sym_points-1),path_coef)
    # TMRc[ens.id] = factor^4 .* Tildef4c((massmu/factor) .* collect(0:sym_points-1),path_coef)

    @info("-----------------------")
end

##=============> amu[NLO] charmed HVP <=============##

HVPcharm = Dict{String, Dict{String, Dict{String, uwreal}}}()
HVPcharm_int = Dict{String, Dict{String, Dict{String, Integrand}}}()

@info("==> amu[NLO] charmed HVP <==")
for (i,ens) in enumerate(ensInfo)

    @info("Starting ensemble: $(ens.id)     ($i/$(length(ensInfo)))")
    @info("Computing HVP ...")

    HVPall, Integrandall = amuHVPNLO("a", corrCharm[ens.id]["ll"], TMRa[ens.id], t0ens=t0[ens.id], pl=false, int=true)
    HVPalc, Integrandalc = amuHVPNLO("a", corrCharm[ens.id]["lc"], TMRa[ens.id], t0ens=t0[ens.id], pl=false, int=true)

    HVPbll, Integrandbll = amuHVPNLO("b", corrCharm[ens.id]["ll"], TMRb[ens.id], t0ens=t0[ens.id], pl=false, int=true)
    HVPblc, Integrandblc = amuHVPNLO("b", corrCharm[ens.id]["lc"], TMRb[ens.id], t0ens=t0[ens.id], pl=false, int=true)

    # HVPcll, Integrandcll = amuHVPNLO("c", corrCharm[ens.id]["ll"], TMRc[ens.id], t0ens=t0[ens.id], pl=false, int=true)
    # HVPclc, Integrandclc = amuHVPNLO("c", corrCharm[ens.id]["lc"], TMRc[ens.id], t0ens=t0[ens.id], pl=false, int=true)

    HVPcharm[ens.id]     = Dict("a" => Dict("ll" => HVPall[1], "lc" => HVPalc[1]), "b" => Dict("ll" => HVPbll[1], "lc" => HVPblc[1]))#, "c" => Dict("ll" => HVPcll[1], "lc" => HVPclc[1]))
    HVPcharm_int[ens.id] = Dict("a" => Dict("ll" => Integrandall[1], "lc" => Integrandalc[1]), "b" => Dict("ll" => Integrandbll[1], "lc" => Integrandblc[1]))#, "c" => Dict("ll" => Integrandcll[1], "lc" => Integrandclc[1]))

    @info("-----------------------")
end

## <<<<===================================>>>> #
## <<<<==========>>>> PLOTS <<<<==========>>>> #
## <<<<===================================>>>> #

##=============> Plotting of integrand + FV corr. <=============##

gamma = "ll"

@info("==> Plotting of integrand + FV corr. for gamma=$gamma <==")
for (i,ens) in enumerate(ensInfo)
    @info("Starting ensemble: $(ens.id)     ($i/$(length(ensInfo)))")
    for diagram in ["a","b"]

        integrand = HVPcharm_int[ens.id][diagram][gamma].obs
        uwerr.(integrand)
        t_a = collect(0:length(integrand)-1)

        label = "charm"
        errorbar(t_a, value.(integrand),ADerrors.err.(integrand), fmt="^", mfc="none", label=label, color="green", capsize=2)

        PyPlot.title(ens.id*": Charm integrand for diagram (" * diagram * ")")
        xlabel("t/a")
        diagram == "a" ? (ylabel(L"$G(t)$ $\tilde{f}^{(a)}_4(m_\mu t)$")) : (ylabel(L"$G(t)$ $\tilde{f}^{(b)}_4(m_\mu t)$"))
        #diagram == "a" ? (legend(label,loc="lower right")) : (legend(label,loc="upper right"))
        axis("tight")
        display(gcf())      #display the figure
        close("all")
    end
    @info("-----------------------")
end

