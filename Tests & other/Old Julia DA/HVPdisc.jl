using ADerrors
using HVPobs

using Plots
using PyPlot
using Colors

using BDIO
using DataStructures
using Statistics

using ProgressBars

using ALPHAio

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

# ensList = ["H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302"]
ensList = ["H102", "C101", "S400", "N203", "N200", "D200", "N302"]
myensList, badensList = ensCheck(ensList, path_HVP, path_rw, path_ms, path_fvc, showbad=true)

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

corrDisc = Dict{String, Dict{String, Vector{Corr}}}()

t0 = Dict{String, uwreal}()

TMRa = Dict{String, Vector{uwreal}}()
TMRb = Dict{String, Vector{uwreal}}()
# TMRc = Dict{String, Matrix{uwreal}}()

@info("==> Reading charm data; renor. and impr. <==")
for (i,ens) in enumerate(ensInfo)

    @info("Starting ensemble: $(ens.id)     ($i/$(length(ensInfo)))")
    @info("Reading corr. ...")
    corrDisc[ens.id] = corr_disc(path_HVP, ens, path_rw = path_rw, frw_bcwd = true, impr = IMPR, impr_set = IMPR_SET, discr = ["ll","lc","cc"], std = STD_DERIV)

    p = joinpath(path_bdio, ens.id, string(ens.id,"a","_amuObs_Set1.bdio"))
    t0[ens.id] = read_BDIO(p, "HVP", "t0")[1]

    @info("Computing TMR Kernels ...")

    sym_points = Int64(length(corrDisc[ens.id]["ll"][1].obs)/2+1)
    
    factor = hbarc * sqrt(t0[ens.id])/t0_ph
    TMRa[ens.id] = factor^2 .* Tildef4a((massmu/factor) .* collect(0:sym_points-1),path_coef)
    TMRb[ens.id] = factor^2 .* Tildef4b((massmu/factor) .* collect(0:sym_points-1),path_coef)
    # TMRc[ens.id] = factor^4 .* Tildef4c((massmu/factor) .* collect(0:sym_points-1),path_coef)

    @info("-----------------------")
end

##=============> amu[NLO] disconnected HVP <=============##

HVPdisc = Dict{String, Dict{String, Dict{String, Vector{uwreal}}}}()
HVPdisc_int = Dict{String, Dict{String, Dict{String, Vector{Integrand}}}}()

@info("==> amu[NLO] charmed HVP <==")
for (i,ens) in enumerate(ensInfo)

    @info("Starting ensemble: $(ens.id)     ($i/$(length(ensInfo)))")
    @info("Computing HVP ...")

    HVPall, Integrandall = amuHVPNLO("a", corrDisc[ens.id]["ll"], TMRa[ens.id], t0ens=t0[ens.id], pl=false, int=true)
    HVPalc, Integrandalc = amuHVPNLO("a", corrDisc[ens.id]["lc"], TMRa[ens.id], t0ens=t0[ens.id], pl=false, int=true)

    HVPbll, Integrandbll = amuHVPNLO("b", corrDisc[ens.id]["ll"], TMRb[ens.id], t0ens=t0[ens.id], pl=false, int=true)
    HVPblc, Integrandblc = amuHVPNLO("b", corrDisc[ens.id]["lc"], TMRb[ens.id], t0ens=t0[ens.id], pl=false, int=true)

    # HVPcll, Integrandcll = amuHVPNLO("c", corrDisc[ens.id]["ll"], TMRc[ens.id], t0ens=t0[ens.id], pl=false, int=true)
    # HVPclc, Integrandclc = amuHVPNLO("c", corrDisc[ens.id]["lc"], TMRc[ens.id], t0ens=t0[ens.id], pl=false, int=true)

    HVPdisc[ens.id]     = Dict("a" => Dict("ll" => HVPall, "lc" => HVPalc), "b" => Dict("ll" => HVPbll, "lc" => HVPblc))#, "c" => Dict("ll" => HVPcll, "lc" => HVPclc))
    HVPdisc_int[ens.id] = Dict("a" => Dict("ll" => Integrandall, "lc" => Integrandalc), "b" => Dict("ll" => Integrandbll, "lc" => Integrandblc))#, "c" => Dict("ll" => Integrandcll, "lc" => Integrandclc))

    @info("-----------------------")
end

## <<<<===================================>>>> #
## <<<<==========>>>> PLOTS <<<<==========>>>> #
## <<<<===================================>>>> #

##=============> Plotting of integrand + FV corr. <=============##

disc = "ll"

@info("==> Plotting of integrand for disc=$disc <==")
for (i,ens) in enumerate(ensInfo[1:2])
    @info("Starting ensemble: $(ens.id)     ($i/$(length(ensInfo)))")
    for diagram in ["a","b"]

        integrand = HVPdisc_int[ens.id][diagram][disc]
        t_a = collect(0:length(integrand[1].obs)-1)

        label=[]
        fmt = ["s","o","^","d"]
        color = ["red","green","blue","gray"]
        for (i,subintegrand) in enumerate(integrand)
            uwerr.(subintegrand.obs)
            # val = Vector{Float64}()
            # err = Vector{Float64}()
            # for subint in subintegrand.obs
            #     if abs(value(subint)<0.0001)
            #         push!(val,value(subint))
            #         push!(err,ADerrors.err(subint))
            #     else
            #         push!(val,0.0)
            #         push!(err,0.0)
            #     end
            # end
            mult = diagram == "a" ? -1 : 1
            errorbar(t_a,mult.*abs.(value.(subintegrand.obs)), ADerrors.err.(subintegrand.obs), fmt=fmt[i], mfc="none", label=subintegrand.gamma, color=color[i], capsize=2)
            # errorbar(t_a, value.(subintegrand.obs),ADerrors.err.(subintegrand.obs), fmt=fmt[i], mfc="none", label=subintegrand.gamma, color=color[i], capsize=2)
            push!(label,subintegrand.gamma)
        end
        yscale("log")
        PyPlot.title(ens.id*": Disconnected integrands for diagram (" * diagram * ")")
        xlabel("t/a")
        diagram == "a" ? (ylabel(L"$|G^{\rm{disc}}(t)|$ $\tilde{f}^{(a)}_4(m_\mu t)$")) : (ylabel(L"$G(t)$ $\tilde{f}^{(b)}_4(m_\mu t)$"))
        ylim(-0.001,0.001)
        legend(label, loc="lower left")
        axis("tight")
        display(gcf())      #display the figure
        close("all")
    end
    @info("-----------------------")
end



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


ens = EnsInfo("N101")

Gllbare, Glcbare = corr_light(path_HVP, ens, path_rw = path_rw, frw_bcwd = true, renorm = false, impr = false,                cons = true, std = true, split = true)    
Gllimpr, Glcimpr = corr_light(path_HVP, ens, path_rw = path_rw, frw_bcwd = true, renorm = false, impr = true, impr_set = "1", cons = true, std = true, split = true)    
Gllimre, Glcimre = corr_light(path_HVP, ens, path_rw = path_rw, frw_bcwd = true, renorm = true,  impr = true, impr_set = "1", cons = true, std = true, split = true)    


uwerr.(Gllbare[1].obs); uwerr.(Glcbare[1].obs); uwerr.(Gllimpr[1].obs); uwerr.(Glcimpr[1].obs); uwerr.(Gllimre[1].obs); uwerr.(Glcimre[1].obs);  

Gllbare[1].obs[2:10]
Gllimpr[1].obs[2:10]
Gllimre[1].obs[2:10]

Glcbare[1].obs[2:10]
Glcimpr[1].obs[2:10]
Glcimre[1].obs[2:10]

##

GCllbare, GClcbare = corr_charm(path_HVP, ens, path_rw = path_rw, frw_bcwd = true, renorm = false, impr = false, std = true,                 cons = true)    
GCllimpr, GClcimpr = corr_charm(path_HVP, ens, path_rw = path_rw, frw_bcwd = true, renorm = false, impr = true,  impr_set = "1", std = true, cons = true)    
GCllimre, GClcimre = corr_charm(path_HVP, ens, path_rw = path_rw, frw_bcwd = true, renorm = true,  impr = true,  impr_set = "1", std = true, cons = true)  

uwerr.(GCllbare.obs); uwerr.(GClcbare.obs); uwerr.(GCllimpr.obs); uwerr.(GClcimpr.obs); uwerr.(GCllimre.obs); uwerr.(GClcimre.obs);  

GCllbare.obs[2:10]
GCllimpr.obs[2:10]
GCllimre.obs[2:10]