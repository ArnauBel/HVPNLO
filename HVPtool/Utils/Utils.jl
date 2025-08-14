module Utils

using HVPobs
using ADerrors
import ADerrors: err

using Statistics, Random

using HDF5
using PyPlot


# include("../Const/Const.jl")

include("../ModAver/ModAver.jl")
using .ModAver

include("../TMR/TMR.jl")
using .TMR

include("../uwConst.jl")


##-- ensCheck function

@doc raw"""
ensCheck(ens::String, path_HVP::String, path_rw::String, path_ms::String, path_fvc::String)

This function is called to check if the data paths fulfill all the necessary requirements for a given ensamble or a list of ensambles to compute an estimation of the HVP.

If a single ensamble is given, the function will return either "true" or "false" depending on if the ensamble fulfills the requirements or not. 
If a list of ensambles is given, the function will return a reduced list with only those ensambles which fulfill the requirements. Finally,  the 'showbad' key can be set to true for the function to also return a reduced list of the ensambles which do not fulfill the requirements.

Examples:
```@example
goodens, badens = ensCheck(["H101", "B450", "N202", "N300", "H102"], path_HVP, path_rw, path_ms, path_fvc, showbad=true)
```
"""
function ensCheck(ens::EnsInfo, ensidNOcharm_List::Vector{String},  ensidNOdisc_List::Vector{String}, path_HVP::String, path_rw::String, path_ms::String, path_fvc::String; data_status::Bool=false)
    ensid = ens.id

    light_req    = isdir(joinpath(path_HVP,ensid,"light"))
    strange_req  = isdir(joinpath(path_HVP,ensid,"light"))
    charm_req    = ensid ∉ ensidNOcharm_List ? isdir(joinpath(path_HVP,ensid,"light")) : true
    disc_req     = (ens.kappa_l != ens.kappa_s && ensid ∉ ensidNOdisc_List) ? isdir(joinpath(path_HVP,"disc",ensid)) : true
    rw_req       = !isempty(filter(x->occursin(ensid, x), readdir(path_rw, join=true))) || !isempty(filter(x->occursin(ensid, x), readdir(joinpath(path_rw,"reweight_deflated"), join=true)))
    ms_req       = !isempty(filter(x->occursin(ensid, x), readdir(path_ms, join=true))) 
    fvc_req      = !isempty(filter(x-> occursin("corr_blat_gsl", x) && occursin(ensid, x), readdir(path_fvc, join=true)))

    HVP_req = ([light_req, strange_req, charm_req] == [true, true, true])
    myBool  = ([HVP_req, disc_req, rw_req, ms_req, fvc_req] == [true, true, true, true, true])

    if data_status && !myBool
        str = "
Data status for ens $ensid: \n
- HVP data: ................. $HVP_req \n
   - light: ................. $light_req \n
   - strange: ............... $strange_req \n"
        if ensid ∉ ensidNOcharm_List
            str *= "
   - charm: ................. $charm_req \n "
        end
        if ens.kappa_l != ens.kappa_s && ensid ∉ ensidNOdisc_List
            str *= "
- Disc data: ................ $disc_req \n "
        end
        str *= "
- Reweighting data: ......... $rw_req \n 
- t0 data: .................. $ms_req \n 
- FVC data: ................. $fvc_req \n "
        println(str)
    end
    return myBool
end
function ensCheck(ensList::Vector{EnsInfo}, ensidNOcharm_List::Vector{String}, ensidNOdisc_List::Vector{String}, path_HVP::String, path_rw::String, path_ms::String, path_fvc::String; data_status::Bool=false)

    mask = [ensCheck(ens, ensidNOcharm_List, ensidNOdisc_List, path_HVP, path_rw, path_ms, path_fvc, data_status = data_status) for ens in ensList]
    
    return ensList[mask], ensList[.!mask]
end
export ensCheck

##-- 2D window structure

@doc raw"""
2D window structure for diagram NLOc
"""
struct Window2D
    func::Function
    function Window2D(str::String)
        delta = 0.15

        if str == "SD"
            d = 0.4
            @. func2Dsd(t,tau) =  1 - 0.5 * (1 + tanh((sqrt(t^2+tau^2)-d)/delta))
            return new(funcsd)
        elseif str == "ID"
            d1 = 0.4
            d2 = 1.0
            @. func2Did(t,tau) = 0.5 * (1 + tanh((sqrt(t^2+tau^2)-d1)/delta)) - 0.5 * (1 + tanh((sqrt(t^2+tau^2)-d2)/delta))
            return new(funcid)
        elseif  str == "LD"
            d = 1.0
            @. func2Dld(t,tau) = 0.5 * (1 + tanh((sqrt(t^2+tau^2)-d)/delta))
            return new(funcld)
        elseif str == "ILD" # intermediate and long distance
            d = 0.4
            @. func2Dild(t,tau) = 0.5 * (1 + tanh((sqrt(t^2+tau^2)-d)/delta))
            return new(funcild)
        else
            error("Window $(str) not defined.")
        end
    end
end
function (a::Window2D)(t,tau)
    return a.func(t,tau)
end
export Window2D


include("mesonMass.jl")
export meff_MA

include("3limpr.jl")
export treelevel_continuum_correlator, compute_HVPtl0

include("apply_corr.jl")
export apply_syst_HVP, apply_syst_HVP!, apply_syst_FVC, apply_syst_FVC!, HVP_VolCorrect, HVP_VolCorrect!, HVP_3limpr!

include("BoundMethod.jl")
export Eeff, corr_bound, bounding_method

include("spec_rec.jl")
export findfirst_uninterrupted, get_spectr_data, corr_n, reconstr_corr

include("resampling.jl")
export jackknife_resampling, bootstrap_resampling, jackknife_err

include("error_control.jl")
export set_fluc_to_zero!, add_t0_err!, get_t0err


end