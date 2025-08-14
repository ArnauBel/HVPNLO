using PyPlot
using HDF5
using Statistics
using Random

include("Const.jl")

include("Fit&MA.jl")

include("TMRKernel.jl")




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


##-- mass plateau searching function with model average, constant and 2-state fits

wpmm = Dict{String, Vector{Float64}}()
# wpmm["H101"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["H102r002"] = [5.0, -2.0, -1.0, -1.0]
# wpmm["H400"]     = [5.0, -1.5, -1.0, -1.0]
# wpmm["N202"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["N200"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["N203"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["N300"]     = [5.0, -1.5, -1.0, -1.0]
# wpmm["J303"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["J304"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["F300"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["J306"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["J307"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J500"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["E300"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["A654"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["E300"]     = [-1.0, 2.0, -1.0, -1.0]


function meff_MA(corr::Corr; 
    pl2state0::Vector{Float64}=[0.5,0.6], plconst0::Vector{Float64}=[0.7,1.0], plf::Float64=1.0, plstep::Int64=1, 
    mdof::Int64=4, 
    state_fit::Bool=true, 
    AIC::Bool=true, 
    returnfitMA::Bool=false, 
    plot::Bool=false, 
    pval::Bool=false,
    fitinfo::Bool=false
    )
    
    obs = corr.obs
    m_obs = meff(obs[1:Int64(length(obs)/2+1)]) # first and last data points are lost (-2) and derivative is taken (-1); length(meff_) = length(obs) - 3
    len = length(m_obs)

    @. const_model(x,p) = p[1] + 0*x

    plconst0_vec = collect(floor(Int64,plconst0[1]*len):plstep:ceil(Int64,min(plf*len-mdof,plconst0[2]*len)))
    pl_f = ceil(Int64,plf*len)

    if isempty(plconst0_vec)
        error("Not enough data to fit with d.o.f. ≥ $mdof. Please decrease $mdof or make increase the plateau search space. ")
    end

    wpm = corr.id in keys(wpmm) ? wpmm : nothing

    fitconst_vec = Vector{FitRes}()
    fitstate_vec = Vector{FitRes}()
    for p0 in plconst0_vec
        # plateau = [p0,pl_f]
        m_data = m_obs[p0:pl_f] 
        fit = fit_routine(const_model,collect(p0:pl_f).+1.5, m_data, 1, pval=pval, wpm=wpm, info=fitinfo, lineprint=fitinfo)
        push!(fitconst_vec,fit)
    end

    if state_fit
        @. state_model(x,p) = p[1] + p[2] * exp(- p[3] * x)
        # @. state_model(x,p) = (p[1] + (p[1] + p[3]) * p[2] * exp(- p[3] * x)) / (1 + p[2] * exp(- p[3] * x))

        plstate0_vec = collect(max(floor(Int64,pl2state0[1]*len),1):plstep:ceil(Int64,pl2state0[end]*len))

        p0erratics = []
        for p0 in plstate0_vec
            try
                m_data = m_obs[p0:pl_f]
                fit = fit_routine(state_model,collect(p0:pl_f).+1.5, m_data, 3, pval=pval, wpm=wpm, info=fitinfo, lineprint=fitinfo)
                push!(fitstate_vec,fit)
            catch
                push!(p0erratics,p0)
                # filter!(x -> x != p0, plstate0_vec)
            end
        end

        plstate0_vec = filter(x -> !(x in p0erratics), plstate0_vec)

        p0_vec  = [plstate0_vec,plconst0_vec]
        fit_vec = [fitstate_vec,fitconst_vec]
    else
        p0_vec  = plconst0_vec
        fit_vec = fitconst_vec
    end

    w = get_w_from_fitres(vcat(fit_vec...), AIC=AIC)
    res_vec = [par[1] for par in getfield.(vcat(fit_vec...),:param)]
    meff_res, meff_sys = model_average(res_vec, w)
    meff_res = meff_res[1]

    if plot
        fig = figure(figsize=(16,12))
        # subplots_adjust(hspace=0.1)
        gs = fig.add_gridspec(4, 1, height_ratios=[4, 1, 1, 1])  # Adjust the height_ratios as needed


        ax1 = fig.add_subplot(gs[1, 1])
        x0 = collect(max(floor(Int64,len*(pl2state0[1]))+1-2,1):len+1) .+ 0.5
        m_vec = m_obs[max(floor(Int64,len*(pl2state0[1]))-2,1):end]; uwerr.(m_vec)
        res = meff_res + uwreal([0.0,meff_sys],"meff MA syst"); uwerr(res)
        # res = meff_res; uwerr(res)

        title("$(corr.id)")
        errorbar(x0, value.(m_vec), err.(m_vec), fmt="d", mfc="none", capsize=2, color="black")
        fill_between(x0, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")
        if state_fit
            maxw_arg = argmax(w[1:length(fit_vec[1])])
            par = getfield.(fit_vec[1],:param)[maxw_arg]
            x_ = collect(max(floor(Int64,len*(pl2state0[1]))+1-2,1):0.1:len+1) .+ 0.5
            y2st_fit = state_model(x_,par); uwerr.(y2st_fit)
            fill_between(x_, value.(y2st_fit)+err.(y2st_fit), value.(y2st_fit)-err.(y2st_fit), alpha=0.3, color="orange")

            axvline(x=plstate0_vec[1]+1+0.1, color="blue", linestyle="--")
            axvline(x=plstate0_vec[end]+2-0.1, color="blue", linestyle="--")
        end
        axvline(x=plconst0_vec[1]+1+0.1, color="red", linestyle="--")
        axvline(x=plconst0_vec[end]+2-0.1, color="red", linestyle="--")
        if plf != 1.0
            axvline(x=ceil(Int64,plf*len)+2, color="orange", linestyle=":")
        end
        axis("tight")
        ylabel(L"$m_{\rm{eff}}$")
        ylim(res.mean-3*res.err,maximum(value.(m_vec[1:3])))
        setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels

        ax2 = fig.add_subplot(gs[2, 1])
        uwerr.(res_vec)
        if state_fit
            errorbar(plstate0_vec .+ 1.5, value.(res_vec[1:length(plstate0_vec)]), err.(res_vec[1:length(plstate0_vec)]), fmt="^", mfc="none", capsize=2, color="blue")
            errorbar(plconst0_vec .+ 1.5, value.(res_vec[end-length(plconst0_vec)+1:end]), err.(res_vec[end-length(plconst0_vec)+1:end]), fmt="^", mfc="none", capsize=2, color="red")
        else
            errorbar(plconst0_vec .+ 1.5, value.(res_vec), err.(res_vec), fmt="^", mfc="none", color="red")
        end
        fill_between(x0, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")
        ylabel(L"m_{D_s}")
        setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels

        ax3 = fig.add_subplot(gs[3, 1])
        fill_between(x0, maximum(w)/2, maximum(w)/2, alpha=0.0, color="white")
        if state_fit
            PyPlot.plot(plstate0_vec .+ 1.5, w[1:length(plstate0_vec)], linestyle="none", marker="o", mfc="none", color="blue")
            PyPlot.plot(plconst0_vec .+ 1.5, w[end-length(plconst0_vec)+1:end], linestyle="none", marker="o", mfc="none", color="red")
        else
            PyPlot.plot(plconst0_vec .+ 1.5, w, linestyle="none", marker="o", mfc="none", color="red")
        end
        ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))

        if pval
            setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels
    
            ax4 = fig.add_subplot(gs[4, 1])
            if state_fit
                pval_vec = [getfield.(fit_vec_,:pval) for fit_vec_ in fit_vec]
                fill_between(x0, maximum(vcat(pval_vec...))/2, maximum(vcat(pval_vec...))/2, alpha=0.0, color="white")
                PyPlot.plot(plstate0_vec .+ 1.5, pval_vec[1], linestyle="none", marker="o", mfc="none", color="blue")
                PyPlot.plot(plconst0_vec .+ 1.5, pval_vec[2], linestyle="none", marker="o", mfc="none", color="red")
            else
                pval_vec = getfield.(fit_vec,:pval)
                fill_between(x0, maximum(vcat(pval_vec...))/2, maximum(vcat(pval_vec...))/2, alpha=0.0, color="white")
                PyPlot.plot(plconst0_vec .+ 1.5, pval_vec, linestyle="none", marker="o", mfc="none", color="red")
            end
            ylabel(L"$\rm{p-values}$")
        end
        xlabel(L"$t/a$")


        tight_layout()
        display(fig)
        close("all")

    end

    if returnfitMA
        return [meff_res, meff_sys], [p0_vec, fit_vec, w]
    else
        return [meff_res, meff_sys]
    end
end


##-- 2D window structure

@doc raw"""

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

# systematic and FVC

function apply_syst_HVP(hvp::Dict,syst::Dict,diag::String,wind::String,ensid::String;systname::String="")::Dict
    hvpkeys  = keys(hvp)
    systkeys = keys(syst)

    HVP = Dict()
    for key in hvpkeys
        if key in systkeys
            if hvp[key] == uwreal
                HVP[key] = hvp[key] + uwreal([0.0,syst[key]],"HVP syst. $systname [$ensid-$diag,$wind,$key]")
            else
                HVP[key] = hvp[key] .+ [uwreal([0.0,syst[key][i]],"HVP syst. $systname [$ensid-$diag,$wind,$key]") for i=1:length(syst[key])]
            end
        else
            HVP[key] = hvp[key]
        end
    end
    return HVP
end

function apply_syst_HVP!(HVP::Dict,syst::Dict,diag::String,wind::String,ensid::String;systname::String="")
    systkeys = keys(syst)

    for key in systkeys
        if HVP[key] == uwreal
            HVP[key] += uwreal([0.0,syst[key]],"HVP syst. $systname [$ensid-$diag,$wind,$key]")
        elseif HVP[key] == Vector{uwreal}
            HVP[key] .+= [uwreal([0.0,syst[key][i]],"HVP syst. $systname [$ensid-$diag,$wind,$key]") for i=1:length(syst[key])]
        end
    end
end

function apply_syst_FVC(fvc::Dict,diag::String,wind::String,ensid::String;factor::Float64=0.25,IMPR_SET::Vector{String}=["1","2"])::Dict
    fvckeys  = keys(fvc)

    FVC = Dict()
    if diag != "NLOc"
        for key in fvckeys
            FVC[key] = fvc[key] .+ [uwreal([0.0,factor*value(fvc[key][i])],"FVC syst.  [$ensid-$diag,$wind,$key]") for i=1:length(fvc[key])]
        end
    else
        for impr_set = IMPR_SET
            for key in fvckeys
                FVC[impr_set][key] = fvc[impr_set][key][end] + uwreal([0.0,factor*value(fvc[impr_set][key][end])],"FVC syst.  [$ensid-$diag,$wind,$key]")
            end
        end
    end
    return FVC
end

function apply_syst_FVC!(FVC::Dict,diag::String,wind::String,ensid::String;factor::Float64=0.25,IMPR_SET::Vector{String}=["1","2"])
    fvckeys  = keys(FVC)

    if diag != "NLOc"
        for key in fvckeys
            FVC[key] = FVC[key] .+ [uwreal([0.0,factor*value(FVC[key][i])],"FVC syst.  [$ensid-$diag,$wind,$key]") for i=1:length(FVC[key])]
        end
    else
        for impr_set = IMPR_SET
            for key in fvckeys
                FVC[impr_set][key] = FVC[impr_set][key] .+ [uwreal([0.0,factor*value(FVC[impr_set][key][i])],"FVC syst.  [$ensid-$diag,$wind,$key]") for i=1:length(FVC[key][impr_set])]
            end
        end
    end
    return FVC
end

function HVP_VolCorrect(HVP::Dict,FVC::Dict,diag::String;IMPR_SET::Vector{String}=["1","2"])::Dict
    hvpkeys = keys(HVP[IMPR_SET[1]])
    DISCR = ["ll","lc"]
    if diag != "NLOc"
        COMP = ["g33"]
        "g88_ll" in hvpkeys ? push!(COMP,"g88") : nothing
        "∆ls_amu_ll" in hvpkeys ? push!(COMP,"∆ls_amu") : nothing
        "∆lc_b_ll" in hvpkeys ? push!(COMP,"∆lc_b") : nothing
    else
        COMP = ["3333","8888","3388","33CC","88CC"]
    end

    amu_ens = deepcopy(HVP)

    for impr_set in IMPR_SET
        for comp in COMP
            for discr in DISCR
                if diag != "NLOc"
                    if typeof(amu_ens[impr_set]["$(comp)_$(discr)"]) == uwreal
                        amu_ens[impr_set]["$(comp)_$(discr)"] = HVP[impr_set]["$(comp)_$(discr)"] + FVC["FVC$(comp)"][end]
                    elseif typeof(amu_ens[impr_set]["$(comp)_$(discr)"]) == Vector{uwreal}
                        amu_ens[impr_set]["$(comp)_$(discr)"] = HVP[impr_set]["$(comp)_$(discr)"] .+ FVC["FVC$(comp)"]
                    else
                        @warn("FVC could not be applied to comp $comp impr. set $impr_set")
                    end
                else
                    amu_ens[impr_set]["g$(comp)_$(discr)"] = HVP[impr_set]["g$(comp)_$(discr)"] + FVC[impr_set]["FVC$(comp)_$(discr)"][end]
                end
            end
        end
    end
    return amu_ens
end

function HVP_VolCorrect!(HVP::Dict,FVC::Dict,diag::String;IMPR_SET::Vector{String}=["1","2"])
    hvpkeys = keys(HVP[IMPR_SET[1]])
    DISCR = ["ll","lc"]
    if diag != "NLOc"
        COMP = ["g33"]
        "g88_ll" in hvpkeys ? push!(COMP,"g88") : nothing
        "∆ls_amu_ll" in hvpkeys ? push!(COMP,"∆ls_amu") : nothing
        "∆lc_b_ll" in hvpkeys ? push!(COMP,"∆lc_b") : nothing
    else
        COMP = ["3333","8888","3388","33CC","88CC"]
    end

    for impr_set in IMPR_SET
        for comp in COMP
            for discr in DISCR
                if diag != "NLOc"
                    if typeof(HVP[impr_set]["$(comp)_$(discr)"]) == uwreal
                        HVP[impr_set]["$(comp)_$(discr)"] += FVC["FVC$(comp)"][end]
                    elseif typeof(HVP[impr_set]["$(comp)_$(discr)"]) == Vector{uwreal}
                        HVP[impr_set]["$(comp)_$(discr)"] .+= FVC["FVC$(comp)"]
                    else
                        @warn("FVC could not be applied to comp $comp impr. set $impr_set")
                    end
                else
                    HVP[impr_set]["g$(comp)_$(discr)"] += FVC[impr_set]["FVC$(comp)_$(discr)"][end]
                end
            end
        end
    end
end

function HVP_3limpr!(HVP::Dict,HVP3l0::Union{Vector{Float64},Float64};IMPR_SET::Vector{String}=["1","2"],meth::String="prod")
    for impr_set in IMPR_SET
        for discr in ["ll","lc"]
            if meth == "prod"
                if typeof(HVP[impr_set]["g33_$discr"]) == Vector{uwreal} && typeof(HVP3l0) == Vector{Float64}
                    HVP[impr_set]["g33_$discr"] .*= HVP3l0 ./ HVP[impr_set]["g33tl_$discr"]
                elseif typeof(HVP[impr_set]["g33_$discr"]) == uwreal && typeof(HVP3l0) == Float64
                    HVP[impr_set]["g33_$discr"]  *= HVP3l0  / HVP[impr_set]["g33tl_$discr"]
                else
                    @warn("3l impr. could not be applied")
                end
            elseif meth == "sum"
                if typeof(HVP[impr_set]["g33_$discr"]) == Vector{uwreal} && typeof(HVP3l0) == Vector{Float64}
                    HVP[impr_set]["g33_$discr"] .+= HVP3l0 .- HVP[impr_set]["g33tl_$discr"]
                elseif typeof(HVP[impr_set]["g33_$discr"]) == uwreal && typeof(HVP3l0) == Float64
                    HVP[impr_set]["g33_$discr"]  += HVP3l0  - HVP[impr_set]["g33tl_$discr"]
                else
                    @warn("3l impr. could not be applied")
                end
            else
                error("Method not recognised, please choose between 'mult' or 'sum'")
            end
        end
    end
end

# manipulate derivatives and errors

function add_t0_err!(obs::Vector{uwreal}, t0phys::uwreal)
    uwerr(t0phys); uwerr.(obs)
    for k in eachindex(obs)
        err_t0 = abs.(mchist(obs[k], "sqrtt0 [fm]") * err(t0phys) / artificial_err)
        obs[k] = obs[k] + uwreal([0.0, err_t0[1]], "sqrtt0 [fm]")
    end
    return nothing
end

function get_t0err(obs::Vector{uwreal}, t0phys::uwreal)
    uwerr(t0phys); uwerr.(obs)
    err_t0 = []
    for k in eachindex(obs)
        push!(err_t0, abs.(mchist(obs[k], "sqrtt0 [fm]")[1] * err(t0phys) / artificial_err))
    end
    return err_t0
end


# tree level iso-vector

function treelevel_continuum_correlator(t)
    return 1. / (2* pi^2 * t^3)
end


function compute_HVPtl0(diag::String,wind::String,Qlist::Union{Vector{Float64},Vector{Int64}},path_coef::String)
    exp_diag = diag == "LO" ? 2 : 3
    Tildef = Dict("LO" => Tildef2, "NLOa" => Tildef4a, "NLOb" => Tildef4b, "NLOa&b" => (x,path) -> Tildef4a(x,path) + Tildef4b(x,path))

    f = nothing  # initialize f

    if wind == "SDsub"
        tl_cont = Float64[]
        for Q in Qlist
            f = x0 -> treelevel_continuum_correlator(x0) .* (Window("SD")(x0) * hbarc^2 * Tildef[diag](massmu/hbarc * x0,path_coef) - (Window("SD")(0) * (16/(Q/hbarc)^2)^2 * π^2 * (massmu/hbarc)^2 * C4[diag](massmu/hbarc * x0) * sin((Q/hbarc/4) * x0)^4))
            res, _ = quadgk(f, 0, 5,rtol=1e-3)
            amu = (alpha/pi)^exp_diag * res * 1e10 / 2
            push!(tl_cont, amu)
        end
    elseif wind == "SD"
        f = x0 -> treelevel_continuum_correlator(x0) .* (Window("SD")(x0) * hbarc^2 * Tildef[diag](massmu/hbarc * x0,path_coef))
        res, _ = quadgk(f, 0, 5,rtol=1e-3)
        tl_cont = (alpha/pi)^exp_diag * res * 1e10 / 2
    else
        error("3l improvement cannot be applied to wind $wind")
    end
    return tl_cont
end


# Bounding method functions


function Eeff(tstar::Int64, obs::Vector{uwreal})

    E_eff = 0.5 * log((obs[tstar+1] / obs[tstar+2]) ^2)
    return E_eff
end

function corr_bound(t::Vector{Int64}, tcut::Int64, obs::Vector{uwreal}, ens::EnsInfo, Eeff::Union{uwreal,Float64}; allowPBC::Bool=false)
    T = 2*t[end]

    Gcut = obs[tcut]
    GPBC(x0) = exp(-Eeff*x0) + exp(-Eeff*(T-x0))
    G(x0) = exp(-Eeff*x0)

    if split(ens.id,"")[end-1]=="5" && allowPBC
        UBarray = (Gcut/GPBC(tcut)) .*  GPBC.(t[tcut:end])
    else
        UBarray = Gcut/G(tcut) .*  G.(t[tcut:end])
    end
    return UBarray[2:end]
end
corr_bound(t::Vector{Int64}, tcut::Int64, corr::Corr, ens::EnsInfo, Eeff::Union{uwreal,Float64}; allowPBC::Bool=false) = corr_bound(t, tcut, corr.obs, ens, Eeff, ; allowPBC=allowPBC)


# function buonding_method(ub::Vector{uwreal},lb::Vector{uwreal};PLAT::Bool=false,AVER::Bool=false)
#     averb = (ub.+lb)./2; uwerr.(averb)
#     x0 = findfirst(abs.(value.(ub).-value.(lb)) .< 0.5.*err.(averb))
#     ∆x = findfirst(abs.(averb[x0:end].-averb[x0]) .> err(averb[x0]))
#     if isnothing(∆x) || ∆x > 5
#         hvp  = averb[x0]
#         syst = 0.0
#         xend = isnothing(∆x) ? length(averb) : (x0-1) + ∆x
#     else
#         xend=x0+3
#         x0 > 3 ? (x0-=3) : (x0=1)
#         hvp  = sum(averb[x0:xend])/length(averb[x0:xend])
#         aux1 = sum(averb[x0:xend].^2)/length(averb[x0:xend])
#         aux2 = hvp^2
#         syst = sqrt(abs(value(aux1 - aux2)))
#     end
#     # plateau_fm = value(aEns).*(collect(x0:xend).+tcut0.-2)

#     returnVec = [hvp,syst]
#     PLAT ? push!(returnVec,[x0,xend]) : nothing
#     AVER ? push!(returnVec,averb) : nothing
#     return returnVec
# end

function buonding_method(ub::Vector{uwreal},lb::Vector{uwreal},aEns::uwreal;PLAT::Bool=false,AVER::Bool=false,tcut0::Union{Int64,Nothing}=nothing)
    averb = (ub.+lb)./2; uwerr.(averb)
    x0   = findfirst(abs.(value.(ub).-value.(lb)) .< err.(averb))
    xend = (4*value(aEns) > 0.25) ? Int64(x0 + round(0.25/value(aEns),RoundUp)) : Int64(x0 + 4)
    hvp  = mean(averb[x0:xend])
    aux1 = mean(averb[x0:xend].^2)
    aux2 = hvp^2
    syst = sqrt(abs(value(aux1 - aux2)))

    if PLAT
        !isnothing(tcut0) ? plateau_fm = value(aEns).*(collect(x0:xend).+tcut0.-2) : error("tcut0 required to output plateau")
    end

    returnVec = [hvp,syst]
    PLAT ? push!(returnVec,plateau_fm) : nothing
    AVER ? push!(returnVec,averb) : nothing
    return returnVec
end


function findfirst_uninterrupted(bools)
    i = 0
    stop = false
    while !stop
        arg = findfirst(bools[i+1:end])
        i = !isnothing(arg) ? (i+arg) : nothing
        if isnothing(i) || isnothing(findfirst(.!bools[i+1:end]))  # stops if there are unterrupted trues or if it reached the end
            stop = true
        end
    end
    return i
end


function get_spectr_data(path::String,ens::EnsInfo)
    ph5 = joinpath(path,"$(ens.id)-pSq0-T1up.h5")
    file = h5open(ph5, "r")

    # file["summary"]

    # Read the energies
    EnKeys = HDF5.keys(file["pSq0-T1up/spectrum"])
    E = uwreal[]
    for n in EnKeys
        En = read(file["pSq0-T1up/spectrum/$n"])
        push!(E,uwreal(jackknife_err(En),n))
    end

    # Read the overlaps
    ZKeys = HDF5.keys(file["pSq0-T1up/overlaps"])
    Z = uwreal[]
    Z_impr = uwreal[]
    for p in ZKeys
        Zp = read(file["pSq0-T1up/overlaps/$p"])
        Zp_impr = read(file["pSq0-T1up/overlaps_imp/$p"])
        push!(Z,uwreal(jackknife_err(Zp),p))
        push!(Z_impr,uwreal(jackknife_err(Zp_impr),p))
    end

    return E, Z, Z_impr
end

# function corr_n(n::Int64,Union{Vector{Int64},Vector{Float64}},E::Vector{uwreal},Z::Vector{uwreal},L::Int64)
#     return Z[n] / (L^3) .* exp.(-E[n] .* t)
# end
function corr_n(n::Int64,t::Union{Vector{Int64},Vector{Float64}},E::Vector{uwreal},Z::Vector{uwreal},L::Int64)
    return (1/2) * Z[n] / (L^3) .* exp.(-E[n] .* t) # the 1/2 factor convers light-light corr to the iso-vector corr
end

function reconstr_corr(
        ens::EnsInfo,
        E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};
        nmax::Union{Nothing,Int64}=nothing,
        impr_set::String="1",
        IMPR::Bool=true,
        RENORM::Bool=true,
        FreeEN::Bool=false,
        total::Bool=false,
    )
    length(Z) != length(Zimpr) ? error("Length of Z and Zimpr are not the same") : nothing
    # nmax chooses the tower of states to be added
    # if nmax is not defined, then the max number of available states are used
    # FreeEN ensures that the last energy level is not used when saturating the corr; usefull for impr BM
    if isnothing(nmax)
        FreeEN ? nmax = minimum([length(E)-1,length(Z)]) : nmax = minimum([length(E),length(Z)])
    else
        if FreeEN
            nmax > length(E)-1 ? error("Input E data not large enough to tower $nmax (-1) states") : nothing
        else
            nmax > length(E) ? error("Input E data not large enough to tower $nmax states") : nothing
        end
        nmax > length(Z) ? error("Input Z data not large enough to tower $nmax states") : nothing
    end
    
    t = collect(1:Int64(HVPobs.Data.get_T(ens.id))/2+1)
    corrVec_PiPi = [corr_n(1,t.-1,E,Z.^2 ,ens.L)]
    corrVec_Impr = [corr_n(1,t.-1,E,Zimpr.*Z,ens.L)]
    for n=2:nmax
        push!(corrVec_PiPi,corrVec_PiPi[n-1] .+ corr_n(n,t.-1,E,Z.^2 ,ens.L))
        push!(corrVec_Impr,corrVec_Impr[n-1] .+ corr_n(n,t.-1,E,Zimpr.*Z,ens.L))
    end

    # corrVec = deepcopy(corrVec_PiPi)
    if IMPR
        if impr_set == "1"
            cv_l = cv_loc(ens.beta)
        elseif impr_set == "1old"
            cv_l = cv_loc_old(ens.beta)
        elseif impr_set == "2"
            cv_l = cv_loc_set2(ens.beta)
        else
            error("Impr Set $impr_set not recoginsed, please choose between '1', '1old' and '2'.")
        end
        # [1:Int64(HVPobs.Data.get_T(ens.id)/2+1)]
        # [improve_corr_vkvk!(corrVec[n], corrVec_JPi[n], 2*cv_l, std=std, treelevel=true) for n in collect(1:nmax)]
        corrVec = [corrVec_PiPi[n] .+ (2*cv_l).*corrVec_Impr[n] for n=1:nmax]
    else
        corrVec = corrVec_PiPi
    end

    if RENORM
        Z3 = get_Z3(ens, impr_set=impr_set)
        [renormalize!(corr, Z3^2) for corr in corrVec]
        if !total
            [renormalize!(corr, Z3^2) for corr in corrVec_PiPi]
            [renormalize!(corr, Z3^2) for corr in corrVec_Impr]
        end
    end

    !total ? (return corrVec) : (return corrVec, corrVec_PiPi, corrVec_Impr)
end
reconstr_corr(ensid::String,E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};nmax::Union{Nothing,Int64}=nothing,impr_set::String="1",IMPR::Bool=true,RENORM::Bool=true,FreeEN::Bool=false,total::Bool=false,) = reconstr_corr(EnsInfo(ensid),E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};nmax=nmax,impr_set=impr_set,IMPR=IMPR,RENORM=RENORM,FreeEN=FreeEN,total=total,)


# Unc functions

function jackknife_resampling(data)
    N = length(data)
    means = zeros(N)
    
    for i in 1:N
        sample = deleteat!(copy(data), i)  # Leave-one-out sample
        means[i] = mean(sample)
    end
    
    mean_jack = mean(means)
    std_jack = sqrt((N - 1) * mean((means .- mean_jack) .^ 2))
    
    return [mean_jack, std_jack]
end

function bootstrap_resampling(data, B=1000)
    N = length(data)
    means = zeros(B)
    
    for i in 1:B
        resample = data[rand(1:N, N)]  # Resample with replacement
        means[i] = mean(resample)
    end
    
    mean_boot = mean(means)
    std_boot = std(means)  # Standard deviation of bootstrap samples
    
    return [mean_boot, std_boot]
end

function jackknife_err(jakknife_sample,∆dof=0)
    N    = length(jakknife_sample)
    Mean = mean(jakknife_sample)
    Err  = sqrt(N*(N-1)/(N-∆dof)) * sqrt(mean((jakknife_sample.-Mean).^2))
    return [Mean, Err]
end