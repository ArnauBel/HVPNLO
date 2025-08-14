# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO

using Suppressor

using Plots
using PyPlot
using Colors

# Path definition

julia_script_directory = @__DIR__

path_heavy = joinpath(julia_script_directory, "..", "..", "LMEData", "s_heavy_data")

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

# using OrderedCollections

include("../tools/HVPobs_ext.jl")
include("../tools/fit&MA.jl")

# Ensamble choice

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302"]
ensInfo = EnsInfo.(ensList)

# physical MDs mass

artificial_err = 10^(-8)
# const MD_ph = uwreal([1968.47*1e-3,artificial_err], "MD_ph [GeV]")
const MD_ph = uwreal([1968.47*1e-3,artificial_err], "MD_ph [GeV]")

filterInfo = (getfield.(ensInfo,:beta) .!= 3.34) .& (getfield.(ensInfo,:beta) .!= 3.85)
ensInfo = ensInfo[filterInfo]

@info("No ensembles with ß = 3.34 or 3.85 can be considered in this analysis. The ensembles considered are: \n - $(getfield.(ensInfo,:id))")

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

NOERR_MDs = true

## Find MDs prime values in lattice units as a function of beta

MD_ph_prime = MD_ph * (t0_ph_Bruno / t0_ph)

uwerr(MD_ph_prime); println(" - MD' shifts to $MD_ph_prime GeV")

aMD_ph_prime = Dict{Float64, uwreal}()
for beta in b_values[2:end-1]
    aMD_ph_prime[beta] = MD_ph_prime * (t0_ph / sqrt(t0sym(beta,Bruno=true))) / hbarc
end

[uwerr(aMD_ph_prime[key]) for key in keys(aMD_ph_prime)]

if NOERR_MDs
    [set_fluc_to_zero!(aMD_ph_prime[key], "sqrtt0 [fm] (Bruno)") for key in keys(aMD_ph_prime)]
    [set_fluc_to_zero!(aMD_ph_prime[key], "t0sym/a2") for key in keys(aMD_ph_prime)]
    [aMD_ph_prime[key] *= 1.0 for key in keys(aMD_ph_prime)] 
end

[uwerr(aMD_ph_prime[key]) for key in keys(aMD_ph_prime)]

aMD_ph_prime

##

ens = EnsInfo("D200")

kappa_c = get_kappa_values()

data_sheavy = read_kappa_charm_all_config(joinpath(path_heavy,ens.id))


pl_step = 1
@. const_model(x,p) = p[1] + 0*x
meff_res = Dict{String,uwreal}(); meff_sys = Dict{String,uwreal}()
plateau_vec = Dict{String,Vector{Vector{Int64}}}()
fit_vec     = Dict{String,Vector{FitRes}}()
for key in ["sh1","sh2","sh3","sh4"]
    println(" - Strange - heavy ($(key[end]))")
    corr_ = corr_obs(data_sheavy[key])
    len_ = Int64(length(corr_.obs)/2+1)
    plateau_0 = collect(Int64(0.5*(len_-1)):pl_step:ceil(Int64,0.65*len_))
    plateau_f = collect(floor(Int64,0.70*len_):pl_step:ceil(Int64,0.85*len_))
    plateau_vec[key] = Vector{Vector{Float64}}()
    fit_vec[key]     = Vector{FitRes}()
    for p0 in plateau_0
        for pf in plateau_f
            plateau = [p0,pf]
            meff_data = meff(corr_.obs)[p0+1:pf+1]
            fit_p = fit_routine(const_model,collect(p0:pf), meff_data, 1)
            push!(plateau_vec[key],plateau)
            push!(fit_vec[key],fit_p)
        end
    end
    w = get_w_from_fitres(fit_vec[key])
    res = [res_vec[1] for res_vec in getfield.(fit_vec[key],:param)]
    meff_res[key], meff_sys[key] = model_average(res, w)
end

kappac_list_inv = 1 ./ kappa_c[ens.id][2:end]
meff_list = [meff_res[key] + uwreal([0.0,value(meff_sys[key])],"meff MA syst") for key in ["sh1","sh2","sh3","sh4"]]

@. lin_model(x,p) = p[1] + p[2] * x

fit = fit_routine(lin_model, kappac_list_inv, meff_list, 2)
par = fit.param

kappac_inv = (aMD_ph_prime[ens.beta] - par[1]) / par[2]

kappac = 1/kappac_inv; uwerr(kappac)
details(kappac)
meff_list

sh = "sh4"
fit_vec[sh][argmax(get_w_from_fitres(fit_vec[sh]))].param
get_w_from_fitres(fit_vec[sh])[argmax(get_w_from_fitres(fit_vec[sh]))]

##

kappac_arr = range(kappac_list_inv[1],kappac_list_inv[end],100)
mDs_arr = lin_model(kappac_arr,par); uwerr.(mDs_arr)
uwerr(aMD_ph_prime[ens.beta]); uwerr(kappac_inv)

##

uwerr(kappac)
title("Ens: $(ens.id); "*L"$\kappa_c$"*" = $(round(value(kappac),digits=6)) ± $(round(err(kappac),digits=6))")
errorbar(kappac_list_inv, value.(meff_list), err.(meff_list), fmt="o", mfc="none", capsize=2, color="blue")
errorbar(value(kappac_inv), value(aMD_ph_prime[ens.beta]), xerr = err(kappac_inv), yerr = 0.0, fmt="^", mfc="none", capsize=2, color="red")
# errorbar(value(kappac_inv), 0.8615, xerr = err(kappac_inv), yerr = 0.0, fmt="", mfc="none", capsize=2, color="red")
fill_between(kappac_arr, value.(mDs_arr)+err.(mDs_arr), value.(mDs_arr)-err.(mDs_arr), alpha=0.4, color="blue")
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
xlabel(L"$1/\kappa_c$")
ylabel(L"$am_{D_s}$")
display(gcf())      #display the figure
close()

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

NOERR_MDs = true
PLOT_meff = true
PLOT_Ktar = true

## Find MDs prime values in lattice units as a function of beta

MD_ph_prime = MD_ph * (t0_ph_Bruno / t0_ph)

uwerr(MD_ph_prime); println(" - MDs shifts to MDs' = $MD_ph_prime GeV")

aMD_ph_prime = Dict{Float64, uwreal}()
for beta in b_values[2:end-1]
    aMD_ph_prime[beta] = MD_ph_prime * (t0_ph / sqrt(t0sym(beta,Bruno=true))) / hbarc
end
aMD_ph_prime
[uwerr(aMD_ph_prime[key]) for key in keys(aMD_ph_prime)]

if NOERR_MDs
    [set_fluc_to_zero!(aMD_ph_prime[key], "sqrtt0 [fm] (Bruno)") for key in keys(aMD_ph_prime)]
    [set_fluc_to_zero!(aMD_ph_prime[key], "t0sym/a2") for key in keys(aMD_ph_prime)]
    [aMD_ph_prime[key] *= 1.0 for key in keys(aMD_ph_prime)] 
end

# Read s-heavy data and find kappa_c target using MDs prime

kappa_c = get_kappa_values()

meff_list = Dict{String,Vector{uwreal}}()
kappa_c_tar = Dict{String,uwreal}()
for ens in ensInfo
    println("- Ens: $(ens.id)")

    data_sheavy = read_kappa_charm_all_config(joinpath(path_heavy,ens.id))
    
    meff_ = Dict{String,uwreal}()
    for key in ["sh1","sh2","sh3","sh4"]
        corr_ = corr_obs(data_sheavy[key])
        len_ = Int64(length(corr_.obs)/2+1)
        meff_[key] = meff(corr_.obs[1:len_],[Int64(0.5*(len_-1)),ceil(Int64,0.75*len_)],pl=PLOT_meff)
    end
    
    kappac_list_inv = 1 ./ kappa_c[ens.id][2:end]
    meff_list[ens.id] = [meff_[key] for key in ["sh1","sh2","sh3","sh4"]]
    
    @. lin_model(x,p) = p[1] + p[2] * x
    
    fit = fit_routine(lin_model, kappac_list_inv, meff_list[ens.id], 2)
    fit.chi2/fit.dof
    par = fit.param
    
    kappac_inv = (aMD_ph_prime[ens.beta] - par[1]) / par[2]
    
    kappa_c_tar[ens.id] = 1/kappac_inv
    
    if PLOT_Ktar
        kappac_arr = range(kappac_list_inv[1],kappac_list_inv[end],100)
        mDs_arr = lin_model(kappac_arr,par); uwerr.(mDs_arr)
        uwerr(aMD_ph_prime[ens.beta]); uwerr(kappac_inv)
        uwerr(kappa_c_tar[ens.id])
        
        title("Ens: $(ens.id); "*L"$\kappa_c$"*" = $(round(value(kappa_c_tar[ens.id]),digits=6)) ± $(round(err(kappa_c_tar[ens.id]),digits=6))")
        errorbar(kappac_list_inv, value.(meff_list[ens.id]), err.(meff_list[ens.id]), fmt="o", mfc="none", capsize=2, color="blue")
        errorbar(value(kappac_inv), value(aMD_ph_prime[ens.beta]), xerr = err(kappac_inv), yerr = 0.0, fmt="^", mfc="none", capsize=2, color="red")
        fill_between(kappac_arr, value.(mDs_arr)+err.(mDs_arr), value.(mDs_arr)-err.(mDs_arr), alpha=0.4, color="blue")
        axis("tight")
        ax = gca()      # get the handle of the current axis (not really used here)
        xlabel(L"$1/\kappa_c$")
        ylabel(L"$am_{D_s}$")
        display(gcf())      #display the figure
        close()
    end
end

meff_list
kappa_c_tar


##

function get_der_from_uwerr(a::uwreal, id_str::String)
    ws = ADerrors.wsg
    id_int = ws.str2id[id_str]
    idx = ADerrors.find_mcid(a, id_int)
    if isnothing(idx)
        error("No error available... maybe run uwerr")
    else
        nd = ws.fluc[ws.map_ids[a.ids[idx]]].nd
        for j in 1:length(a.prop)
            if (a.prop[j] && ((ws.map_nob[j] == a.ids[idx])))
                ws.fluc[j].delta[:] .= 0.0
            end
        end
        return nothing
    end
    
end

##

function meff_MA_(corr::Corr; 
    pl2state0::Vector{Float64}=[0.5,0.6], plconst0::Float64=0.7, plf::Float64=1.0, plstep::Int64=1, 
    mdof::Int64=4, 
    state_fit::Bool=true, 
    AIC::Bool=true, 
    plot::Bool=false, 
    return_fitMA::Bool=false, 
    pval::Bool=false
    )
    obs = corr.obs
    m_obs = meff(obs[1:Int64(length(obs)/2+1)])  # first and last data points are lost (-2) and derivative is taken (-1); length(meff_) = length(obs) - 3
    len = length(m_obs)

    @. const_model(x,p) = p[1] + 0*x

    plconst0_vec = collect(floor(Int64,plconst0*len):plstep:ceil(Int64,plf*len-mdof))
    pl_f = ceil(Int64,plf*len)

    if isempty(plconst0_vec)
        error("Not enough data to fit with d.o.f. ≥ $mdof. Please decrease $mdof or make increase the plateau search space. ")
    end

    fitconst_vec = Vector{FitRes}()
    fitstate_vec = Vector{FitRes}()
    for p0 in plconst0_vec
        # plateau = [p0,pl_f]
        m_data = m_obs[p0:pl_f] 
        fit = fit_routine(const_model,collect(p0:pl_f), m_data, 1, pval=pval)
        push!(fitconst_vec,fit)
    end

    if state_fit
        @. state_model(x,p) = p[1] + p[2] * exp(- p[3] * x)

        plstate0_vec = collect(floor(Int64,pl2state0[1]*len):plstep:ceil(Int64,pl2state0[end]*len-mdof))

        for p0 in plstate0_vec
            # plateau = [p0,pl_f]
            m_data = m_obs[p0:pl_f]
            fit = fit_routine(state_model,collect(p0:pl_f), m_data, 3, pval=pval)
            push!(fitstate_vec,fit)
        end

        p0_vec  = [plstate0_vec,plconst0_vec]
        fit_vec = [fitstate_vec,fitconst_vec]
    else
        p0_vec  = plconst0_vec
        fit_vec = fitconst_vec
    end


    w = get_w_from_fitres(vcat(fit_vec...), AIC=AIC)
    res_vec = [par[1] for par in getfield.(vcat(fit_vec...),:param)]
    meff_res, meff_sys = model_average(res_vec, w)

    if plot
        fig = figure(figsize=(16,12))
        # subplots_adjust(hspace=0.1)
        gs = fig.add_gridspec(3, 1, height_ratios=[4, 1, 1])  # Adjust the height_ratios as needed


        ax1 = fig.add_subplot(gs[1, 1])
        x0 = collect(floor(Int64,len*(pl2state0[1]))-4+1:len+1) .+ 0.5
        m_vec = m_obs[floor(Int64,len*(pl2state0[1]))-4:end]; uwerr.(m_vec)
        res = meff_res + uwreal([0.0,meff_sys],"meff MA syst"); uwerr(res)
        # res = meff_res; uwerr(res)

        title("$(corr.id)")
        errorbar(x0, value.(m_vec), err.(m_vec), fmt="d", mfc="none", capsize=2, color="black")
        fill_between(x0, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")
        if state_fit
            axvline(x=plstate0_vec[1]+1, color="blue", linestyle="--")
            axvline(x=plstate0_vec[end]+2, color="blue", linestyle="--")
        end
        xlim((x0[1]-3,x0[end]+3))
        axvline(x=plconst0_vec[1]+1, color="red", linestyle="--")
        axvline(x=plconst0_vec[end]+2, color="red", linestyle="--")
        axis("tight")
        ylabel(L"$m_{\rm{eff}}$")
        setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels

        ax2 = fig.add_subplot(gs[2, 1])
        if state_fit
            PyPlot.plot(plstate0_vec .+ 1.5, w[1:length(plstate0_vec)], linestyle="none", marker="o", mfc="none", color="blue")
            PyPlot.plot(plconst0_vec .+ 1.5, w[length(plstate0_vec)+1:end], linestyle="none", marker="o", mfc="none", color="red")
        else
            PyPlot.plot(plconst0_vec .+ 1.5, w, linestyle="none", marker="o", mfc="none", color="red")
        end
        xlim((x0[1]-3,x0[end]+3))
        ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))

        if pval
            setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels

    
            ax3 = fig.add_subplot(gs[3, 1])
            if state_fit
                pval_vec = [getfield.(fit_vec_,:pval) for fit_vec_ in fit_vec]
                PyPlot.plot(plstate0_vec .+ 1.5, pval_vec[1], linestyle="none", marker="^", mfc="none", color="blue")
                PyPlot.plot(plconst0_vec .+ 1.5, pval_vec[2], linestyle="none", marker="^", mfc="none", color="red")
            else
                pval_vec = getfield.(fit_vec,:pval)
                PyPlot.plot(plconst0_vec .+ 1.5, pval_vec, linestyle="none", marker="^", mfc="none", color="red")
            end
            xlim((x0[1]-3,x0[end]+3))
            ylabel(L"$\rm{p-values}$")
        end
        xlabel(L"$t/a$")


        tight_layout()
        display(fig)
        close("all")

    end

    if return_fitMA
        return [meff_res, meff_sys], [p0_vec, fit_vec, w]
    else
        return [meff_res, meff_sys]
    end
end

##

ens = EnsInfo("D200")

kappa_c = get_kappa_values()

data_sheavy = read_kappa_charm_all_config(joinpath(path_heavy,ens.id))

corr = corr_obs(data_sheavy["sh3"])

amDs_res, fitMA = meff_MA(corr; state_fit=true, pl2state0=[0.15,0.5], plconst0=0.7, plf=1.0, plstep=1, mdof=4, AIC=true, plot=true, fitMA=true, pval=true)

amDs = amDs_res[1] + uwreal([0.0,amDs_res[2]],"amDs syst"); uwerr(amDs)
amDs
