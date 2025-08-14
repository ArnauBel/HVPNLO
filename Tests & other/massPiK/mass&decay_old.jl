# Import packages

include("HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO

using Plots
using PyPlot
using Colors

# Path definition

julia_script_directory = @__DIR__

path_hvp_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "HVPData"),
    "SSD"   => joinpath(julia_script_directory, "..", "ObsExternal","PortableSSD","HVPData"),
    "clust" => joinpath(julia_script_directory, "..", "ObsExternal", "mogon_mount", "HVPdata")
)

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_meson = joinpath(julia_script_directory, "..", "HVPData", "Meson_data")

path_rw_   = joinpath(path_hvp_dict["local"], "reweight")
path_ms    = joinpath(path_hvp_dict["local"], "ms_t0_dat")

path_fvcPI    = joinpath(path_hvp_dict["local"], "FSE_HP", "inf", "JKMPI_Mvmd")  # _Mvmd
path_fvcK     = joinpath(path_hvp_dict["local"], "FSE_HP", "inf", "JKMK")
path_fvcPIref = joinpath(path_hvp_dict["local"], "FSE_HP", "ref", "JKMPI_Mvmd")
path_fvcKref  = joinpath(path_hvp_dict["local"], "FSE_HP", "ref", "JKMK")

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

# path to ObsBDIO

path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")

# Ensamble choice

# Ens with problems: H105, C101, N202

mPi_fitinfo = Dict{String,Dict{String,Any}}(
    "A653" => Dict{String,Any}("2state" => true, "plat" => [0.60,1.00], "mdof" => 4),
    "A654" => Dict{String,Any}("2state" => true, "plat" => [0.50,1.00], "mdof" => 4),
    "D150" => Dict{String,Any}("2state" => true, "plat" => [0.20,0.40], "mdof" => 4),
    "B450" => Dict{String,Any}("2state" => true, "plat" => [0.30,0.40], "mdof" => 4),
    "N452" => Dict{String,Any}("2state" => true, "plat" => [0.50,0.60], "mdof" => 4),
    "N451" => Dict{String,Any}("2state" => true, "plat" => [0.30,0.40], "mdof" => 4),
    "D450" => Dict{String,Any}("2state" => true, "plat" => [0.80,0.85], "mdof" => 4),
    "D451" => Dict{String,Any}("2state" => true, "plat" => [0.70,0.75], "mdof" => 4),
    "D452" => Dict{String,Any}("2state" => true, "plat" => [0.40,0.50], "mdof" => 4),
    "D251" => Dict{String,Any}("2state" => true, "plat" => [0.30,0.40], "mdof" => 4),
    "E250" => Dict{String,Any}("2state" => true, "plat" => [0.30,0.50], "mdof" => 4),

    "H101" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "H102" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "N101" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.75], "mdof" => 5),
    # "C101" => Dict{String,Any}("2state" => false, "plat" => , "mdof" => ),
    "C102" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.75], "mdof" => 5),
    "S400" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.75], "mdof" => 5),
    "H200" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    # "N202" => Dict{String,Any}("2state" => false, "plat" => , "mdof" => ),
    "N203" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "N200" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 20),
    "D200" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.65], "mdof" => 10),
    "D201" => Dict{String,Any}("2state" => false, "plat" => [0.30,0.60], "mdof" => 10),
    "N300" => Dict{String,Any}("2state" => false, "plat" => [0.30,0.70], "mdof" => 10),
    "J307" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "N302" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "J306" => Dict{String,Any}("2state" => false, "plat" => [0.10,0.75], "mdof" => 80),
    "J303" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "J304" => Dict{String,Any}("2state" => false, "plat" => [0.15,0.80], "mdof" => 90),
    "E300" => Dict{String,Any}("2state" => false, "plat" => [0.10,0.70], "mdof" => 40),
    "F300" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 100),
    "J500" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "J501" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 30),
)

mK_fitinfo = Dict{String,Dict{String,Any}}(
    "A654" => Dict{String,Any}("2state" => true, "plat" => [0.50,1.00], "mdof" => 4),
    "D150" => Dict{String,Any}("2state" => true, "plat" => [0.30,0.50], "mdof" => 4),
    "N452" => Dict{String,Any}("2state" => true, "plat" => [0.50,0.60], "mdof" => 4),
    "N451" => Dict{String,Any}("2state" => true, "plat" => [0.30,0.40], "mdof" => 4),
    "D450" => Dict{String,Any}("2state" => true, "plat" => [0.20,0.40], "mdof" => 4),
    "D451" => Dict{String,Any}("2state" => false, "plat" => [0.35,0.45], "mdof" => 4),
    "D452" => Dict{String,Any}("2state" => true, "plat" => [0.80,0.85], "mdof" => 4),
    "D251" => Dict{String,Any}("2state" => true, "plat" => [0.50,0.60], "mdof" => 4),
    "E250" => Dict{String,Any}("2state" => true, "plat" => [0.30,0.50], "mdof" => 4),

    "H102" => Dict{String,Any}("2state" => false,"plat" => [0.20,0.80], "mdof" => 10),
    "N101" => Dict{String,Any}("2state" => false,"plat" => [0.25,0.75], "mdof" => 5),
    # "C101" => Dict{String,Any}("2state" => false, "plat" => , "mdof" => ),
    "C102" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.75], "mdof" => 5),
    "S400" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.75], "mdof" => 5),
    "N203" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "N200" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 20),
    "D200" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.65], "mdof" => 10),
    "D201" => Dict{String,Any}("2state" => false, "plat" => [0.30,0.60], "mdof" => 10),
    "N302" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "J306" => Dict{String,Any}("2state" => false, "plat" => [0.10,0.75], "mdof" => 80),
    "J303" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "J304" => Dict{String,Any}("2state" => false, "plat" => [0.15,0.80], "mdof" => 90),
    "E300" => Dict{String,Any}("2state" => false, "plat" => [0.10,0.70], "mdof" => 40),
    "F300" => Dict{String,Any}("2state" => false, "plat" => [0.15,0.70], "mdof" => 70),
    "J501" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 30),
)

# fPi_fitinfo = Dict{Bool,Dict{String,Dict{String,Any}}}(
#     "" => Dict{String,Any}("plat" => [,]),
# )


wpmm = Dict{String, Vector{Float64}}()
wpmm["D450"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["H101"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["H102r002"] = [5.0, -2.0, -1.0, -1.0]
wpmm["H400"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["N202"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N200"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N203"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N300"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["N302"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["J303"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J304"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["F300"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J306"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J307"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J500"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["E300"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["A654"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["E300"]     = [5.0, -2.0, -1.0, -1.0]

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> Data reading <==========================##

ens = "E250"

ens = EnsInfo(ens)

# use deflated when possible: 
path_rw = !isempty(filter(x->occursin(ens.id, x), readdir(joinpath(path_rw_,"reweight_deflated"), join=true))) ? joinpath(path_rw_,"reweight_deflated") : path_rw_

@info("Reading data")

println("   - Pion correlator")
corr_pi_pp  = get_mesons_corr(path_meson, ens, "ll", "PP", path_rw=path_rw, frw_bcwd=false, L=1 )
corr_pi_a0p = get_mesons_corr(path_meson, ens, "ll", "A0P", path_rw=path_rw, frw_bcwd=false, L=1 )

if ens.kappa_l != ens.kappa_s
    println("   - Kaon correlator")
    corr_k_pp  = get_mesons_corr(path_meson, ens, "ls", "PP", path_rw=path_rw, frw_bcwd=false, L=1 )
    corr_k_a0p = get_mesons_corr(path_meson, ens, "ls", "A0P", path_rw=path_rw, frw_bcwd=false, L=1 )
end


@info("Ready")

## pbc : 

t = collect(1:Int(length(corr_pi_pp.obs)/2 + 1))

obs_pp = -corr_pi_pp.obs[t]; uwerr.(obs_pp)
obs_a0p = corr_pi_a0p.obs[t]; uwerr.(obs_a0p)

fig = figure(figsize=(8,6))
errorbar(t, value.(obs_pp), err.(obs_pp), capsize=2, fmt="o", mfc="none", color="black", label = "Pi-PP 2pt. function")
errorbar(t, value.(obs_a0p), err.(obs_a0p), capsize=2, fmt="o", mfc="none", color="blue", label = "Pi-A0P 2pt. function")
# ylim(-0.02,0.05)
xlabel("t/a")
legend()
tight_layout()
display(gcf())
close()

## obc :

t = collect(1:Int(length(corr_pi_pp.obs)))

obs_pp  = -corr_pi_pp.obs[t]
obs_a0p = corr_pi_a0p.obs[t]

m_obs = meff(corr_pi_pp.obs); uwerr.(m_obs)
der_a0p = (obs_a0p[3:end] .- obs_a0p[1:end-2]) / 2 
der2_pp = (obs_pp[1:end-4] .- 2*obs_pp[3:end-2] .+ obs_pp[5:end])/4
der_a0p = der_a0p[2:end-1] .+ ca(ens.beta) .* der2_pp

aux = - der_a0p ./ ( 2. .* obs_pp[3:end-2]); uwerr.(aux)

mps = 0.18302
# 0.009202
m_pcac = 0.010988190841355997

rave = ZP(ens.beta)/ mps^2 * m_pcac

ba = 1 + 0.0472 * (6/ens.beta)
Za_l_sub(ens.beta) * (1 + ba * m_pcac  ) * rave

fig = figure(figsize=(8,6))
errorbar(t[2:end-2].+1/2, value.(m_obs), err.(m_obs), capsize=2, fmt="o", mfc="none", color="gray", label = "Pi eff. mass")
errorbar(t[3:end-2].+1/2, value.(aux), err.(aux), capsize=2, fmt="o", mfc="none", color="brown", label = "R ratio")

yscale("log")
# ylim(-0.02,0.05)
xlabel("t/a")
legend()
tight_layout()
display(gcf())
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> mPi & mK fit <==========================##

MESON = "Pion"  #  Pion  Kaon

# if MESON == "Pion"
#     STATE = mPi_fitinfo[ens.id]["2state"]
#     m_pl  = mPi_fitinfo[ens.id]["plat"]
#     f_pl  = fPi_fitinfo[ens.id]["plat"]
#     mdof  = mPi_fitinfo[ens.id]["mdof"]
# elseif MESON == "Kaon"
#     STATE = mK_fitinfo[ens.id]["2state"]
#     pl    = mK_fitinfo[ens.id]["plat"]
#     mdof  = mK_fitinfo[ens.id]["mdof"]
# end

STATE = false
m_pl = [0.30,0.50]
f_pl = [0.30,0.60]
mdof = 4

PVAL  = false

PLOT  = true
WRITE = false

if MESON == "Pion"
    corr_pp = deepcopy(corr_pi_pp)
    corr_a0p = deepcopy(corr_pi_a0p)
elseif MESON == "Kaon"
    corr_pp = deepcopy(corr_k_pp)
    corr_a0p = deepcopy(corr_k_a0p)
else
    error("Meson type $MESON not recognized.")
end

obs_pp = corr_pp.obs
obs_a0p = corr_a0p.obs
T = HVPobs.Data.get_T(corr_pp.id)

bc = ens.bc

AIC = true  # bc == "obc"


if bc == "obc"
    # m_obs = meff(obs[1:Int64(T/2+1)]) # first and last data points are lost (-2) and derivative is taken (-1); length(meff_) = length(obs) - 3
    m_obs = meff(obs_pp) # first and last data points are lost (-2) and derivative is taken (-1); length(meff_) = length(obs) - 3
    len = length(m_obs)

    der_a0p = (obs_a0p[3:end] .- obs_a0p[1:end-2]) / 2 
    der2_pp = (obs_pp[1:end-4] - 2*obs_pp[3:end-2] + obs_pp[5:end])/4
    der_a0p = der_a0p[2:end-1] + ca(ens.beta) * der2_pp

    aux = der_a0p ./ ( 2. .* obs_pp[3:end-2])

    @. const_model(x,p) = p[1] + 0*x

    # plconst_vec_left  = collect(floor(Int64,pl[1]*len):2:ceil(Int64,len/2)-ceil(Int64,mdof/2))
    # plconst_vec_right = collect(ceil(Int64,len/2)+floor(Int64,mdof/2):2:ceil(Int64,pl[2]*len))

    plconst_vec_left  = collect(floor(Int64,pl[1]*len):2:ceil(Int64,len*(pl[1]+pl[2])/2)-ceil(Int64,mdof/2))
    plconst_vec_right = collect(ceil(Int64,len*(pl[1]+pl[2])/2)+floor(Int64,mdof/2):2:ceil(Int64,pl[2]*len))

    if isempty(plconst_vec_left) || isempty(plconst_vec_right)
        error("Not enough data to fit with d.o.f. ≥ $mdof. Decrease $mdof or increase the plateau search range. ")
    end

    wpm = corr_pp.id in keys(wpmm) ? wpmm : nothing

    fitconst_vec = Vector{FitRes}()
    pl_vec = Vector{Vector{Float64}}()
    for p0 in plconst_vec_left
        for pf in plconst_vec_right
            # plateau = [p0,pl_f]
            m_data = m_obs[p0:pf] 
            fit = fit_routine(const_model,collect(p0:pf), m_data, 1, pval=PVAL, wpm=wpm, info=false, lineprint=false)
            push!(fitconst_vec,fit)
            push!(pl_vec,[p0,pf])
        end
    end

    fit_vec = fitconst_vec

elseif bc == "pbc"
    len = length(obs_pp)/2+1

    der_pp = (obs_pp[3:end] - obs_pp[1:end-2])/2.
    der_pp[1] = obs_pp[3] - obs_pp[2]
    push!(der_pp, obs_pp[end] - obs_pp[end-1])

    obs_a0p[2:end] .+= ca(ens.beta) .* der_pp

    der_a0p = (obs_a0p[3:end] .- obs_a0p[1:end-2]) / 2 
    der2_pp = (obs_pp[1:end-4] - 2*obs_pp[3:end-2] + obs_pp[5:end])/4
    der_a0p = der_a0p[2:end-1] + ca(ens.beta) * der2_pp

    aux = der_a0p ./ ( 2. .* obs_pp[3:end-2])


    obs_pp = -obs_pp[1:Int(len)]
    obs_a0p = obs_a0p[1:Int(len)]

    # cosh_model = 0.0
    if STATE
        @. cosh_model(x,p) = p[2] * (exp(-p[1]*x) + exp(-p[1]*(T-x))) +  p[4] * (exp(-p[3]*x) + exp(-p[3]*(T-x)))
        @. cosh_model_inv(x,p) = p[2] * (exp(-p[1]*x) - exp(-p[1]*(T-x))) +  p[4] * (exp(-p[3]*x) - exp(-p[3]*(T-x)))
        np = 4
    else
        @. cosh_model(x,p) = p[2] * (exp(-p[1]*x) + exp(-p[1]*(T-x)))
        @. cosh_model_inv(x,p) = p[2] * (exp(-p[1]*x) - exp(-p[1]*(T-x)))
        np = 2
    end

    plcosh_vec = collect(floor(Int64,m_pl[1]*len):1:ceil(Int64,min(len-mdof,m_pl[2]*len)))
    plcoshinv_vec = collect(floor(Int64,f_pl[1]*len):1:ceil(Int64,min(len-mdof,f_pl[2]*len)))
    pl_f = ceil(Int64,len)

    if isempty(plcosh_vec)
        error("Not enough data to fit with d.o.f. ≥ $mdof. Decrease $mdof or increase the plateau search range. ")
    end

    wpm = corr_pp.id in keys(wpmm) ? wpmm : nothing

    fitcosh_vec_pp = Vector{FitRes}()
    fitcosh_vec_a0p = Vector{FitRes}()
    for p0 in plcosh_vec
        data_pp = obs_pp[p0:end] 

        fit_pp = fit_routine(cosh_model,collect(p0:pl_f).-1, data_pp, np, pval=PVAL, wpm=wpm, info=false, lineprint=false)
        push!(fitcosh_vec_pp,fit_pp)
    end
    for p0 in plcoshinv_vec
        data_a0p = obs_a0p[p0:end]

        idx = findall(x->x==0.0, value.(data_a0p))
        [data_a0p[i] = uwreal([0.0, 1e-10], "artificialerr")  for i in idx] # artificial error is added so that the corr matrix doesn't explode

        fit_a0p = fit_routine(cosh_model_inv,collect(p0:pl_f).-1, data_a0p, np, pval=PVAL, wpm=wpm, info=false, lineprint=false)
        push!(fitcosh_vec_a0p,fit_a0p)
    end

    m_p0_vec  = plcosh_vec
    f_p0_vec  = plcoshinv_vec
    fit_vec_pp  = fitcosh_vec_pp
    fit_vec_a0p = fitcosh_vec_a0p
end

w_pp = get_w_from_fitres(vcat(fit_vec_pp...), AIC=AIC)
w_a0p = get_w_from_fitres(vcat(fit_vec_a0p...), AIC=AIC)

if bc ==  "pbc" && STATE
    m_res_vec = [[par[1],par[3]][argmin(value.([par[1],par[3]]))] for par in getfield.(vcat(fit_vec_pp...),:param)]
    cpp_res_vec  = [[par[2],par[4]][argmin(value.([par[1],par[3]]))] for par in getfield.(vcat(fit_vec_pp...),:param)]
    ca0p_res_vec = [[par[2],par[4]][argmin(value.([par[1],par[3]]))] for par in getfield.(vcat(fit_vec_a0p...),:param)]
else
    m_res_vec = [par[1] for par in getfield.(vcat(fit_vec_pp...),:param)]
    cpp_res_vec  = [par[2] for par in getfield.(vcat(fit_vec_pp...),:param)]
    ca0p_res_vec = [par[2] for par in getfield.(vcat(fit_vec_a0p...),:param)]
end

m_res, m_sys = model_average(m_res_vec, w_pp)
cpp_res, cpp_sys   = model_average(cpp_res_vec, w_pp)
ca0p_res, ca0p_sys = model_average(ca0p_res_vec, w_a0p)

m_res = m_res[1] + uwreal([0.0,m_sys],"meff MA syst"); uwerr(m_res)
cpp_res  = cpp_res[1] + uwreal([0.0,cpp_sys],"cpp MA syst"); uwerr(cpp_res)
ca0p_res = ca0p_res[1] + uwreal([0.0,ca0p_sys],"ca0p MA syst"); uwerr(ca0p_res)

f_res = Za_l_sub(ens.beta) * (sqrt(2) / sqrt(m_res) * ca0p_res / sqrt(cpp_res))

bestW = 20

##

if PLOT
    fig = figure(figsize=(16,12))
    # subplots_adjust(hspace=0.1)
    gs = fig.add_gridspec(4, 1, height_ratios=[4, 1, 1, 1])  # Adjust the height_ratios as needed

    if bc == "obc"
        ax1 = fig.add_subplot(gs[1, 1])
        title("$(corr_pp.id) ($MESON)")

        x0 = max(Int64(plconst_vec_left[1]-5),1):min(Int64(plconst_vec_right[end]+5),len)
        m_vec = m_obs[x0]; uwerr.(m_vec)

        errorbar(collect(x0).+1.5, value.(m_vec), err.(m_vec), fmt="o", capsize=2, color="black")
        fill_between([plconst_vec_left[1]-0.5,plconst_vec_right[end]+0.5].+1.5, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")
        
        WARG = sortperm(w,rev=true)[1:(length(w) >= bestW ? bestW : length(w))]
        uwerr.(res_vec[WARG])
        for (i,warg) in enumerate(WARG)
            y_penal = i/length(WARG)
            fill_between(pl_vec[warg].+1.5, value.(res_vec[warg])+y_penal*err.(res_vec[warg]), value.(res_vec[warg])-y_penal*err.(res_vec[warg]), alpha=w[warg], color="blue")
        end

        errorbar([plconst_vec_left[1],plconst_vec_left[end]].+1.5  , 0.5*(res.mean + maximum(value.(m_vec))).*[1,1], 0.01*(-res.mean + maximum(value.(m_vec))).*[1,1], fmt="", color="limegreen")
        errorbar([plconst_vec_right[1],plconst_vec_right[end]].+1.5, 0.5*(res.mean + maximum(value.(m_vec))).*[1,1], 0.01*(-res.mean + maximum(value.(m_vec))).*[1,1], fmt="", color="firebrick")

        axvline(x=plconst_vec_left[1]+1.0, color="cyan", linestyle="--")
        axvline(x=plconst_vec_right[end]+2.0, color="cyan", linestyle="--")

        axis("tight")
        xlabel(L"t/a")
        if MESON == "Pion"
            ylabel(L"$m^{\rm{eff}}_\pi(t)$")
        elseif MESON == "Kaon"
            ylabel(L"$m^{\rm{eff}}_K(t)$")
        end
        # ylim(res.mean-5*res.err,maximum(value.(m_vec) .+ 1.5 .* err.(m_vec)))

        ax2 = fig.add_subplot(gs[2, 1])
        # pl_vec[WARG]
        errorbar(collect(range(x0[1],x0[end],length(WARG))), value.(res_vec[WARG]), err.(res_vec[WARG]), fmt="d", mfc="none", color="blue")
        fill_between([x0[1],x0[end]], value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")

        if MESON == "Pion"
            ylabel(L"$m_\pi$")
        elseif MESON == "Kaon"
            ylabel(L"$m_K$")
        end
        setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels

        ax3 = fig.add_subplot(gs[3, 1])

        PyPlot.plot(collect(range(x0[1],x0[end],length(WARG))), w[WARG], linestyle="none", marker="o", mfc="none", color="blue")
        if AIC
            ylabel(latexstring("\\rm{w}\\ [\\rm{AIC}]"))
        else
            ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))
        end

        if PVAL
            setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels

            ax4 = fig.add_subplot(gs[4, 1])

            pval_vec = getfield.(fit_vec,:pval)

            PyPlot.plot(collect(range(x0[1],x0[end],length(WARG))), pval_vec[WARG], linestyle="none", marker="o", mfc="none", color="blue")

            ylabel(L"$\rm{p-values}$")
            ax4.set_xticks(collect(range(x0[1],x0[end],length(WARG))))
            ax4.set_xticklabels(string.([Int.(pl) for pl in pl_vec[WARG]]))
        else
            ax3.set_xticks(collect(range(x0[1],x0[end],length(WARG))))
            ax3.set_xticklabels(string.([Int.(pl) for pl in pl_vec[WARG]]))
        end
        xlabel("Best fits")


        tight_layout()
        display(fig)
        close("all")

    elseif bc == "pbc"
        for curr in ["A0P"]

            if curr == "PP"
                pl = m_pl
                obs = obs_pp
                w = w_pp
                p0_vec = m_p0_vec
                res_vec = m_res_vec
                res = m_res
                fit_vec = fit_vec_pp
            elseif curr == "A0P"
                pl = f_pl
                obs = obs_a0p
                w = w_a0p
                p0_vec = f_p0_vec
                res_vec = ca0p_res_vec
                res = ca0p_res
                fit_vec = fit_vec_a0p
            end

            ax1 = fig.add_subplot(gs[1, 1])
            title("$(corr_pp.id) ($MESON)")

            x0 = collect(max(floor(Int64,len*(pl[1])),1)-2:len) .- 1
            vec = obs[Int(x0[1])+1:end]; uwerr.(vec)

            errorbar(x0, value.(vec), err.(vec), fmt="o", capsize=2, color="black")
            maxw_arg = argmax(w[1:length(fit_vec)])
            par = getfield.(fit_vec,:param)[maxw_arg]

            x_ = collect(max(floor(Int64,len*(pl[1])),1)-0.5:0.1:len).-1
            if curr == "PP"
                ycosh_fit = cosh_model(x_,par); uwerr.(ycosh_fit)
            elseif curr == "A0P"
                ycosh_fit = cosh_model_inv(x_,par); uwerr.(ycosh_fit)
            end

            fill_between(x_, value.(ycosh_fit)+err.(ycosh_fit), value.(ycosh_fit)-err.(ycosh_fit), alpha=0.3, color="orange")

            axvline(x=p0_vec[1]-1.5, color="red", linestyle="--")
            axvline(x=p0_vec[end]-0.5, color="red", linestyle="--")

            axis("tight")
            ylabel(L"$C(t)$")

            p0_vec = Float64.(p0_vec) .- (1.5 + 1.0)
            
            if curr == "PP"
                ax21 = fig.add_subplot(gs[2, 1])
                uwerr.(res_vec)

                errorbar(p0_vec .+ 1.5, value.(res_vec), err.(res_vec), fmt="d", mfc="none", color="blue")
                fill_between(x0, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")

                if MESON == "Pion"
                    ylabel(L"$m^{\rm{PP}}_\pi$")
                elseif MESON == "Kaon"
                    ylabel(L"$m^{\rm{PP}}_K$")
                end
                setp(ax21.get_xticklabels(),visible=false) # Disable x tick labels

                ax22 = fig.add_subplot(gs[3, 1])
                uwerr.(cpp_res_vec)

                errorbar(p0_vec .+ 1.5, value.(cpp_res_vec), err.(cpp_res_vec), fmt="d", mfc="none", color="blue")
                fill_between(x0, value.(cpp_res)+err.(cpp_res), value.(cpp_res)-err.(cpp_res), alpha=0.4, color="green")

                if MESON == "Pion"
                    ylabel(L"$C^{\rm{PP}}_\pi$")
                elseif MESON == "Kaon"
                    ylabel(L"$C^{\rm{PP}}_K$")
                end
                setp(ax22.get_xticklabels(),visible=false) # Disable x tick labels

                pl_idx = 3
            elseif curr == "A0P"
                ax2 = fig.add_subplot(gs[2, 1])
                uwerr.(res_vec)

                errorbar(p0_vec .+ 1.5, value.(res_vec), err.(res_vec), fmt="d", mfc="none", color="blue")
                fill_between(x0, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")

                # ylabel(L"m_{D_s}")
                if MESON == "Pion"
                    ylabel(L"$C^{\rm{A0P}}_\pi$")
                elseif MESON == "Kaon"
                    ylabel(L"$C^{\rm{A0P}}_K$")
                end
                setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels

                pl_idx = 2
            end

            ax3 = fig.add_subplot(gs[pl_idx+1, 1])

            fill_between(x0, maximum(w)/2, maximum(w)/2, alpha=0.0, color="white")
            PyPlot.plot(p0_vec .+ 1.5, w, linestyle="none", marker="o", mfc="none", color="blue")

            ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))

            if PVAL
                setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels

                ax4 = fig.add_subplot(gs[pl_idx+2, 1])

                pval_vec = getfield.(fit_vec,:pval)
                fill_between(x0, maximum(vcat(pval_vec...))/2, maximum(vcat(pval_vec...))/2, alpha=0.0, color="white")
                PyPlot.plot(p0_vec .+ 1.5, pval_vec, linestyle="none", marker="o", mfc="none", color="blue")

                ylabel(L"$\rm{p-values}$")
            end
            xlabel(L"$t/a$")


            tight_layout()
            display(fig)
            close()
        end
    end
end
if WRITE


end

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> fPi fit <==========================##

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> Reading test <==========================##

if PLOT
    fig = figure(figsize=(16,12))
    # subplots_adjust(hspace=0.1)
    gs = fig.add_gridspec(4, 1, height_ratios=[4, 1, 1, 1])  # Adjust the height_ratios as needed

    if bc == "obc"
        ax1 = fig.add_subplot(gs[1, 1])
        title("$(corr_pp.id) ($MESON)")

        x0 = max(Int64(plconst_vec_left[1]-5),1):min(Int64(plconst_vec_right[end]+5),len)
        m_vec = m_obs[x0]; uwerr.(m_vec)

        errorbar(collect(x0).+1.5, value.(m_vec), err.(m_vec), fmt="o", capsize=2, color="black")
        fill_between([plconst_vec_left[1]-0.5,plconst_vec_right[end]+0.5].+1.5, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")
        
        WARG = sortperm(w,rev=true)[1:(length(w) >= bestW ? bestW : length(w))]
        uwerr.(res_vec[WARG])
        for (i,warg) in enumerate(WARG)
            y_penal = i/length(WARG)
            fill_between(pl_vec[warg].+1.5, value.(res_vec[warg])+y_penal*err.(res_vec[warg]), value.(res_vec[warg])-y_penal*err.(res_vec[warg]), alpha=w[warg], color="blue")
        end

        errorbar([plconst_vec_left[1],plconst_vec_left[end]].+1.5  , 0.5*(res.mean + maximum(value.(m_vec))).*[1,1], 0.01*(-res.mean + maximum(value.(m_vec))).*[1,1], fmt="", color="limegreen")
        errorbar([plconst_vec_right[1],plconst_vec_right[end]].+1.5, 0.5*(res.mean + maximum(value.(m_vec))).*[1,1], 0.01*(-res.mean + maximum(value.(m_vec))).*[1,1], fmt="", color="firebrick")

        axvline(x=plconst_vec_left[1]+1.0, color="cyan", linestyle="--")
        axvline(x=plconst_vec_right[end]+2.0, color="cyan", linestyle="--")

        axis("tight")
        xlabel(L"t/a")
        if MESON == "Pion"
            ylabel(L"$m^{\rm{eff}}_\pi(t)$")
        elseif MESON == "Kaon"
            ylabel(L"$m^{\rm{eff}}_K(t)$")
        end
        # ylim(res.mean-5*res.err,maximum(value.(m_vec) .+ 1.5 .* err.(m_vec)))

        ax2 = fig.add_subplot(gs[2, 1])
        # pl_vec[WARG]
        errorbar(collect(range(x0[1],x0[end],length(WARG))), value.(res_vec[WARG]), err.(res_vec[WARG]), fmt="d", mfc="none", color="blue")
        fill_between([x0[1],x0[end]], value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")

        if MESON == "Pion"
            ylabel(L"$m_\pi$")
        elseif MESON == "Kaon"
            ylabel(L"$m_K$")
        end
        setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels

        ax3 = fig.add_subplot(gs[3, 1])

        PyPlot.plot(collect(range(x0[1],x0[end],length(WARG))), w[WARG], linestyle="none", marker="o", mfc="none", color="blue")
        if AIC
            ylabel(latexstring("\\rm{w}\\ [\\rm{AIC}]"))
        else
            ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))
        end

        if PVAL
            setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels

            ax4 = fig.add_subplot(gs[4, 1])

            pval_vec = getfield.(fit_vec,:pval)

            PyPlot.plot(collect(range(x0[1],x0[end],length(WARG))), pval_vec[WARG], linestyle="none", marker="o", mfc="none", color="blue")

            ylabel(L"$\rm{p-values}$")
            ax4.set_xticks(collect(range(x0[1],x0[end],length(WARG))))
            ax4.set_xticklabels(string.([Int.(pl) for pl in pl_vec[WARG]]))
        else
            ax3.set_xticks(collect(range(x0[1],x0[end],length(WARG))))
            ax3.set_xticklabels(string.([Int.(pl) for pl in pl_vec[WARG]]))
        end
        xlabel("Best fits")


        tight_layout()
        display(fig)
        close("all")

    elseif bc == "pbc"
        ax1 = fig.add_subplot(gs[1, 1])
        title("$(corr_pp.id) ($MESON)")

        m_x0 = collect(max(floor(Int64,len*(m_pl[1])),1)-2:len) .- 1
        m_vec = obs_pp[Int(m_x0[1])+1:end]; uwerr.(m_vec)
        f_x0 = collect(max(floor(Int64,len*(f_pl[1])),1)-2:len) .- 1
        f_vec = obs_a0p[Int(f_x0[1])+1:end]; uwerr.(f_vec)


        errorbar(m_x0, value.(m_vec), err.(m_vec), fmt="o", capsize=2, color="black", label="corr. PP")
        errorbar(f_x0, value.(f_vec), err.(f_vec), fmt="o", capsize=2, color="gray", label="corr. A0P")

        maxw_arg_pp = argmax(w_pp[1:length(fit_vec_pp)])
        par_pp = getfield.(fit_vec_pp,:param)[maxw_arg_pp]
        maxw_arg_a0p = argmax(w_a0p[1:length(fit_vec_a0p)])
        par_a0p = getfield.(fit_vec_a0p,:param)[maxw_arg_a0p]

        m_x = collect(max(floor(Int64,len*(m_pl[1])),1)-0.5:0.1:len).-1
        f_x = collect(max(floor(Int64,len*(f_pl[1])),1)-0.5:0.1:len).-1
        ycosh_fit = cosh_model(m_x,par_pp); uwerr.(ycosh_fit)
        ycoshinv_fit = cosh_model_inv(f_x,par_a0p); uwerr.(ycoshinv_fit)
        fill_between(m_x, value.(ycosh_fit)+err.(ycosh_fit), value.(ycosh_fit)-err.(ycosh_fit), alpha=0.3, color="orange")
        fill_between(f_x, value.(ycoshinv_fit)+err.(ycoshinv_fit), value.(ycoshinv_fit)-err.(ycoshinv_fit), alpha=0.3, color="red")

        axvline(x=m_p0_vec[1]-1.5, color="orange", linestyle="--")
        axvline(x=m_p0_vec[end]-0.5, color="orange", linestyle="--")
        axvline(x=f_p0_vec[1]-1.5, color="red", linestyle="--")
        axvline(x=f_p0_vec[end]-0.5, color="red", linestyle="--")

        axis("tight")
        ylabel(L"$C(t)$")

        m_p0_vec = Float64.(m_p0_vec) .- (1.5 + 1.0)
        f_p0_vec = Float64.(f_p0_vec) .- (1.5 + 1.0)
        
        ax2 = fig.add_subplot(gs[2, 1])
        uwerr.(m_res_vec)

        errorbar(m_p0_vec .+ 1.5, value.(m_res_vec), err.(m_res_vec), fmt="d", mfc="none", color="blue")
        fill_between(m_x0, value.(m_res)+err.(m_res), value.(m_res)-err.(m_res), alpha=0.4, color="green")

        # ylabel(L"m_{D_s}")
        if MESON == "Pion"
            ylabel(L"$m_\pi$")
        elseif MESON == "Kaon"
            ylabel(L"$m_K$")
        end
        setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels

        ax3 = fig.add_subplot(gs[3, 1])

        fill_between(m_x0, maximum(w_pp)/2, maximum(w_pp)/2, alpha=0.0, color="white")
        PyPlot.plot(m_p0_vec .+ 1.5, w_pp, linestyle="none", marker="o", mfc="none", color="blue")

        ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))

        if PVAL
            setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels

            ax4 = fig.add_subplot(gs[4, 1])

            pval_vec = getfield.(fit_vec,:pval)
            fill_between(x0, maximum(vcat(pval_vec...))/2, maximum(vcat(pval_vec...))/2, alpha=0.0, color="white")
            PyPlot.plot(p0_vec .+ 1.5, pval_vec, linestyle="none", marker="o", mfc="none", color="blue")

            ylabel(L"$\rm{p-values}$")
        end
        xlabel(L"$t/a$")


        tight_layout()
        display(fig)
        close("all")

    end
end