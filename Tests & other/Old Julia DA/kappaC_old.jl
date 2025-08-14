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

path_heavy = joinpath(julia_script_directory, "..", "HVPData", "s_heavy_data")

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

# using OrderedCollections

# include("HVPtools/Const.jl")

# include("HVPtools/Utils.jl")

# include("HVPtools/Fit&MA.jl")

# path to ObsBDIO

path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")

# Ensamble choice

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302


# Already included: "H101","H102","N101","C101","B450","S400","N202","N203","N200","D200","E250","N300","N302","J303","E300","A654","D452","J500","","","","","","","","","",""
# Without charm:    "J501","N451","D150","D451","J304","C102","D251","D201","J306","J307","F300","H200"

# To include: ""

# physical MDs mass

ensList = ["A653"] # already included
ensInfo = EnsInfo.(ensList)

# filterInfo = (getfield.(ensInfo,:beta) .!= 3.34) .& (getfield.(ensInfo,:beta) .!= 3.85)
# ensInfo = ensInfo[filterInfo]

# @info("No ensembles with ß = 3.34 or 3.85 can be considered in this analysis. The ensembles considered are: \n - $(getfield.(ensInfo,:id))")


# dict struct.: AIC >> ens >> fit kwargs
mDs_fitinfo = Dict{Bool,Dict{String,Dict{String,Any}}}(
    false => Dict{String,Dict{String,Any}}(
        "H101" => Dict{String,Any}("const" => [0.50,0.54], "2state" => [0.33,0.37], "plf" => 1.00, "mdof" => 4),
        "H102" => Dict{String,Any}("const" => [0.60,0.70], "2state" => [0.45,0.50], "plf" => 1.00, "mdof" => 4),
        "N101" => Dict{String,Any}("const" => [0.62,0.75], "2state" => [0.30,0.35], "plf" => 1.00, "mdof" => 4),
        "C101" => Dict{String,Any}("const" => [0.72,0.78], "2state" => [0.45,0.50], "plf" => 1.00, "mdof" => 4),
        "B450" => Dict{String,Any}("const" => [0.77,1.00], "2state" => [0.40,0.42], "plf" => 0.85, "mdof" => 2),
        "S400" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.22,0.30], "plf" => 1.00, "mdof" => 4),
        "N202" => Dict{String,Any}("const" => [0.50,0.60], "2state" => [0.20,0.30], "plf" => 1.00, "mdof" => 4),
        "N203" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.20,0.30], "plf" => 1.00, "mdof" => 4),
        "N200" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.25,0.30], "plf" => 1.00, "mdof" => 4),
        "D200" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.24,0.30], "plf" => 1.00, "mdof" => 4),
        "E250" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.25,0.30], "plf" => 0.90, "mdof" => 4),
        "N300" => Dict{String,Any}("const" => [0.55,0.60], "2state" => [0.17,0.25], "plf" => 1.00, "mdof" => 4),
        "N302" => Dict{String,Any}("const" => [0.60,0.70], "2state" => [0.20,0.30], "plf" => 1.00, "mdof" => 4),
        "J303" => Dict{String,Any}("const" => [0.55,0.60], "2state" => [0.13,0.17], "plf" => 1.00, "mdof" => 4),
    ),
    true => Dict{String,Dict{String,Any}}(
        "H101" => Dict{String,Any}("const" => [0.50,0.54], "2state" => [0.40,0.43], "plf" => 1.00, "mdof" => 4),
        "H102" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.24,0.40], "plf" => 1.00, "mdof" => 4),
        "N101" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.30,0.35], "plf" => 1.00, "mdof" => 4),
        "C101" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.27,0.40], "plf" => 1.00, "mdof" => 4),
        "B450" => Dict{String,Any}("const" => [0.75,1.00], "2state" => [0.40,0.43], "plf" => 0.85, "mdof" => 2),
        "S400" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.22,0.30], "plf" => 1.00, "mdof" => 4),
        "N202" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.15,0.30], "plf" => 1.00, "mdof" => 4),
        "N203" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.20,0.30], "plf" => 1.00, "mdof" => 4),
        "N200" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.17,0.30], "plf" => 1.00, "mdof" => 4),
        "D200" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.24,0.30], "plf" => 0.85, "mdof" => 4),
        "E250" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.27,0.30], "plf" => 0.89, "mdof" => 4),
        # "N300" => Dict{String,Any}("const" => [0.60,0.75], "2state" => [0.15,0.20], "plf" => 1.00, "mdof" => 4), # paper Kc(tar)
        "N300" => Dict{String,Any}("const" => [0.55,0.60], "2state" => [0.18,0.23], "plf" => 1.00, "mdof" => 4),
        "N302" => Dict{String,Any}("const" => [0.58,0.60], "2state" => [0.23,0.27], "plf" => 1.00, "mdof" => 4),
        "J303" => Dict{String,Any}("const" => [0.40,0.45], "2state" => [0.16,0.20], "plf" => 1.00, "mdof" => 4),
        "E300" => Dict{String,Any}("const" => [0.40,0.50], "2state" => [0.15,0.25], "plf" => 1.00, "mdof" => 4),
        "J500" => Dict{String,Any}("const" => [0.45,0.50], "2state" => [0.10,0.20], "plf" => 1.00, "mdof" => 4), # beta = 3.85
        "A654" => Dict{String,Any}("const" => [0.70,0.80], "2state" => [0.40,0.50], "plf" => 0.85, "mdof" => 2), # beta = 3.34
        "D452" => Dict{String,Any}("const" => [0.45,0.55], "2state" => [0.20,0.30], "plf" => 0.90, "mdof" => 4),
        "D450" => Dict{String,Any}("const" => [0.45,0.50], "2state" => [0.20,0.30], "plf" => 0.90, "mdof" => 4),
        "A653" => Dict{String,Any}("const" => [0.80,0.85], "2state" => [0.45,0.55], "plf" => 0.85, "mdof" => 2), # beta = 3.34

    )
)


## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> mDs fit ensemble per ensemble <==========================##

# Set important parameters

NOERR_MDs = true
PLOT_meff = true
PLOT_Ktar = true

AIC = true

# Find MDs prime values in lattice units as a function of beta

MD_ph

MD_ph_prime = MD_ph * (sqrtt0_ph_Bruno / sqrtt0_ph)

sqrtt0_ph_Bruno

uwerr(MD_ph_prime); println("- MDs shifts to MDs' = $MD_ph_prime GeV")

aMDSU3_prime = Dict{String, uwreal}()
for beta in b_values
    aMDSU3_prime["$beta"] = MD_ph_prime * (sqrtt0_ph / sqrt(t0sym(beta,Bruno=true))) / hbarc
end
aMDSU3_prime
[uwerr(aMDSU3_prime[key]) for key in keys(aMDSU3_prime)]

if NOERR_MDs
    [set_fluc_to_zero!(aMDSU3_prime[key], "sqrtt0 [fm] (Bruno)") for key in keys(aMDSU3_prime)]
    [set_fluc_to_zero!(aMDSU3_prime[key], "t0sym/a2") for key in keys(aMDSU3_prime)]
    [aMDSU3_prime[key] *= 1.0 for key in keys(aMDSU3_prime)] 
end

##-- Ensemble choice and plateau parameters

ens_id = "E250"

# pl2state0 = [0.45,0.55]
# plconst0  = [0.80,0.85]
# plf       = 0.85
# mdof      = 2

pl2state0 = mDs_fitinfo[AIC][ens_id]["2state"]
plconst0  = mDs_fitinfo[AIC][ens_id]["const"]
plf       = mDs_fitinfo[AIC][ens_id]["plf"]
mdof      = mDs_fitinfo[AIC][ens_id]["mdof"]

# Read s-heavy data and find kappaC target using MDs prime

ens = EnsInfo(ens_id)

kappaC_sim = get_kappa_values()

amDs   = Dict{String,Vector{uwreal}}()
kappaC_tar = Dict{String,uwreal}()

println("- Ens: $(ens.id)")

data_sheavy = read_kappa_charm_all_config(joinpath(path_heavy,ens.id))

println("   - Computing amDs in simulations...")

corr = corr_obs(data_sheavy["sh1"])
corr.obs
amDs[ens.id] = Vector{uwreal}()

for sh in ["sh1","sh2","sh3","sh4"]
    println("      - strange-heavy $(sh[end])")

    corr = corr_obs(data_sheavy[sh])

    amDs_res, amDs_syst = meff_MA(corr; 
        pl2state0=pl2state0, plconst0=plconst0, plf=plf, 
        plstep=1, mdof=mdof, 
        state_fit=true, AIC=AIC, 
        returnfitMA=false, fitinfo=false,
        plot=PLOT_meff, pval=PLOT_meff
    )

    push!(amDs[ens.id], amDs_res + uwreal([0.0,amDs_syst],"amDs syst"))
end

println("   - Computing KappaC target...")

kappaCinv_sim = 1 ./ kappaC_sim[ens.id][2:end]

@. lin_model(x,p) = p[1] + p[2] * x

fit = fit_routine(lin_model, kappaCinv_sim, amDs[ens.id], 2, lineprint=false)
par = fit.param

kappaCinv_tar = (aMDSU3_prime["$(ens.beta)"] - par[1]) / par[2]
kappaC_tar[ens.id] = 1/kappaCinv_tar

if PLOT_Ktar
    println("   - Ploting")

    kappaC_arr = range(kappaCinv_sim[1],kappaCinv_sim[end],100)
    mDs_arr = lin_model(kappaC_arr,par); uwerr.(mDs_arr)
    uwerr.(amDs[ens.id]); uwerr(aMDSU3_prime["$(ens.beta)"]); uwerr(kappaCinv_tar); uwerr(kappaC_tar[ens.id])
    
    title("Ens: $(ens.id); "*L"$\kappa_c$"*" = $(round(value(kappaC_tar[ens.id]),digits=6)) ± $(round(err(kappaC_tar[ens.id]),digits=6))")
    errorbar(kappaCinv_sim, value.(amDs[ens.id]), err.(amDs[ens.id]), fmt="o", mfc="none", capsize=2, color="blue")
    errorbar(value(kappaCinv_tar), value(aMDSU3_prime["$(ens.beta)"]), xerr = err(kappaCinv_tar), yerr = 0.0, fmt="^", mfc="none", capsize=2, color="red")
    fill_between(kappaC_arr, value.(mDs_arr)+err.(mDs_arr), value.(mDs_arr)-err.(mDs_arr), alpha=0.4, color="blue")
    axis("tight")
    ax = gca()      # get the handle of the current axis (not really used here)
    xlabel(L"$1/\kappa_c$")
    ylabel(L"$am_{D_s}$")
    display(gcf())      #display the figure
    close()
end



## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> Loop for all ensembles (+ BDIO print) <==========================##

# Set important parameters

NOERR_MDs = true
PLOT_meff = true
PLOT_Ktar = true

AIC = true

PRINTmDS  = false
PRINTBDIO = true
OVERWRITE = false

# Find MDs prime values in lattice units as a function of beta

println("- Computing mDs prime in the physical points and in the sym points...")

MD_ph_prime = MD_ph * (sqrtt0_ph_Bruno / sqrtt0_ph)  # shift in the mass, usefull to compute the derivative

uwerr(MD_ph_prime); println("- MDs shifts to MDs' = $MD_ph_prime GeV")

aMDSU3_prime = Dict{String, uwreal}()
for beta in b_values
    aMDSU3_prime["$beta"] = MD_ph_prime * (sqrtt0_ph / sqrt(t0sym(beta,Bruno=true))) / hbarc
end
aMDSU3_prime
[uwerr(aMDSU3_prime[key]) for key in keys(aMDSU3_prime)]

if NOERR_MDs
    [set_fluc_to_zero!(aMDSU3_prime[key], "sqrtt0 [fm] (Bruno)") for key in keys(aMDSU3_prime)]
    [set_fluc_to_zero!(aMDSU3_prime[key], "t0sym/a2") for key in keys(aMDSU3_prime)]
    [aMDSU3_prime[key] *= 1.0 for key in keys(aMDSU3_prime)] 
end

if PRINTmDS
    pBDIOkappaCtar = joinpath(path_bdio,"kappaC_tar","mDs_prime")
    if ispath(pBDIOkappaCtar)
        OVERWRITE ? rm(pBDIOkappaCtar, recursive=true) : error("This information already exist; set 'OVERWRITE = true' to overwrite")
    end

    io = IOBuffer()
    write(io, "MDs 'prime' physical and sym points")
    fb = ALPHAdobs_create(pBDIOkappaCtar, io)

    ALPHAdobs_write(fb, MD_ph_prime)
    ALPHAdobs_write(fb, aMDSU3_prime)

    ALPHAdobs_close(fb)
end

# Read s-heavy data and find kappaC target using MDs prime

println("- Reading s-heavy data and interpolating for kappaC target...")

kappaC_sim = get_kappa_values()

amDs   = Dict{String,Array{uwreal}}()
kappaC_tar = Dict{String,uwreal}()
for ens in ensInfo
    println("   - Ens: $(ens.id)")

    data_sheavy = read_kappa_charm_all_config(joinpath(path_heavy,ens.id))

    println("      - Computing amDs in simulations...")

    amDs[ens.id] = []
    for sh in ["sh1","sh2","sh3","sh4"]
        println("         - strange-heavy $(sh[end])")

        corr = corr_obs(data_sheavy[sh])

        amDs_res, amDs_syst = meff_MA(corr; 
            pl2state0=mDs_fitinfo[AIC][ens.id]["2state"], plconst0=mDs_fitinfo[AIC][ens.id]["const"], plf=mDs_fitinfo[AIC][ens.id]["plf"], 
            plstep=1, mdof=mDs_fitinfo[AIC][ens.id]["mdof"], 
            state_fit=true, AIC=AIC, 
            returnfitMA=false, # fitinfo=false,
            plot=PLOT_meff, pval=PLOT_meff
        )

        push!(amDs[ens.id], amDs_res + uwreal([0.0,amDs_syst],"amDs syst"))
    end

    println("      - Computing KappaC target...")

    kappaCinv_sim = 1 ./ kappaC_sim[ens.id][2:end]
    
    @. lin_model(x,p) = p[1] + p[2] * x
    
    fit = fit_routine(lin_model, kappaCinv_sim, amDs[ens.id], 2) #, showinfo=false)
    par = fit.param
    
    kappaCinv_tar = (aMDSU3_prime["$(ens.beta)"] - par[1]) / par[2]
    kappaC_tar[ens.id] = 1/kappaCinv_tar

    uwerr(kappaC_tar[ens.id]); println("         ⟹ kappaC (tar) = $(kappaC_tar[ens.id])")
    
    if PLOT_Ktar
        println("      - Ploting...")

        kappaC_arr = range(kappaCinv_sim[1],kappaCinv_sim[end],100)
        mDs_arr = lin_model(kappaC_arr,par); uwerr.(mDs_arr)
        uwerr.(amDs[ens.id]); uwerr(aMDSU3_prime["$(ens.beta)"]); uwerr(kappaCinv_tar); uwerr(kappaC_tar[ens.id])
        
        title("Ens: $(ens.id); "*L"$\kappa_c$"*" = $(round(value(kappaC_tar[ens.id]),digits=6)) ± $(round(err(kappaC_tar[ens.id]),digits=6))")
        errorbar(kappaCinv_sim, value.(amDs[ens.id]), err.(amDs[ens.id]), fmt="o", mfc="none", capsize=2, color="blue")
        errorbar(value(kappaCinv_tar), value(aMDSU3_prime["$(ens.beta)"]), xerr = err(kappaCinv_tar), yerr = 0.0, fmt="^", mfc="none", capsize=2, color="red")
        fill_between(kappaC_arr, value.(mDs_arr)+err.(mDs_arr), value.(mDs_arr)-err.(mDs_arr), alpha=0.4, color="blue")
        axis("tight")
        ax = gca()      # get the handle of the current axis (not really used here)
        xlabel(L"$1/\kappa_c$")
        ylabel(L"$am_{D_s}$")
        display(gcf())      #display the figure
        close()
    end

    if PRINTBDIO
        pens = joinpath(path_bdio,"kappaC_tar",ens.id)
        !ispath(pens) ? mkdir(pens) : nothing

        pBDIO_mDs = joinpath(pens,"$(ens.id)_mDs")
        pBDIO_kappaC = joinpath(pens,"$(ens.id)_kappaC")
        if ispath(pBDIO_mDs) || ispath(pBDIO_kappaC)
            if OVERWRITE
                rm(pBDIO_mDs, recursive=true)
                rm(pBDIO_kappaC, recursive=true)
            else
                error("This information already exist; set 'OVERWRITE = true' to overwrite")
            end
        end

        io = IOBuffer()
        write(io, "$(ens.id) mDs 'prime'")
        fb = ALPHAdobs_create(pBDIO_mDs, io)

        extra = Dict{String, Any}("ens" => ens.id)
        ALPHAdobs_write(fb, amDs[ens.id], extra=extra)
        
        ALPHAdobs_close(fb)


        io = IOBuffer()
        write(io, "$(ens.id) kappaC target")
        fb = ALPHAdobs_create(pBDIO_kappaC, io)

        extra = Dict{String, Any}("ens" => ens.id)
        ALPHAdobs_write(fb, kappaC_tar[ens.id], extra=extra)
        
        ALPHAdobs_close(fb)
    end
end

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> Reading test <==========================##

ensid = "H102"
data  = "kappaC tar"  # mDs prime   kappaC tar   mDs sim


if data == "mDs prime"
    fb = BDIO_open(joinpath(path_bdio,"kappaC_tar","mDs_prime"),"r")
    res = Dict{String,Any}()
    i = 0
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        nobs = d["nobs"]
        dims = d["dimensions"]
        i=i+1
        if i==1
            res["mDs ph"] = ALPHAdobs_read_next(fb)
        else
            ks = collect(d["keys"])
            res["amDs SU3"] = ALPHAdobs_read_next(fb, keys=ks)
        end
    end
elseif data == "kappaC tar"
    fb = BDIO_open(joinpath(path_bdio,"kappaC_tar",ensid,"$(ensid)_kappaC"),"r")
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        nobs = d["nobs"]
        dims = d["dimensions"]
        res = ALPHAdobs_read_next(fb)
    end
elseif data == "mDs sim"
    fb = BDIO_open(joinpath(path_bdio,"kappaC_tar",ensid,"$(ensid)_mDs"),"r")
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        nobs = d["nobs"]
        dims = d["dimensions"]
        sz = tuple(d["size"]...)
        res = ALPHAdobs_read_next(fb, size=sz)
    end
end

uwerr.(res); res
