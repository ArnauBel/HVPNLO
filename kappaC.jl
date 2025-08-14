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

# include uwreal constants

include("HVPtool/uwConst.jl")

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

# Path definition

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

# Ensamble choice

# Already included: "H101","H102","N101","C101","B450","S400","N202","N203","N200","D200","E250","N300","N302","J303","E300","A654","D452","J500","","","","","","","","","",""
# Without charm:    "J501","N451","D150","D451","J304","C102","D251","D201","J306","J307","F300", ¿"H200"?

# To include: ""

# physical MDs mass

ensList = ["A653","A654","B450","C101","D200","D450","D452","E250","E300","H101","H102","J303","J500","N101","N200","N202","N203","N300","N302","S400"] 
ensInfo = EnsInfo.(ensList)

# filterInfo = (getfield.(ensInfo,:beta) .!= 3.34) .& (getfield.(ensInfo,:beta) .!= 3.85)
# ensInfo = ensInfo[filterInfo]

# @info("No ensembles with ß = 3.34 or 3.85 can be considered in this analysis. The ensembles considered are: \n - $(getfield.(ensInfo,:id))")

wpmm = Dict{String, Vector{Float64}}()
wpmm["A653"] = [5.0, -2.0, -1.0, -1.0]
wpmm["A654"] = [5.0, -2.0, -1.0, -1.0]
wpmm["D450"] = [5.0, -2.0, -1.0, -1.0]
wpmm["E300"] = [5.0, -2.0, -1.0, -1.0]
wpmm["J303"] = [5.0, -2.0, -1.0, -1.0]
wpmm["J500"] = [5.0, -2.0, -1.0, -1.0]


@info("Ready")

## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##

##==========================> Data reading and testing <==========================##

ensid = "D200"
sh    = "sh1"

ens   = EnsInfo(ensid)

@info("Reading data")

println("   - Pion correlator")

data_sheavy = read_kappa_charm_all_config(joinpath(path_heavy,ensid))
corr = corr_obs(data_sheavy[sh])

@info("Ready")

## corr :

obs = corr.obs; uwerr.(obs)

t = collect(1:length(obs))

fig = figure(figsize=(8,6))
errorbar(t, value.(obs), err.(obs), capsize=2, fmt="o", mfc="none", color="black", label = "s-heavy 2pt. function")
yscale("log")
xlabel("t/a")
legend()
tight_layout()
display(gcf())
close()

## effective mass

obs = corr.obs

if ens.bc == "obc"
    m_obs = meff(obs[1:end-2]); uwerr.(m_obs)
elseif ens.bc == "pbc"
    m_obs = meff(obs[1:Int(HVPobs.Data.get_T(ens.id)/2 + 1)]); uwerr.(m_obs)
end

t = collect(2.:length(m_obs)+1.)

fig = figure(figsize=(8,6))
errorbar(t.+1/2, value.(m_obs), err.(m_obs), capsize=2, fmt="o", mfc="none", color="gray", label = "Ds eff. mass")

yscale("log")
# ylim(-0.02,0.05)
xlabel("t/a")
legend()
tight_layout()
display(gcf())
close()

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> Loop for all ensembles (+ BDIO print) <==========================##

# Set important parameters

ensid = ""

mDs_fitinfo = Dict{}(
    "A653" => Dict{String,Any}("2state" => true , "plat" => [0.40,0.60], "pl_f" => 1.00, "mdof" => 4), # beta = 3.34
    "A654" => Dict{String,Any}("2state" => true , "plat" => [0.40,0.60], "pl_f" => 1.00, "mdof" => 4), # beta = 3.34
    "B450" => Dict{String,Any}("2state" => false, "plat" => [0.70,0.90], "pl_f" => 1.00, "mdof" => 4),
    "D450" => Dict{String,Any}("2state" => true , "plat" => [0.25,0.40], "pl_f" => 1.00, "mdof" => 4),
    "D452" => Dict{String,Any}("2state" => true , "plat" => [0.25,0.40], "pl_f" => 1.00, "mdof" => 4),
    "E250" => Dict{String,Any}("2state" => true , "plat" => [0.25,0.40], "pl_f" => 1.00, "mdof" => 50),

    "C101" => Dict{String,Any}("2state" => false, "plat" => [0.35,0.65], "pl_f" => 0.70, "mdof" => 4),
    "D200" => Dict{String,Any}("2state" => true , "plat" => [0.10,0.20], "pl_f" => 0.70, "mdof" => 4),
    "E300" => Dict{String,Any}("2state" => true , "plat" => [0.05,0.10], "pl_f" => 0.55, "mdof" => 4),
    "H101" => Dict{String,Any}("2state" => true , "plat" => [0.20,0.25], "pl_f" => 0.60, "mdof" => 4),
    "H102" => Dict{String,Any}("2state" => true , "plat" => [0.10,0.15], "pl_f" => 0.60, "mdof" => 4),
    "N101" => Dict{String,Any}("2state" => true , "plat" => [0.15,0.17], "pl_f" => 0.50, "mdof" => 4),
    "J303" => Dict{String,Any}("2state" => true , "plat" => [0.07,0.12], "pl_f" => 0.52, "mdof" => 4),
    "J500" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.48], "pl_f" => 0.48, "mdof" => 4), # beta = 3.85
    "N200" => Dict{String,Any}("2state" => true , "plat" => [0.05,0.15], "pl_f" => 0.65, "mdof" => 4),
    "N202" => Dict{String,Any}("2state" => true , "plat" => [0.05,0.10], "pl_f" => 0.65, "mdof" => 4),
    "N203" => Dict{String,Any}("2state" => true , "plat" => [0.08,0.13], "pl_f" => 0.65, "mdof" => 4),
    "N300" => Dict{String,Any}("2state" => true , "plat" => [0.05,0.15], "pl_f" => 0.60, "mdof" => 4),
    "N302" => Dict{String,Any}("2state" => false, "plat" => [0.28,0.35], "pl_f" => 0.70, "mdof" => 4),
    "S400" => Dict{String,Any}("2state" => true , "plat" => [0.13,0.18], "pl_f" => 0.60, "mdof" => 4),
)

NOERR_MDs  = true  # always in true!

AIC        = true  # always

PVAL       = false

PLOT_mDs   = false
PLOT_Ktar  = false

WRITE_mDS  = false
WRITE_Ktar = false
OVERWRITE  = false

path_bdio = path_bdio_dict["local"]

# Find MDs prime values in lattice units as a function of beta

println("- Computing mDs prime in the physical points and in the sym points...")

MD_ph_prime = MD_ph * (sqrtt0_ph_Bruno / sqrtt0_ph)  # shift in the mass, usefull to compute the derivative

uwerr(MD_ph_prime); println("   ⟹ MDs' = $(print_uwreal(MD_ph_prime)) GeV   [shifted mass] \n")

MDs_SU3_prime = Dict{String, uwreal}()
for beta in b_values
    MDs_SU3_prime["$beta"] = MD_ph_prime * (sqrtt0_ph / sqrt(t0sym(beta,Bruno=true))) / hbarc
end
MDs_SU3_prime
[uwerr(MDs_SU3_prime[key]) for key in keys(MDs_SU3_prime)]

if NOERR_MDs
    [set_fluc_to_zero!(MDs_SU3_prime[key], "sqrtt0 [fm] (Bruno)") for key in keys(MDs_SU3_prime)]
    [set_fluc_to_zero!(MDs_SU3_prime[key], "t0sym/a2") for key in keys(MDs_SU3_prime)]
    [MDs_SU3_prime[key] *= 1.0 for key in keys(MDs_SU3_prime)] 
end

if WRITE_mDS
    p = create_path(path_bdio,["mass&dec","mDs_prime"],OVERWRITE=OVERWRITE)

    io = IOBuffer()
    write(io, "MDs 'prime' physical and sym points")
    fb = ALPHAdobs_create(p, io)

    ALPHAdobs_write(fb, MD_ph_prime)
    ALPHAdobs_write(fb, MDs_SU3_prime)

    ALPHAdobs_close(fb)
end

@. fit_model(x,p) = 0.0  # the function needs to be initialized (maybe not?)

# println("- Reading s-heavy data and interpolating for kappaC target...")

kappaC_sim = get_kappa_values()

mDs = Dict{String,Array{uwreal}}()
kappaC_tar = Dict{String,uwreal}()

for ens in [EnsInfo(ensid)]  # ensInfo  [EnsInfo(ensid)]
    println("- Ens: $(ens.id)")

    data_sheavy = read_kappa_charm_all_config(joinpath(path_heavy,ens.id))

    println("   - Computing mDs...")

    bc = ens.bc
    T  = HVPobs.Data.get_T(ens.id)

    pl    = mDs_fitinfo[ens.id]["plat"]
    plf   = mDs_fitinfo[ens.id]["pl_f"]
    mdof  = mDs_fitinfo[ens.id]["mdof"]
    STATE = mDs_fitinfo[ens.id]["2state"]

    mDs[ens.id] = []
    for sh in ["sh1","sh2","sh3","sh4"]
        println("      - strange-heavy $(sh[end])")

        corr = corr_obs(data_sheavy[sh])

        obs = corr.obs

        println("         - Fitting...")
        
        if bc == "pbc"
            obs = obs[1:Int(T/2+1)]
            if STATE
                @. fit_model(x,p) = p[1] + p[2] * exp(-p[3]*x) + log((1 + exp(-p[1] * (T - 2*x))) / (1 + exp(-p[1] * (T - 2*(x+1)))))
                np = 3
            else
                @. fit_model(x,p) = p[1] + 0*x + log((1 + exp(-p[1] * (T - 2*x))) / (1 + exp(-p[1] * (T - 2*(x+1)))))
                np = 1
            end
        elseif bc == "obc"
            if STATE
                @. fit_model(x,p) = p[1] + p[2] * exp(-p[3]*x)
                np = 3
            else
                @. fit_model(x,p) = p[1] + 0*x
                np = 1
            end
        end

        m_obs = meff(obs)
        len   = length(m_obs)

        pl_f   = min(ceil(Int64,plf*len),len)
        p0_vec = collect(floor(Int64,pl[1]*len):1:ceil(Int64,min(pl_f-mdof,pl[2]*len)))

        if isempty(p0_vec)
            error("Not enough data to fit with d.o.f. ≥ $mdof. Decrease $mdof or increase the plateau search range. ")
        end

        wpm = ens.id in keys(wpmm) ? wpmm : nothing

        fit_vec = Vector{FitRes}()
        for p0 in p0_vec
            m_data = m_obs[p0:pl_f] 

            fit = fit_routine(fit_model,collect(p0:pl_f), m_data, np, pval=PVAL, wpm=wpm, info=false, lineprint=false)
            push!(fit_vec,fit)
        end

        println("         - Computing average...")

        w = get_w_from_fitres(vcat(fit_vec...), AIC=AIC)

        mDs_res_vec = [par[1] for par in getfield.(vcat(fit_vec...),:param)]

        mDs_res, mDs_syst = model_average(mDs_res_vec, w)
        mDs_res = mDs_res[1] + uwreal([0.0,mDs_syst],"mPP MA syst")

        push!(mDs[ens.id], mDs_res + uwreal([0.0,mDs_syst],"mDs syst"))

        if PLOT_mDs
            println("         - Plotting (mass)...")

            bestW = 20

            res     = mDs_res    ; uwerr.(res)
            res_vec = mDs_res_vec; uwerr.(res_vec)

            fig = figure(figsize=(16,12))

            gs = fig.add_gridspec(4, 1, height_ratios=[4, 1, 1, 1])  # Adjust the height_ratios as needed

            ax1 = fig.add_subplot(gs[1, 1])
            title("$(ens.id) (Ds)")

            x0 = collect(max(floor(Int64,len*(pl[1]))-2,1):len*plf+3) # .- 1
            vec = m_obs[Int.(x0)]; uwerr.(vec)

            errorbar(x0, value.(vec), err.(vec), fmt="o", capsize=2, color="black")
            maxw_arg = argmax(w[1:length(fit_vec)])
            par = getfield.(fit_vec,:param)[maxw_arg]

            if plf == 1.0
                x_fit = collect(max(floor(Int64,len*(pl[1])),1)-0.4:0.1:pl_f)
            else
                x_fit = collect(max(floor(Int64,len*(pl[1])),1)-0.4:0.1:pl_f+0.6)
            end
            y_fit = fit_model(x_fit,par); uwerr.(y_fit)

            fill_between(x_fit, value.(y_fit)+err.(y_fit), value.(y_fit)-err.(y_fit), alpha=0.3, color="orange")

            axvline(x=p0_vec[1]-0.4, color="red", linestyle="--")
            axvline(x=p0_vec[end]+0.4, color="red", linestyle="--")
            if plf != 1.0
                axvline(x=pl_f+0.6, color="orange", linestyle="dotted")
            end

            axis("tight")
            yscale("log")
            ylim(res.mean-8*res.err,maximum([value(vec[1]),res.mean+10*res.err]))
            # ylim(minimum(value.(vec)),maximum(value.(vec)))
            ylabel(L"$m^{\rm{eff}}_{\rm{s-h}}(t)$")
            
            setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels

            ax2 = fig.add_subplot(gs[2, 1])

            errorbar(p0_vec .- 1., value.(res_vec), err.(res_vec), fmt="d", mfc="none", color="blue")
            fill_between(x0, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")

            ylabel(L"$m_{\rm{Ds}}$")

            setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels


            ax3 = fig.add_subplot(gs[3, 1])

            fill_between(x0, maximum(w)/2, maximum(w)/2, alpha=0.0, color="white")
            PyPlot.plot(p0_vec .- 1., w, linestyle="none", marker="o", mfc="none", color="blue")

            ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))

            if PVAL
                setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels

                ax4 = fig.add_subplot(gs[4, 1])

                pval_vec = getfield.(fit_vec,:pval)
                fill_between(x0, maximum(vcat(pval_vec...))/2, maximum(vcat(pval_vec...))/2, alpha=0.0, color="white")
                PyPlot.plot(p0_vec .- 1., pval_vec, linestyle="none", marker="o", mfc="none", color="blue")

                ylabel(L"$\rm{p-values}$")
            end
            xlabel(L"$t/a$")

            tight_layout()
            display(fig)
            close()
        end
    end

    println("   - Computing KappaC target...")

    kappaCinv_sim = 1 ./ kappaC_sim[ens.id][2:end]

    @. lin_model(x,p) = p[1] + p[2] * x

    # fit = fit_routine(lin_model, kappaCinv_sim, mDs[ens.id], 2, info=false, lineprint=false)

    fit = lin_fit(kappaCinv_sim, mDs[ens.id], lineprint=false)
    par = fit[1]

    kappaCinv_tar = (MDs_SU3_prime["$(ens.beta)"] - par[1]) / par[2]
    kappaC_tar[ens.id] = 1/kappaCinv_tar

    uwerr(kappaC_tar[ens.id]); println("         ⟹ kappaC (tar) = $(print_uwreal(kappaC_tar[ens.id]))")

    if PLOT_Ktar
        println("   - Plotting...")

        kappaC_arr = range(kappaCinv_sim[1],kappaCinv_sim[end],100)
        mDs_arr = lin_model(kappaC_arr,par); uwerr.(mDs_arr)
        uwerr.(mDs[ens.id]); uwerr(MDs_SU3_prime["$(ens.beta)"]); uwerr(kappaCinv_tar); uwerr(kappaC_tar[ens.id])
        
        title("Ens: $(ens.id); "*L"$\kappa_c$"*" = $(print_uwreal(kappaC_tar[ens.id]))")
        errorbar(kappaCinv_sim, value.(mDs[ens.id]), err.(mDs[ens.id]), fmt="o", mfc="none", capsize=2, color="blue")
        errorbar(value(kappaCinv_tar), value(MDs_SU3_prime["$(ens.beta)"]), xerr = err(kappaCinv_tar), yerr = 0.0, fmt="^", mfc="none", capsize=2, color="red")
        fill_between(kappaC_arr, value.(mDs_arr)+err.(mDs_arr), value.(mDs_arr)-err.(mDs_arr), alpha=0.4, color="blue")
        axis("tight")
        ax = gca()      # get the handle of the current axis (not really used here)
        xlabel(L"$1/\kappa_c$")
        ylabel(L"$am_{D_s}$")
        display(gcf())      #display the figure
        close()
    end

    if WRITE_Ktar
        println("   - Printing BDIO...")

        p_mDs = create_path(path_bdio,["mass&dec",ens.id,"$(ens.id)_mDsKappa"],OVERWRITE=OVERWRITE)

        io = IOBuffer()
        write(io, "$(ens.id) mDs 'prime' & KappaC target")
        fb = ALPHAdobs_create(p_mDs, io)

        extra = Dict{String, Any}("ens" => ens.id)
        ALPHAdobs_write(fb, mDs[ens.id], extra=extra)
        ALPHAdobs_write(fb, kappaC_tar[ens.id], extra=extra)

        ALPHAdobs_close(fb)
    end
end # end ens loop

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> Reading test <==========================##

ensid = "A653"

path_bdio_r = path_bdio_dict["local"]

mDs_ph_prime, Ds_dict = BDIOread_mDs_kappaC(path_bdio_r,ensid)
