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

# include("const.jl")

# Plot specifications

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

##=============> Path definition, ensamble choice <=============##

julia_script_directory = @__DIR__

path_HVP  = joinpath(julia_script_directory, "..", "LatticeData", "HVP_data")
path_rw   = joinpath(julia_script_directory, "..", "LatticeData", "rwf_deflated")
path_ms   = joinpath(julia_script_directory, "..", "LatticeData", "ms_t0_dat")
path_fvc  = joinpath(julia_script_directory, "..", "LatticeData", "JKMPI")

path_coef = joinpath(julia_script_directory, "..", "Coefficients")

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

#ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302"]
#myensList, badensList = ensCheck(ensList, path_HVP, path_rw, path_ms, path_fvc, showbad=true)
subensList1 = ["H101", "B450", "N202", "N300"]
subensList2 = ["H102", "N101", "C101", "S400"]
subensList3 = ["N203", "N200", "D200", "N302"]
myensList, badensList = ensCheck(subensList3, path_HVP, path_rw, path_ms, path_fvc, showbad=true)

ensInfo = EnsInfo.(myensList)

isempty(badensList) ? @info("Enough information has been found for all ensembles") : @info("Not enough information has been found concerning ensables $(join(badensList, ", "))\n")

## <<<<=============================================>>>> #
## <<<<==========>>>> DATA EXTRACTION <<<<==========>>>> #
## <<<<=============================================>>>> #

IMPR      = true
IMPR_SET  = "2" # either "1" or "2"
RENORM    = true
STD_DERIV = false

@info("Impr. set: $IMPR_SET")

##=============> Reading, renormalization and impr. <=============##

corr = Dict{String, Dict{String, Vector{Corr}}}()
fvcHP = Dict{String, Vector{Corr}}()

t0 = Dict{String, uwreal}()

TMRa = Dict{String, Vector{uwreal}}()
TMRb = Dict{String, Vector{uwreal}}()
TMRc = Dict{String, Matrix{uwreal}}()

@info("==> Reading light data; renorm. and impr. <==")
for (i,ens) in enumerate(ensInfo)

    @info("Starting ensemble: $(ens.id)     ($i/$(length(ensInfo)))")
    @info("Reading corr. ...")
    Gll, Glc = corr_light(path_HVP, ens, path_rw = path_rw, frw_bcwd = true, renorm = RENORM, impr = IMPR, impr_set = IMPR_SET, cons = true, std = STD_DERIV, split = true)    
    corr[ens.id] = Dict("ll" => Gll, "lc" => Glc)

    @info("Reading FV correc. ...")
    fvcHPens = ∆GHP(path_fvc, ens, nmin=1, nmax=-1)
    fvcHP[ens.id] = fvcHPens

    @info("Reading t0 ...")
    t0ens = get_t0(path_ms, ens, path_rw = path_rw, pl = false)
    t0[ens.id] = t0ens

    @info("Computing TMR Kernels ...")

    sym_points = Int64(length(corr[ens.id]["ll"][1].obs)/2+1)
    
    factor = hbarc * sqrt(t0[ens.id])/t0_ph  
    TMRa[ens.id] = factor^2 .* Tildef4a((massmu/factor) .* collect(0:sym_points-1),path_coef)
    TMRb[ens.id] = factor^2 .* Tildef4b((massmu/factor) .* collect(0:sym_points-1),path_coef)
    TMRc[ens.id] = factor^4 .* Tildef4c((massmu/factor) .* collect(0:sym_points-1),path_coef)

    @info("-----------------------")
end

##=============> amu[NLO] HVP and FV corr. <=============##

HVP = Dict{String, Dict{String, Dict{String, Vector{uwreal}}}}()
HVP_int = Dict{String, Dict{String, Dict{String, Vector{Integrand}}}}()

HVP∆G = Dict{String, Dict{String, Vector{uwreal}}}()
HVP∆G_int = Dict{String, Dict{String, Vector{Integrand}}}()

@info("==> amu[NLO] HVP and FV corr. <==")
for (i,ens) in enumerate(ensInfo)

    @info("Starting ensemble: $(ens.id)     ($i/$(length(ensInfo)))")
    @info("Computing HVP ...")

    HVPall, Integrandall = amuHVPNLO("a", corr[ens.id]["ll"], TMRa[ens.id], t0ens=t0[ens.id], pl=false, int=true)
    HVPalc, Integrandalc = amuHVPNLO("a", corr[ens.id]["lc"], TMRa[ens.id], t0ens=t0[ens.id], pl=false, int=true)

    HVPbll, Integrandbll = amuHVPNLO("b", corr[ens.id]["ll"], TMRb[ens.id], t0ens=t0[ens.id], pl=false, int=true)
    HVPblc, Integrandblc = amuHVPNLO("b", corr[ens.id]["lc"], TMRb[ens.id], t0ens=t0[ens.id], pl=false, int=true)

    HVPcll, Integrandcll = amuHVPNLO("c", corr[ens.id]["ll"], TMRc[ens.id], t0ens=t0[ens.id], pl=false, int=true)
    HVPclc, Integrandclc = amuHVPNLO("c", corr[ens.id]["lc"], TMRc[ens.id], t0ens=t0[ens.id], pl=false, int=true)

    HVP[ens.id] = Dict("a" => Dict("ll" => HVPall, "lc" => HVPalc), "b" => Dict("ll" => HVPbll, "lc" => HVPblc), "c" => Dict("ll" => HVPcll, "lc" => HVPclc))
    HVP_int[ens.id] = Dict("a" => Dict("ll" => Integrandall, "lc" => Integrandalc), "b" => Dict("ll" => Integrandbll, "lc" => Integrandblc), "c" => Dict("ll" => Integrandcll, "lc" => Integrandclc))

    @info("Computing FV corr. ...")

    HVPa∆G, Integranda∆G = amu∆G("a", fvcHP[ens.id], TMRa[ens.id], pl=false, int=true)

    HVPb∆G, Integrandb∆G = amu∆G("b", fvcHP[ens.id], TMRb[ens.id], pl=false, int=true)

    HVPc∆G, Integrandc∆G = amu∆G("c", fvcHP[ens.id], TMRc[ens.id], corr=corr[ens.id]["ll"][end], pl=false, int=true)

    HVP∆G[ens.id] = Dict("a" => HVPa∆G, "b" => HVPb∆G, "c" => HVPc∆G)
    HVP∆G_int[ens.id] = Dict("a" => Integranda∆G, "b" => Integrandb∆G, "c" => Integrandc∆G)

    @info("-----------------------")
end

## <<<<==================================>>>> #
## <<<<==========>>>> BDIO <<<<==========>>>> #
## <<<<==================================>>>> #

##=============> Saving all observables in a BDIO file <=============##

path_bdio = joinpath(julia_script_directory,"..", "obsBDIO")
for (k, ensid) in enumerate(getfield.(ensInfo,:id))
    for diagram in ["a","b","c"]
        pens = joinpath(path_bdio, ensid)
        !ispath(pens) ? mkdir(pens) : nothing
        p = joinpath(pens, string(ensid,diagram,"_amuObs_Set$(IMPR_SET).bdio"))

        fb = BDIO_open(p, "w", ensid)
        uinfo = 0 

        write_uwreal(t0[ensid], fb, uinfo)
        uinfo +=1

        for i in collect(1:3)
            write_uwreal(HVP[ensid][diagram]["ll"][i], fb, uinfo)
            uinfo +=1
            write_uwreal(HVP[ensid][diagram]["lc"][i], fb, uinfo)
            uinfo +=1
        end
        
        [write_uwreal(v, fb, uinfo) for v in HVP∆G[ensid][diagram]]
        uinfo +=1

        BDIO_close!(fb)
    end
end

##for (k, ensid) in enumerate(getfield.(ensInfo,:id))
#    for diagram in ["a","b"]
#        pens = joinpath(path_bdio, ensid)
#        !ispath(pens) ? mkdir(pens) : nothing
#        p = joinpath(pens, string(ensid,diagram,"_amuObs.bdio"))
#
#        fb = BDIO_open(p, "w", ensid)
#        uinfo = 0 
#
#        for i in collect(1:3)
#            [write_uwreal(v, fb, uinfo) for v in HVP_int[ensid][diagram]["ll"][i].obs]
#            uinfo +=1
#            [write_uwreal(v, fb, uinfo) for v in HVP_int[ensid][diagram]["lc"][i].obs]
#            uinfo +=1
#        end
#
#        for i in collect(1:6)
#            [write_uwreal(v, fb, uinfo) for v in HVP∆G_int[ensid][diagram][i].obs]
#            uinfo +=1
#        end
#
#        BDIO_close!(fb)
#    end
#end

## <<<<===================================>>>> #
## <<<<==========>>>> PLOTS <<<<==========>>>> #
## <<<<===================================>>>> #

##=============> Plotting of integrand + FV corr. <=============##

gamma = "ll"

@info("==> Plotting of integrand + FV corr. for gamma=$gamma <==")
for (i,ens) in enumerate(ensInfo)
    @info("Starting ensemble: $(ens.id)     ($i/$(length(ensInfo)))")
    for diagram in ["a","b"]

        integrand = HVP_int[ens.id][diagram][gamma]
        ∆integrand = HVP∆G_int[ens.id][diagram]
        t_a = collect(0:length(integrand[1].obs)-1)

        label=[]
        for subintegrand in integrand
            uwerr.(subintegrand.obs)
            if subintegrand.gamma == "G08"*gamma
                if !all(x -> x == 0.0, value.(subintegrand.obs))
                    fill_between(t_a, -value.(subintegrand.obs)-ADerrors.err.(subintegrand.obs), -value.(subintegrand.obs)+ADerrors.err.(subintegrand.obs), alpha=1, label="-"*label_func(subintegrand.gamma))
                    push!(label,"-"*label_func(subintegrand.gamma))
                else
                    nothing
                end
            else
                fill_between(t_a, value.(subintegrand.obs)-ADerrors.err.(subintegrand.obs), value.(subintegrand.obs)+ADerrors.err.(subintegrand.obs), alpha=0.8, label=label_func(subintegrand.gamma))
                push!(label,label_func(subintegrand.gamma))
            end
        end
        PyPlot.title(ens.id*": Integrands for diagram (" * diagram * ")")
        xlabel("t/a")
        diagram == "a" ? (ylabel(L"$G(t)$ $\tilde{f}^{(a)}_4(m_\mu t)$")) : (ylabel(L"$G(t)$ $\tilde{f}^{(b)}_4(m_\mu t)$"))
        diagram == "a" ? (legend(label, loc="lower right")) : (legend(label, loc="upper right"))
        axis("tight")
        display(gcf())      #display the figure
        close("all")

        label=[]
        diagram == "a" ? (background=-ADerrors.err.(integrand[end].obs)) : (background=ADerrors.err.(integrand[end].obs))
        fill_between(t_a, 0, background, alpha=0.2, color="gray", label=L"$\delta G^{(R)}(t)$")
        push!(label,L"$\delta G^{(R)}(t)$")

        for ∆subintegrand in ∆integrand
            uwerr.(∆subintegrand.obs)
            fill_between(t_a, value.(∆subintegrand.obs)-ADerrors.err.(∆subintegrand.obs), value.(∆subintegrand.obs)+ADerrors.err.(∆subintegrand.obs), alpha=0.4, label=label_func(∆subintegrand.gamma))
            push!(label,label_func(∆subintegrand.gamma))
        end
        axis("tight")
        PyPlot.title(ens.id*": FV corrections for diagram (" * diagram * ")")
        xlabel("t/a")
        diagram == "a" ? (ylabel(L"$\Delta G(t)$ $\tilde{f}^{(a)}_4(m_\mu t)$")) : (ylabel(L"$\Delta G(t)$ $\tilde{f}^{(b)}_4(m_\mu t)$"))
        diagram == "a" ? (legend(label, loc="lower right")) : (legend(label, loc="upper right"))
        #grid("on")
        display(gcf())      #display the figure
        close("all")
    end
    @info("-----------------------")
end

##=============> Bounding Method plotting <=============##

tcut = 20
@info("==> Bounding Method plotting <==")
for (i,ens) in enumerate(ensInfo)

    @info("Ensemble: $(ens.id)     ($i/$(length(ensInfo)))")
    @info("Bounded correlator plots")
    mycorr33 = corr[ens.id]["ll"][1].obs
    mycorr88 = corr[ens.id]["ll"][2].obs
    mycorr08 = corr[ens.id]["ll"][3].obs
    uwerr.(mycorr33)
    uwerr.(mycorr88)
    uwerr.(mycorr08)

    sym_points = Int64(length(mycorr33)/2+1)
    t = collect(0:sym_points-1)
    
    E_eff33=Eeff(tcut, mycorr33)
    E_eff88=Eeff(tcut, mycorr88)
    E_eff08=Eeff(tcut, mycorr08)
    mpi = m_ens[ens.id]["m_pi"] 
    mrho = m_ens[ens.id]["m_rho"]
    L = ens.L
    E2pi = 2*sqrt(mpi^2 + (2*pi/L)^2)
    
    mrho < E2pi ? (UB33 = Gb(t, tcut, mycorr33, ens, mrho)) : (UB33 = Gb(t, tcut, mycorr33, ens, E2pi))
    LB33 = Gb(t, tcut, mycorr33, ens, E_eff33)
    UB88 = Gb(t, tcut, mycorr88, ens, mrho)
    LB88 = Gb(t, tcut, mycorr88, ens, E_eff88)
    uwerr.(UB33)
    uwerr.(LB33)
    uwerr.(UB88)
    uwerr.(LB88)

    legends = ["G33","G33 UB","G33 LB","G88","G88 UB","G88 LB"]
    errorbar(t, -value.(mycorr33[1:sym_points]), ADerrors.err.(mycorr33[1:sym_points]), fmt="s", label=legends[3], color="green", capsize=2) 
    errorbar(t[tcut+1:end], value.(UB33), ADerrors.err.(UB33), fmt="s", mfc="none", label=legends[1], color="limegreen", capsize=2)
    errorbar(t[tcut+1:end], value.(LB33), ADerrors.err.(LB33), fmt="s", mfc="none", label=legends[2], color="darkolivegreen", capsize=2)
  
    errorbar(t, -value.(mycorr88[1:sym_points]), ADerrors.err.(mycorr88[1:sym_points]), fmt="o", label=legends[6], color="blue", capsize=2)
    errorbar(t[tcut+1:end], value.(UB88), ADerrors.err.(UB88), fmt="o", mfc="none", label=legends[4], color="skyblue", capsize=2)
    errorbar(t[tcut+1:end], value.(LB88), ADerrors.err.(LB88), fmt="o", mfc="none", label=legends[5], color="dodgerblue", capsize=2)
    if value(mycorr08[1]) != 0.0
        UB08 = Gb(t, tcut, mycorr08, ens, mrho)
        LB08 = Gb(t, tcut, mycorr08, ens, E_eff08)
        uwerr.(UB08)
        uwerr.(LB08)

        legends08 = ["-G08","-G08 UB","-G08 LB"]
        errorbar(t, value.(mycorr08[1:sym_points]), ADerrors.err.(mycorr08[1:sym_points]), fmt="d", label=legends08[3], color="red", capsize=2)
        errorbar(t[tcut+1:end], -value.(UB08), ADerrors.err.(UB08), fmt="d", mfc="none", label=legends08[1], color="indianred", capsize=2)
        errorbar(t[tcut+1:end], -value.(LB08), ADerrors.err.(LB08), fmt="d", mfc="none", label=legends08[2], color="firebrick", capsize=2)
        legends = vcat(legends,legends08)
    end
    axis("tight")
    ax = gca()      # get the handle of the current axis (not really used here)
    yscale("log")
    xlabel("t/a")
    value(mycorr08[1]) != 0.0 ? (ylim((0.5*minimum(filter(x -> x >= 0, value.(mycorr08))), 2*maximum(-value.(mycorr33))))) : (ylim((0.5*minimum(filter(x -> x >= 0, -value.(mycorr33))), 2*maximum(-value.(mycorr33)))))
    ylabel("G(t)")
    PyPlot.title("Boundings for ensemble $(ens.id) and tcut = $tcut")
    legend(legends, loc  = "best")
    display(gcf())      #display the figure
    close()

    @info("Bounding Method applied")

    amu_BM("a", corr[ens.id]["ll"], TMRa[ens.id], bound_impr=true, t_step=1, t0ens=t0[ens.id], pl=true)
    amu_BM("b", corr[ens.id]["ll"], TMRb[ens.id], bound_impr=true, t_step=1, t0ens=t0[ens.id], pl=true)
    amu_BM("c", corr[ens.id]["ll"], TMRc[ens.id], bound_impr=true, t_step=1, t0ens=t0[ens.id], pl=true)

    @info("-----------------------")
end

## <<<<==================================>>>> #
## <<<<==========>>>> FITS <<<<==========>>>> #
## <<<<==================================>>>> # 

##=============> Continuum limit; Global chiral-continuum fit <=============##

phi2_ph = 8*(t0_ph*mpi_ph/hc)^2; uwerr(phi2_ph)
phi4_ph = 8*t0_ph^2*((mK_ph/hc)^2 + 0.5*(mpi_ph/hc)^2); uwerr(phi4_ph)     # for TrM=cst ensembles

#grouped_ensInfo = OrderedDict{Float64, Vector{EnsInfo}}()
#for info in ensInfo
#    if haskey(grouped_ensInfo, info.beta)
#        push!(grouped_ensInfo[info.beta], info)
#    else
#        grouped_ensInfo[info.beta] = [info]
#    end
#end

# The considered model is the simpler one can propose, a more concise fitting routine will be carried out in a separate Julia file

isov_basemodel(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+ p[4] .* x[:,1].^3/2. .+ p[5] .* (x[:,2].^2 .- value.(phi2_ph).^2)
# .+ p[4] .* (x[:,3] .- value.(phi4_ph))

ydata_ll   = Dict{String, Vector{uwreal}}()
ydata_lc   = Dict{String, Vector{uwreal}}()
fit_ll     = Dict{String, FitRes}()
fit_lc     = Dict{String, FitRes}()
res_val_ll = Dict{String, uwreal}()
res_val_lc = Dict{String, uwreal}()

xdata = hcat([[1 / (8*t0[ens.id]), 8*t0[ens.id]*m_ens[ens.id]["m_pi"]^2] for ens in ensInfo]...)
xdata = [xdata[1,:] xdata[2,:]]

for diagram in ["a","b"]
    ydata_ll[diagram] = [HVP[ens.id][diagram]["ll"][end]+HVP∆G[ens.id][diagram][end] for ens in ensInfo]; uwerr.(ydata_ll[diagram])
    ydata_lc[diagram] = [HVP[ens.id][diagram]["lc"][end]+HVP∆G[ens.id][diagram][end] for ens in ensInfo]; uwerr.(ydata_lc[diagram])

    fit_ll[diagram] = fit_routine(isov_basemodel,value.(xdata),ydata_ll[diagram],5)
    fit_lc[diagram] = fit_routine(isov_basemodel,value.(xdata),ydata_lc[diagram],5)
    res_val_ll[diagram] = fit_ll[diagram].param[1]; uwerr(res_val_ll[diagram])
    res_val_lc[diagram] = fit_lc[diagram].param[1]; uwerr(res_val_lc[diagram])
end

## == Projection to continuum Plot ==

for diagram in ["a","b"]
    xdata_proj = [xdata[:,1] fill(phi2_ph, length(xdata[:,1]))]; uwerr.(xdata_proj)
    ydata_proj_ll = isov_basemodel(xdata_proj, fit_ll[diagram].param); uwerr.(ydata_proj_ll)
    ydata_proj_lc = isov_basemodel(xdata_proj, fit_lc[diagram].param); uwerr.(ydata_proj_lc)

    xarr = [Float64.(range(0.0, maximum(value.(xdata[:,1])), length=100)) fill(phi2_ph, 100)]
    yarr_ll = isov_basemodel(xarr, fit_ll[diagram].param); uwerr.(yarr_ll)
    yarr_lc = isov_basemodel(xarr, fit_lc[diagram].param); uwerr.(yarr_lc)

    errorbar(value.(xdata_proj[:,1]), xerr=ADerrors.err.(xdata_proj[:,1]), value.(ydata_proj_ll), yerr=ADerrors.err.(ydata_proj_ll), fmt="s", capsize=2, color="blue", mfc="none", label="local-local")
    errorbar(value.(xdata_proj[:,1]), xerr=ADerrors.err.(xdata_proj[:,1]), value.(ydata_proj_lc), yerr=ADerrors.err.(ydata_proj_lc), fmt="o", capsize=2, color="indianred", mfc="none", label="local-conserved")
    errorbar(0.0, value(res_val_ll[diagram]), ADerrors.err(res_val_ll[diagram]), fmt="s", capsize=2, color="purple", mfc="none")
    errorbar(0.0, value(res_val_lc[diagram]), ADerrors.err(res_val_lc[diagram]), fmt="o", capsize=2, color="red", mfc="none")

    fill_between(xarr[:,1], value.(yarr_ll)-ADerrors.err.(yarr_ll), value.(yarr_ll)+ADerrors.err.(yarr_ll), alpha=0.2, color="blue")
    fill_between(xarr[:,1], value.(yarr_lc)-ADerrors.err.(yarr_lc), value.(yarr_lc)+ADerrors.err.(yarr_lc), alpha=0.2, color="indianred")

    PyPlot.title("Projection to continuum extrapolation")
    axvline(ls="dashed", color="black", lw=0.2, alpha=0.7) 
    xlabel(L"$a^2/8t_0$")
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{NLO}_{$diagram}]"))
    legend(["local-local","local-conserved"], loc="best")
    tight_layout()
    display(gcf())
    close()
end

## == Projection to chiral Plot ==

fmt = ["o","^","s","d","*"]
color = ["orange","red","purple","blue","green"]

for diagram in ["a","b"]
    label = []
    for (k,b) in  enumerate(sort(unique(getfield.(ensInfo, :beta))))
        push!(label,L"$\beta=$"*"$b")
        n_ = findall(x->x.beta == b, ensInfo)
        a2_aux = mean(value.(xdata[:,1][n_]))
        errorbar(value.(xdata[n_,2]), value.(ydata_ll[diagram][n_]), ADerrors.err.(ydata_ll[diagram][n_]), fmt=fmt[k], capsize=2, color=color[k], mfc="none", label=label[k])
        xxx = [fill(a2_aux,100) Float64.(range(0.0, 0.8, length=100))]
        yyy_ll = isov_basemodel(xxx, fit_ll[diagram].param)
        PyPlot.plot(xxx[:,2],value.(yyy_ll),ls="--",color=color[k],lw=0.5)
    end
    xxx_ph = [fill(0.,100) Float64.(range(0.0, 0.8, length=100))]
    yyy_ph_ll = isov_basemodel(xxx_ph, fit_ll[diagram].param); uwerr.(yyy_ph_ll)
    errorbar(value(phi2_ph), value(res_val_ll[diagram]), ADerrors.err(yyy_ph_ll[1]), fmt="^", capsize=2, color="black",label="ph.")
    push!(label,"ph.")
    PyPlot.plot(xxx_ph[:,2],value.(yyy_ph_ll),ls="--",color="gray",lw=0.5)
    fill_between(xxx_ph[:,2], value.(yyy_ph_ll)-ADerrors.err.(yyy_ph_ll), value.(yyy_ph_ll)+ADerrors.err.(yyy_ph_ll), alpha=0.2, color="gray")
    PyPlot.title("Projection to chiral extrapolation (local-local)")
    axvline(x=value(phi2_ph), ls="dashed", color="black", lw=0.2, alpha=0.7) 
    xlabel(L"$\Phi_2$")
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{NLO}_{$diagram}]"))
    legend(label,loc="best")
    tight_layout()
    display(gcf())
    close()

    label = []
    for (k,b) in  enumerate(sort(unique(getfield.(ensInfo, :beta))))
        push!(label,L"$\beta=$"*"$b")
        n_ = findall(x->x.beta == b, ensInfo)
        a2_aux = mean(value.(xdata[:,1][n_]))
        errorbar(value.(xdata[n_,2]), value.(ydata_lc[diagram][n_]), ADerrors.err.(ydata_lc[diagram][n_]), fmt=fmt[k], capsize=2, color=color[k], mfc="none", label=label[k])
        xxx = [fill(a2_aux,100) Float64.(range(0.0, 0.8, length=100))]
        yyy_lc = isov_basemodel(xxx, fit_lc[diagram].param)
        PyPlot.plot(xxx[:,2],value.(yyy_lc),ls="--",color=color[k],lw=0.5)
    end
    xxx_ph = [fill(0.,100) Float64.(range(0.0, 0.8, length=100))]
    yyy_ph_lc = isov_basemodel(xxx_ph, fit_lc[diagram].param); uwerr.(yyy_ph_lc)
    errorbar(value(phi2_ph), value(res_val_lc[diagram]), ADerrors.err(yyy_ph_lc[1]), fmt="^", capsize=2, color="black",label="ph.")
    push!(label,"ph.")
    PyPlot.plot(xxx_ph[:,2],value.(yyy_ph_lc),color="gray",lw=0.5)
    fill_between(xxx_ph[:,2], value.(yyy_ph_lc)-ADerrors.err.(yyy_ph_lc), value.(yyy_ph_lc)+ADerrors.err.(yyy_ph_lc), alpha=0.2, color="gray")
    PyPlot.title("Projection to chiral extrapolation (local-conserved)")
    axvline(x=value(phi2_ph), ls="dashed", color="black", lw=0.2, alpha=0.7) 
    xlabel(L"$\Phi_2$")
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{NLO}_{$diagram}]"))
    legend(label,loc="best")
    tight_layout()
    display(gcf())
    close()
end
