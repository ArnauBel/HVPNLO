# Import packages

using HVPobs

using ALPHAio, ADerrors
import ADerrors: err

using Statistics

using BDIO
using JLD2

using Plots
using PyPlot
using Colors

using ProgressBars
using Suppressor

# Include usefull functions

include("../HVPtools/Utils.jl")

include("../HVPtools/Fit&MA.jl")

# Include Isovector Model (needed constants already included inside)

include("../isovModel.jl")

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

# BDIO path definition

julia_script_directory = @__DIR__

STD_DERIV = true

if STD_DERIV
    @info("Standard derivative has been chosen !!")
    path_bdio = joinpath(julia_script_directory, "..", "ObsBDIOstd")
else
    path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")
end

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> Import data for Plots: corr, TMR and FVC <==========================##

# Provide ens, impr. set and discretisation

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ens      = EnsInfo("N101")
impr_set = "1"

@info("Ens $(ens.id) and impr. set $impr_set")

println("- Reading t0...")

t0 = uwreal(0.0) 
fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_t0"),"r")
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    t0 = ALPHAdobs_read_next(fb, size=sz, keys=ks)["t0"][1]
end

println("- Reading TMR...")

fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_TMR"),"r") 
TMR = Dict{String, Any}()
partial_res = Vector{Dict}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    push!(partial_res,ALPHAdobs_read_next(fb, size=sz, keys=ks))
end
BDIO_close!(fb)
for dict in partial_res
    merge!(TMR, dict)
end
TMRa = TMR["TMRa"]
TMRb = TMR["TMRb"]
TMRc = TMR["TMRc"]

println("- Reading corr...")

fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_corr_set"*impr_set),"r")
corr = Dict{String, Array{uwreal}}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    corr = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

println("- Reading FVC...")

FVC_corr = Dict{String,Vector{uwreal}}()
fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_FVC"),"r") 
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    FVC_corr = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

@info("All data for 'corr, TMR and FVC' read. ")

##==========================> Import data for Plots: HVP, HVP info and FVC <==========================##

ensid = "C101"
wind  = "NW"
diag  = "NLOa&b"
impr_set = "1"

ens = EnsInfo(ensid)

pens = joinpath(path_bdio,"HVP&FVC",wind,ensid)

fb = BDIO_open(joinpath(pens,"HVP","$(ensid)_HVP$(diag)_set$(impr_set)"),"r")
val = Dict{String, Dict{String, uwreal}}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    ks = collect(d["keys"])
    val["HVP"] =  ALPHAdobs_read_next(fb, keys=ks)
end
BDIO_close!(fb)
if wind ∉ ["SD","ID"]
    info = load(joinpath(pens,"HVP","$(ensid)_HVP$(diag)_BMinfo_set$(impr_set).jld2"), "HVPinfo")
    HVP = merge(val,info) 
else
    HVP = val
end

if diag == "NLOc"
    fb = BDIO_open(joinpath(pens,"FVC","$(ensid)_FVC$(diag)_set$(impr_set)"),"r")
else
    fb = BDIO_open(joinpath(pens,"FVC","$(ensid)_FVC$(diag)"),"r")    
end
FVC = Dict{String, Vector{uwreal}}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    FVC = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

##==========================> Results, Fits & MA <==========================##

DictComptoKey = Dict{String,Vector{String}}(
    "33"      => ["g33_ll","g33_lc"],
    "88"      => ["g88_ll","g88_lc"],
    "cc conn" => ["gcc_ll_conn","gcc_lc_conn"],
    "cc disc" => ["gcc_cc_disc"],
    "c8 disc" => ["gc8_cc_disc"],

    "3333"    => ["g3333_ll","g3333_lc"],
    "8888"    => ["g8888_ll","g8888_lc"],
    "CCCC"    => ["gCCCC_ll","gCCCC_lc"],
    "3388"    => ["g3388_ll","g3388_lc"],
    "33CC"    => ["g33CC_ll","g33CC_lc"],
    "88CC"    => ["g88CC_ll","g88CC_lc"]
)

wind = "ID"
diag = "LO"
comp = "33"

type_basemodel = "phi4"      #  phi4  simple
type_DA        = "All-(a4,a2loga)"  # All  All-(a4)  All-(a4,a2loga)  All-(a4,a2phi4)  ||  Only(phi2inv,logphi2)  Only(a3cutoff,phi2inv,logphi2)

IMPR_SET = ["1old","2"] # ["1","2"]

fittype_dict = joinpath(path_bdio,"Fit&MA",wind,"base[$type_basemodel]",type_DA)

@info("Comp. $(comp)")

mykeys = DictComptoKey[comp]

println("- Reading FitRes...")

fitres = load(joinpath(fittype_dict,"Fit",diag,comp,"FitRes.jld2"), "FitRes")

println("- Reading res & param...")

res   = Dict{String, Any}()
param = Dict{String, Any}()

modelinfo = load(joinpath(fittype_dict,"Fit",diag,comp,"ModelInfo.jld2"), "info")

fb = BDIO_open(joinpath(fittype_dict,"Fit",diag,comp,"param"),"r")
partial_res = Vector{Dict}()
full_dict = Dict{String, Any}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    push!(partial_res,ALPHAdobs_read_next(fb, size=sz, keys=ks))
end
BDIO_close!(fb)
for dict in partial_res
    merge!(full_dict, dict)
end
res   = Dict{String, Any}()
param = Dict{String, Any}()
for impr_set in IMPR_SET
    res[impr_set]   = Dict{String, Vector{uwreal}}()
    param[impr_set] = Dict{String, Vector{Vector{uwreal}}}()
    for key in mykeys
        res[impr_set][key]   = []
        param[impr_set][key] = []
        for i in collect(1:modelinfo["length"])
            push!(res[impr_set][key]  , full_dict["$(diag)_$(key)_set$(impr_set):[$i]"][1])
            push!(param[impr_set][key], full_dict["$(diag)_$(key)_set$(impr_set):[$i]"])
        end
    end
end

println("- Reading 'x' & 'y' data...")

fb = BDIO_open(joinpath(fittype_dict,"Fit",diag,comp,"xydata"),"r")
partial_res = Vector{Dict}()
full_dict = Dict{String, Any}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    push!(partial_res,ALPHAdobs_read_next(fb, size=sz, keys=ks))
end
BDIO_close!(fb)
for dict in partial_res
    merge!(full_dict, dict)
end
xydata = Dict{String, Any}()
xydata["xdata"] = full_dict["xdata"]
xydata["ydata"] = Dict{String, Dict}()
for impr_set in IMPR_SET
    xydata["ydata"][impr_set] = Dict{String, Vector{uwreal}}()
    for key in mykeys
        xydata["ydata"][impr_set][key] = full_dict["$(key)_set$(impr_set)"]
    end
end


println("- Reading MA...")

MA_jld2 = load(joinpath(fittype_dict,"MA",diag,comp,"MA_info.jld2"), "MAinfo")

fb = BDIO_open(joinpath(fittype_dict,"MA",diag,comp,"MA"),"r")
i=0; reskeys = ["res_tot","amu"]
MA = Dict{String,Any}()
while ALPHAdobs_next_p(fb)
    i += 1
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    MA[reskeys[i]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)[reskeys[i]][i == 1 ? (1:end) : 1]
end
BDIO_close!(fb)

MA["syst"] = MA_jld2["syst"]
MA["weight_tot"] = MA_jld2["weight_tot"]


println("- Reading info")

myinfo = load(joinpath(fittype_dict,"Fit",diag,comp,"ModelInfo.jld2"), "info")

println("- Reading list of ensembles...")

ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302","E250","J303","E300"] #J500
ensInfo = EnsInfo.(ensList)

println("- Importing isov. model...")

mdof = comp in ["cc disc","c8 disc"] ? 3 : 4   # impose a minimum to the dof for all fits 
nens = myinfo["nens"]
f_tot_isov, n_par_tot_isov, label_tot_isov = call_models(type_basemodel,type_DA,mdof,nens)

@info("All data for 'HVP, HVP info, FVC and results (MA)' read. ")

##======================================================================##
##==========================> (1) Corr Plots <==========================##
##======================================================================##

#-- Plot 1: Corr 

@info("Plot 1: Corr")

gamma_list = ["g33_ll","g88_ll_conn","gcc_ll_conn","g08_ll_conn","gc8_ll_disc"]

println("   - Ploting for $gamma_list")

points = collect(1:Int64(length(corr["g33_ll"])/2+1))

fmt_list    = ["o","s","d","^","*","x","·"]
color_list  = ["red","orange","green","blue","purple","brown","gray"]
legend_list = Vector{String}()
ymax = 0; ymin = 10
for (i,gamma) in enumerate(gamma_list)
    if gamma == "g08_cc_disc"
        obs = corr["g08_cc_disc"][points].+corr["g80_cc_disc"][points]; uwerr.(obs)
    elseif gamma == "gc8_cc_disc"
        obs = corr["gc8_cc_disc"][points].+corr["g8c_cc_disc"][points]; uwerr.(obs)
    else
        obs = corr[gamma][points]; uwerr.(obs)
    end
    errorbar(points, abs.(value.(obs)), err.(obs), fmt=fmt_list[i], mfc="none", label=gamma, color=color_list[i], capsize=2)
    push!(legend_list,gamma)
    ymax = maximum(vcat(abs.(value.(obs)),ymax))
    length(gamma_list)==1 || gamma!="gcc_ll_conn" ? ymin = minimum(vcat(abs.(value.(obs)),ymin)) : nothing
end
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
yscale("log")
ylim(ymin*1e-1,ymax*1e1)
xlabel("t/a")
ylabel(latexstring("|G(t)|"))
legend(legend_list, loc  = "best")
display(gcf())      #display the figure
close()

##=================================================================================##
##==========================> (2) Bounding Method Plots <==========================##
##=================================================================================##


discr = "ll"


if ens.kappa_l == ens.kappa_s
    obs33 = corr["g33_"*discr]; uwerr.(obs33)
    obs88 = corr["g88_"*discr*"_conn"]; uwerr.(obs88)
    obsCC = corr["gcc_"*discr*"_conn"]; uwerr.(obsCC)

    obs88R = obs88
    obs_ls = obs33 .+ (1/3).*obs88
    obs = obs33 .+ (1/3).*obs88 .+ (4/9).*obsCC; uwerr.(obs)
else
    obs33 = corr["g33_"*discr]; uwerr.(obs33)
    obs88 = corr["g88_"*discr*"_conn"] .+ corr["g88_"*discr*"_disc"]; uwerr.(obs88)
    if discr == "ll"
        obs08 = (2).*corr["g08_"*discr*"_conn"] .+ corr["g08_"*discr*"_disc"] .+ corr["g80_"*discr*"_disc"]; uwerr.(obs08)
    elseif discr == "lc"
        obs08 = corr["g08_"*discr*"_conn"] .+ corr["g80_"*discr*"_disc"]; uwerr.(obs08)
    end
    obsCCconn = corr["gcc_"*discr*"_conn"]; uwerr.(obsCCconn)
    obsCCdisc = corr["gcc_cc_disc"]; uwerr.(obsCCdisc)
    obsC8 = corr["gc8_cc_disc"] .+ corr["g8c_cc_disc"]; uwerr.(obsC8)

    obs88R = obs88.+obs08
    obs_ls = obs33 .+ (1/3).*(obs88.+obs08)
    obs = obs33 .+ (1/3).*(obs88.+obs08) .+ (4/9).*obsCCconn .+ (4/9).*obsCCdisc .+ (2/(3*sqrt(3))).*obsC8   # .+ (4/(3*sqrt(3))).*obsC8
end

t = collect(1:Int64(length(obs33)/2+1))
tcut = 10

Eeff33=Eeff(tcut, obs33)
Eeff88=Eeff(tcut, obs88)
ens.kappa_l != ens.kappa_s ? Eeff08=Eeff(tcut, obs08) : nothing
mpi = m_ens[ens.id]["m_pi"] 
mrho = m_ens[ens.id]["m_rho"]
L = ens.L
E2pi = 2*sqrt(mpi^2 + (2*pi/L)^2)
E3pi = 2*sqrt(mpi^2 + (2*pi/L)^2) + sqrt(mpi^2 + 2(2*pi/L)^2)

#-- Plot 2.1: Correlator bounding

@info("Plot 2.1: Correlator bounding")

UB33 = mrho < E2pi ? Gbounding(t, tcut, obs33, ens, mrho) : Gbounding(t, tcut, obs33, ens, E2pi); uwerr.(UB33)
LB33 = Gbounding(t, tcut, obs33, ens, Eeff33); uwerr.(LB33)
UB88 = mrho < E3pi ? Gbounding(t, tcut, obs88, ens, mrho) : Gbounding(t, tcut, obs88, ens, E3pi); uwerr.(UB88)
LB88 = Gbounding(t, tcut, obs88, ens, Eeff88); uwerr.(LB88)


legends = ["G33","G33 UB","G33 LB","G88","G88 UB","G88 LB"]
errorbar(t, value.(obs33[t]), err.(obs33[t]), fmt="s", label=legends[3], color="green", capsize=2) 
errorbar(t[tcut+1:end], value.(UB33), err.(UB33), fmt="s", mfc="none", label=legends[1], color="limegreen", capsize=2)
errorbar(t[tcut+1:end], value.(LB33), err.(LB33), fmt="s", mfc="none", label=legends[2], color="darkolivegreen", capsize=2)

errorbar(t, value.(obs88[t]), err.(obs88[t]), fmt="o", label=legends[6], color="blue", capsize=2)
errorbar(t[tcut+1:end], value.(UB88), err.(UB88), fmt="o", mfc="none", label=legends[4], color="skyblue", capsize=2)
errorbar(t[tcut+1:end], value.(LB88), err.(LB88), fmt="o", mfc="none", label=legends[5], color="dodgerblue", capsize=2)
if ens.kappa_l != ens.kappa_s
    UB08 = Gbounding(t, tcut, obs08, ens, mrho   ); uwerr.(UB08)
    LB08 = Gbounding(t, tcut, obs08, ens, Eeff08); uwerr.(LB08)

    legends08 = ["-G08","-G08 UB","-G08 LB"]
    errorbar(t, -value.(obs08[t]), err.(obs08[t]), fmt="d", label=legends08[3], color="red", capsize=2)
    errorbar(t[tcut+1:end], -value.(UB08), err.(UB08), fmt="d", mfc="none", label=legends08[1], color="indianred", capsize=2)
    errorbar(t[tcut+1:end], -value.(LB08), err.(LB08), fmt="d", mfc="none", label=legends08[2], color="firebrick", capsize=2)
    legends = vcat(legends,legends08)
end
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
yscale("log")
xlabel("t/a")
ymax = maximum(abs.(value.(obs33)))
ymin = minimum(abs.(value.(obs08)))
ylim(ymin*1e-1,ymax*1e1)
ylabel("G(t)")
PyPlot.title("Boundings for ensemble $(ens.id) and tcut = $tcut")
legend(legends, loc  = "best")
display(gcf())      #display the figure
close()

#-- Plot 2.2: Product of correlators bounding

@info("Plot 2.2: Product of correlators bounding")

obs3333 = obs33.*obs33; uwerr.(obs3333)
obs8888 = obs88R.*obs88R; uwerr.(obs8888)
obs3388 = obs33.*obs88; uwerr.(obs3388)

Eeff3333=Eeff(tcut, obs3333)
Eeff8888=Eeff(tcut, obs8888)
Eeff3388=Eeff(tcut, obs3388)


UB3333 = mrho < E2pi ? Gbounding(t, tcut, obs3333, ens, 2*mrho) : Gbounding(t, tcut, obs3333, ens, 2*E2pi); uwerr.(UB3333)
LB3333 = Gbounding(t, tcut, obs3333, ens, Eeff3333); uwerr.(LB3333)
UB8888 = mrho < E3pi ? Gbounding(t, tcut, obs8888, ens, 2*mrho) : Gbounding(t, tcut, obs8888, ens, 2*E3pi); uwerr.(UB8888)
LB8888 = Gbounding(t, tcut, obs8888, ens, Eeff8888); uwerr.(LB8888)
UB3388 = mrho < E2pi ? Gbounding(t, tcut, obs3388, ens, 2*mrho) : (mrho < E3pi ? Gbounding(t, tcut, obs3388, ens, E2pi+mrho) : Gbounding(t, tcut, obs3388, ens, E2pi+E3pi)); uwerr.(UB3388)
LB3388 = Gbounding(t, tcut, obs3388, ens, Eeff3388); uwerr.(LB3388)

legends3333 = ["(G33xG33)","(G33xG33) UB","(G33xG33) LB"]

errorbar(t, value.(obs3333[t]), err.(obs3333[t]), fmt="s", label=legends3333[1], color="green", capsize=2) 
errorbar(t[tcut+1:end], value.(UB3333), err.(UB3333), fmt="s", mfc="none", label=legends3333[2], color="limegreen", capsize=2)
errorbar(t[tcut+1:end], value.(LB3333), err.(LB3333), fmt="s", mfc="none", label=legends3333[3], color="darkolivegreen", capsize=2)

axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
yscale("log")
xlabel("t/a")
ymax = maximum(abs.(value.(obs3333)))
ymin = minimum(abs.(value.(obs3333)))
ylim(ymin*1e-1,ymax*1e1)
ylabel(L"G^{33}(t)\times G^{33}(t)")
PyPlot.title("Boundings for ensemble $(ens.id) and tcut = $tcut")
legend(legends, loc  = "best")
display(gcf())      #display the figure
close()

legends8888 = ["(G88xG88)","(G88xG88) UB","(G88xG88) LB"]

errorbar(t, value.(obs8888[t]), err.(obs8888[t]), fmt="o", label=legends8888[1], color="blue", capsize=2) 
errorbar(t[tcut+1:end], value.(UB8888), err.(UB8888), fmt="o", mfc="none", label=legends8888[2], color="skyblue", capsize=2)
errorbar(t[tcut+1:end], value.(LB8888), err.(LB8888), fmt="o", mfc="none", label=legends8888[3], color="dodgerblue", capsize=2)

axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
yscale("log")
xlabel("t/a")
ymax = maximum(abs.(value.(obs8888)))
ymin = minimum(abs.(value.(obs8888)))
ylim(ymin*1e-1,ymax*1e1)
ylabel(L"G^{88}(t)\times G^{88}(t)")
PyPlot.title("Boundings for ensemble $(ens.id) and tcut = $tcut")
legend(legends, loc  = "best")
display(gcf())      #display the figure
close()

legends3388 = ["(G33xG88)","(G33xG88) UB","(G33xG88) LB"]

errorbar(t, value.(obs3388[t]), err.(obs3388[t]), fmt="d", label=legends3388[1], color="red", capsize=2) 
errorbar(t[tcut+1:end], value.(UB3388), err.(UB3388), fmt="d", mfc="none", label=legends3388[2], color="indianred", capsize=2)
errorbar(t[tcut+1:end], value.(LB3388), err.(LB3388), fmt="d", mfc="none", label=legends3388[3], color="firebrick", capsize=2)

axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
yscale("log")
xlabel("t/a")
ymax = maximum(abs.(value.(obs3388)))
ymin = minimum(abs.(value.(obs3388)))
ylim(ymin*1e-1,ymax*1e1)
ylabel(L"G^{33}(t)\times G^{88}(t)")
PyPlot.title("Boundings for ensemble $(ens.id) and tcut = $tcut")
legend(legends, loc  = "best")
display(gcf())      #display the figure
close()

## The charm charm is so short distance that no bounding method is required:

obs33CC = obs33.*obsCCconn; uwerr.(obs33CC)

errorbar(t, value.(obs33CC[t]), err.(obs33CC[t]), fmt="o", mfc="none", label="(G33xGCC)", color="orange", capsize=2) 

axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
yscale("log")
xlabel("t/a")
ymax = maximum(abs.(value.(obs33CC)))
ymin = minimum(abs.(value.(obs33CC)))
ylim(ymin*1e-1,ymax*1e1)
ylabel(L"G^{33}\times G^{\rm{CC}}")
PyPlot.title("Boundings for ensemble $(ens.id) and tcut = $tcut")
display(gcf())      #display the figure
close()


#-- Plot 2.3: Bounding Method (33 & 88) [a & b]


@info("Plot 2.3: Bounding Method (33 & 88) [a & b]")

tcut0 = 10
tstep = 1

LBOUND_IMPR = true; tcut_fix = 1.2  # fm

DIAG = ["b"] # ["a","b"] # ["a","b","c"]

for diag in DIAG
println("   - Starting diagram $diag")

    aens = t0_ph / sqrt(t0)
    tcut_fm =  value.(aens .* (collect(tcut0:tstep:t[end-1]).-1))

    ub33 = Vector{uwreal}()
    ub88 = Vector{uwreal}()

    lb33 = Vector{uwreal}()
    lb88 = Vector{uwreal}()

    if LBOUND_IMPR
        lb_impr33 = Vector{uwreal}()
        lb_impr88 = Vector{uwreal}()
    end

    obs88R = ens.kappa_l == ens.kappa_s ? obs88 : obs88.+obs08

    if diag in ["a", "b"]
        TMR = diag == "a" ? TMRa : TMRb 
        int33 = obs33[t] .* TMR
        int88 = (obs88R[t]) .* TMR

        println("      - Computing bounded estimations...")
        for tcut in tcut0:tstep:t[end-1]

            UB33 = mrho < E2pi ? Gbounding(t, tcut, obs33, ens, mrho) : Gbounding(t, tcut, obs33, ens, E2pi)
            UB88 = Gbounding(t, tcut, obs88R, ens, mrho)

            UBInt33 = UB33 .* TMR[tcut+1:end]
            UBInt88 = UB88 .* TMR[tcut+1:end]

            lb_amuNLO33 = (alpha/pi)^3 * sum(int33[1:tcut]) * 1e10
            lb_amuNLO88 = (alpha/pi)^3 * sum(int88[1:tcut]) * 1e10

            ub_amuNLO33 = (alpha/pi)^3 * (sum(int33[1:tcut])+sum(UBInt33)) * 1e10
            ub_amuNLO88 = (alpha/pi)^3 * (sum(int88[1:tcut])+sum(UBInt88)) * 1e10

            push!(lb33, lb_amuNLO33)
            push!(lb88, lb_amuNLO88)

            push!(ub33, ub_amuNLO33)
            push!(ub88, ub_amuNLO88)

            if LBOUND_IMPR
                if aens*tcut < tcut_fix  # we fix the eff energy at some point
                    Eeff33=Eeff(tcut, obs33)
                    Eeff88=Eeff(tcut, obs88R)
                end

                LB33 = Gbounding(t, tcut, obs33, ens, Eeff33)
                LB88 = Gbounding(t, tcut, obs88R, ens, Eeff88)

                LBint33 = LB33 .* TMR[tcut+1:end]
                LBint88 = LB88 .* TMR[tcut+1:end]

                lb_impr_amuNLO33 = (alpha/pi)^3 * (sum(int33[1:tcut])+sum(LBint33)) * 1e10
                lb_impr_amuNLO88 = (alpha/pi)^3 * (sum(int88[1:tcut])+sum(LBint88)) * 1e10

                push!(lb_impr33, lb_impr_amuNLO33)
                push!(lb_impr88, lb_impr_amuNLO88)
            end
        end
    elseif diag == "c"
        int33 = (obs33[t] .* hcat(obs33[t]...)) .* TMRc
        int88 = ((obs88R[t]) .* hcat(obs88R[t]...)) .* TMRc

        println("      - Computing bounded estimations...")
        for tcut in tcut0:tstep:t[end-1]

            UB33 =  mrho < E2pi ? Gbounding(t, tcut, obs33, ens, mrho) : Gbounding(t, tcut, obs33, ens, E2pi)
            UB88 = Gbounding(t, tcut, obs88R, ens, mrho)

            UBInt33_1 = (UB33 .* hcat(UB33...)) .* TMRc[tcut+1:end,tcut+1:end]
            UBInt88_1 = (UB88 .* hcat(UB88...)) .* TMRc[tcut+1:end,tcut+1:end]
            UBInt33_2 = (obs33[1:tcut] .* hcat(UB33...)) .* TMRc[1:tcut,tcut+1:end]
            UBInt88_2 = ((obs88R[1:tcut]) .* hcat(UB88...)) .* TMRc[1:tcut,tcut+1:end]

            lb_amuNLO33 = (alpha/pi)^3 * sum(int33[1:tcut,1:tcut]) * 1e10
            lb_amuNLO88 = (alpha/pi)^3 * sum(int88[1:tcut,1:tcut]) * 1e10
            ub_amuNLO33 = (alpha/pi)^3 * (sum(int33[1:tcut,1:tcut]) + sum(UBInt33_1) + 2*sum(UBInt33_2)) * 1e10
            ub_amuNLO88 = (alpha/pi)^3 * (sum(int88[1:tcut,1:tcut]) + sum(UBInt88_1) + 2*sum(UBInt88_2)) * 1e10

            push!(lb33, lb_amuNLO33)
            push!(lb88, lb_amuNLO88)
            push!(ub33, ub_amuNLO33)
            push!(ub88, ub_amuNLO88)

            if LBOUND_IMPR
                if aens*(tcut-1) < tcut_fix  # we fix the eff energy at some point
                    Eeff33=Eeff(tcut, obs33)
                    Eeff88=Eeff(tcut, obs88R)
                end

                LB33 = Gbounding(t, tcut, obs33, ens, Eeff33)
                LB88 = Gbounding(t, tcut, obs88R, ens, Eeff88)

                LBint33_1 = (LB33 .* hcat(LB33...)) .* TMRc[tcut+1:end,tcut+1:end]
                LBint88_1 = (LB88 .* hcat(LB88...)) .* TMRc[tcut+1:end,tcut+1:end]
                LBint33_2 = (obs33[1:tcut] .* hcat(LB33...)) .* TMRc[1:tcut,tcut+1:end]
                LBint88_2 = (obs88[1:tcut] .* hcat(LB88...)) .* TMRc[1:tcut,tcut+1:end]
                
                lb_impr_amuNLO33 = (alpha/pi)^3 * (sum(int33[1:tcut,1:tcut]) + sum(LBint33_1) + 2*sum(LBint33_2)) * 1e10
                lb_impr_amuNLO88 = (alpha/pi)^3 * (sum(int88[1:tcut,1:tcut]) + sum(LBint88_1) + 2*sum(LBint88_2)) * 1e10

                push!(lb_impr33, lb_impr_amuNLO33)
                push!(lb_impr88, lb_impr_amuNLO88)
            end
        end
    end
    ub = [ub33,ub88]
    lb = [lb33,lb88]
    LBOUND_IMPR ? lb_impr = [lb_impr33,lb_impr88] : nothing

    LB = LBOUND_IMPR ? lb_impr : lb

    for i in collect(1:2)
        component = ["33","88"]

        println("      - Applying Bounding Method for component $(component[i])...")

        averb = (ub[i].+lb[i])./2; uwerr.(averb)
        x0    = findfirst(abs.(value.(ub[i]).-value.(lb[i])) .< 0.75.*err.(averb))
        xend_x0 = findfirst(abs.(averb[x0:end].-averb[x0]) .> 0.5*err(averb[x0]))
        xend_ = isnothing(xend_x0) ? x0+15  : x0-1 + xend_x0
        if xend_-x0 > 5
            xend = xend_>length(averb) ? length(averb) : xend_
        else
            xend=x0+3
            x0 > 3 ? x0-=3 : x0=1
        end
        # print("      - x = $x0 and xend = $xend")

        amu = sum(averb[x0:xend])/length(averb[x0:xend]); uwerr(amu)
        # amu = sum(ub[x0:xend]+lb[x0:xend])/(2*length(averb[x0:xend])); uwerr(amu)

        aux1 = sum(averb[x0:xend].^2)/length(averb[x0:xend])
        aux2 = amu^2
        syst = sqrt(abs(value(aux1 - aux2)))

        println("      - Plotting...")

        label = ["BM soultion",L"Upper bound ($E_0$)","Lower bound (zero)"]

        plateau_fm = value(aens).*(collect(x0:xend).+tcut0.-2)
        fill_between(plateau_fm, value(amu)-sqrt(err(amu)^2+syst^2), value(amu)+sqrt(err(amu)^2+syst^2), alpha=0.4, color="gray", label="BM soultion")

        uwerr.(ub[i]); uwerr.(lb[i])
        errorbar(tcut_fm, value.(ub[i]), err.(ub[i]), fmt="o", mfc="none", label=label[1], color="red", capsize=2)
        errorbar(tcut_fm, value.(lb[i]), err.(lb[i]), fmt="d", mfc="none", label=label[2], color="green", capsize=2)
        if LBOUND_IMPR
            uwerr.(lb_impr[i])
            errorbar(tcut_fm, value.(lb_impr[i]), err.(lb_impr[i]), fmt="d", mfc="none", label=L"Lower bound ($E_{eff}$)", color="limegreen", capsize=2)
            push!(label, L"Lower bound ($E_{eff}$)")
        end
        errorbar(tcut_fm, value.(averb), err.(averb), fmt="s", mfc="none", label="Average", color="purple", capsize=2)
        push!(label, "Average")
            
        axis("tight")
        PyPlot.title(ens.id*": Bounding Method (impr. set $impr_set; discr. $discr)")
        xlabel(L"$t_{cut}$ [fm]")
        ylabel(latexstring("a_\\mu^{\\rm{HVP}}[\\rm{NLO}_$diag^{$(component[i])}]"))
        ylim((value(amu)-6*err(amu), value(amu)+6*err(amu)))
        legend(label, loc="best")
        display(gcf())      #display the figure
        close("all")

        println("      ⇒ Result = $(value(amu)) ± $(err(amu)) ± $syst.")
        # println("$(value.(ub[i][x0-6:x0]))\n")
        # println("$(value.(lb_impr[i][x0-6:x0]))\n")
        # println("$(err.(averb[x0-6:x0]))\n")
    end
end


#-- Plot 2.4: Bounding Method (33 & 88) [c]

@info("Plot 2.4: Bounding Method (33 & 88) [c]")

tcut0 = 10
tstep = 1

tcut_fix = 1.2  # fm

aens = t0_ph / sqrt(t0)
tcut_fm =  value.(aens .* (collect(tcut0:tstep:t[end-1]).-1))

obs88 = obs88R

int3333 = (obs33[t].*hcat(obs33[t]...)) .* TMRc
int8888 = (obs88[t].*hcat(obs88[t]...)) .* TMRc
int3388 = (obs33[t].*hcat(obs88[t]...)) .* TMRc

ub3333 = Vector{uwreal}(); lb3333 = Vector{uwreal}()
ub8888 = Vector{uwreal}(); lb8888 = Vector{uwreal}()
ub3388 = Vector{uwreal}(); lb3388 = Vector{uwreal}()

println("      - Computing bounded estimations...")

Eeff33 = uwreal(0.0); Eeff88 = uwreal(0.0)
if ens.kappa_l == ens.kappa_s
    for tcut in ProgressBar(tcut0:tstep:t[end-1])

        if aens*tcut < tcut_fix  # we fix the eff energy at some point
            Eeff33=Eeff(tcut, obs33)
        end
        UB33 = mrho < E2pi ? Gbounding(t, tcut, obs33, ens, mrho) : Gbounding(t, tcut, obs33, ens, E2pi)
        LB33 = Gbounding(t, tcut, obs33, ens, Eeff33)


        UBInt3333_1 = (UB33 .* hcat(UB33...)) .* TMRc[tcut+1:end,tcut+1:end]
        UBInt3333_2 = (obs33[1:tcut] .* hcat(UB33...)) .* TMRc[1:tcut,tcut+1:end]

        LBInt3333_1 = (LB33 .* hcat(LB33...)) .* TMRc[tcut+1:end,tcut+1:end]
        LBInt3333_2 = (obs33[1:tcut] .* hcat(LB33...)) .* TMRc[1:tcut,tcut+1:end]

        ub_amuNLO3333 = (alpha/pi)^3 * (sum(int3333[1:tcut,1:tcut]) + sum(UBInt3333_1) + 2*sum(UBInt3333_2)) * 1e10

        lb_amuNLO3333 = (alpha/pi)^3 * (sum(int3333[1:tcut,1:tcut]) + sum(LBInt3333_1) + 2*sum(LBInt3333_2)) * 1e10

        push!(ub3333, ub_amuNLO3333); push!(lb3333, lb_amuNLO3333)
        push!(ub8888, ub_amuNLO3333); push!(lb8888, lb_amuNLO3333)
        push!(ub3388, ub_amuNLO3333); push!(lb3388, lb_amuNLO3333)
    end
else
    for tcut in ProgressBar(tcut0:tstep:t[end-1])

        if aens*tcut < tcut_fix  # we fix the eff energy at some point
            Eeff33=Eeff(tcut, obs33)
            Eeff88=Eeff(tcut, obs88)
        end
        UB33 = mrho < E2pi ? Gbounding(t, tcut, obs33, ens, mrho) : Gbounding(t, tcut, obs33, ens, E2pi)
        LB33 = Gbounding(t, tcut, obs33, ens, Eeff33)
        UB88 = mrho < E3pi ? Gbounding(t, tcut, obs88, ens, mrho) : Gbounding(t, tcut, obs88, ens, E3pi)
        LB88 = Gbounding(t, tcut, obs88, ens, Eeff88)

        UBInt3333_1 = (UB33 .* hcat(UB33...)) .* TMRc[tcut+1:end,tcut+1:end]
        UBInt8888_1 = (UB88 .* hcat(UB88...)) .* TMRc[tcut+1:end,tcut+1:end]
        UBInt3388_1 = (UB33 .* hcat(UB88...)) .* TMRc[tcut+1:end,tcut+1:end]
        UBInt3333_2 = (obs33[1:tcut] .* hcat(UB33...)) .* TMRc[1:tcut,tcut+1:end]
        UBInt8888_2 = (obs88[1:tcut] .* hcat(UB88...)) .* TMRc[1:tcut,tcut+1:end]
        UBInt3388_2 = (obs33[1:tcut] .* hcat(UB88...)) .* TMRc[1:tcut,tcut+1:end]
        UBInt8833_2 = (obs88[1:tcut] .* hcat(UB33...)) .* TMRc[1:tcut,tcut+1:end]

        LBInt3333_1 = (LB33 .* hcat(LB33...)) .* TMRc[tcut+1:end,tcut+1:end]
        LBInt8888_1 = (LB88 .* hcat(LB88...)) .* TMRc[tcut+1:end,tcut+1:end]
        LBInt3388_1 = (LB33 .* hcat(LB88...)) .* TMRc[tcut+1:end,tcut+1:end]
        LBInt3333_2 = (obs33[1:tcut] .* hcat(LB33...)) .* TMRc[1:tcut,tcut+1:end]
        LBInt8888_2 = (obs88[1:tcut] .* hcat(LB88...)) .* TMRc[1:tcut,tcut+1:end]
        LBInt3388_2 = (obs33[1:tcut] .* hcat(LB88...)) .* TMRc[1:tcut,tcut+1:end]
        LBInt8833_2 = (obs88[1:tcut] .* hcat(LB33...)) .* TMRc[1:tcut,tcut+1:end]

        ub_amuNLO3333 = (alpha/pi)^3 * (sum(int3333[1:tcut,1:tcut]) + sum(UBInt3333_1) + 2*sum(UBInt3333_2)) * 1e10
        ub_amuNLO8888 = (alpha/pi)^3 * (sum(int8888[1:tcut,1:tcut]) + sum(UBInt8888_1) + 2*sum(UBInt8888_2)) * 1e10
        ub_amuNLO3388 = (alpha/pi)^3 * (sum(int3388[1:tcut,1:tcut]) + sum(UBInt3388_1) + sum(UBInt3388_2) + sum(UBInt8833_2)) * 1e10

        lb_amuNLO3333 = (alpha/pi)^3 * (sum(int3333[1:tcut,1:tcut]) + sum(LBInt3333_1) + 2*sum(LBInt3333_2)) * 1e10
        lb_amuNLO8888 = (alpha/pi)^3 * (sum(int8888[1:tcut,1:tcut]) + sum(LBInt8888_1) + 2*sum(LBInt8888_2)) * 1e10
        lb_amuNLO3388 = (alpha/pi)^3 * (sum(int3388[1:tcut,1:tcut]) + sum(LBInt3388_1) + sum(LBInt3388_2) + sum(LBInt8833_2)) * 1e10

        push!(ub3333, ub_amuNLO3333); push!(lb3333, lb_amuNLO3333)
        push!(ub8888, ub_amuNLO8888); push!(lb8888, lb_amuNLO8888)
        push!(ub3388, ub_amuNLO3388); push!(lb3388, lb_amuNLO3388)
    end
end

ub = [ub3333,ub8888,ub3388]
lb = [lb3333,lb8888,lb3388]

for (i,comp) in enumerate(["3333","8888","3388"])
    println("      - Applying Bounding Method for component $(comp)...")

    averb = (ub[i].+lb[i])./2; uwerr.(averb)
    x0    = findfirst(abs.(value.(ub[i]).-value.(lb[i])) .< 0.75.*err.(averb))
    xend_x0 = findfirst(abs.(averb[x0:end].-averb[x0]) .> 0.5*err(averb[x0]))
    xend_ = isnothing(xend_x0) ? x0+15  : x0-1 + xend_x0
    if xend_-x0 > 5
        xend = xend_>length(averb) ? length(averb) : xend_
    else
        xend=x0+3
        x0 > 3 ? x0-=3 : x0=1
    end
    plateau_fm = value(aens).*(collect(x0:xend).+tcut0.-2)
    amu = sum(averb[x0:xend])/length(averb[x0:xend])

    aux1 = sum(averb[x0:xend].^2)/length(averb[x0:xend])
    aux2 = amu^2
    syst = sqrt(abs(value(aux1 - aux2)))

    tmax = length(ub[i])
    t_ = 1:tmax

    uwerr.(ub[i]); uwerr.(lb[i])
    uwerr(amu)

    label = ["BM soultion",L"Upper bound ($E_0$)",L"Lower bound ($E_{eff}$)","Average"]
    fill_between(plateau_fm, value(amu)-sqrt(err(amu)^2+syst^2), value(amu)+sqrt(err(amu)^2+syst^2), alpha=0.4, color="gray", label=label[1])
    errorbar(tcut_fm, value.(ub[i][t_]), err.(ub[i][t_]), fmt="o", mfc="none", label=label[1], color="red", capsize=2)
    errorbar(tcut_fm, value.(lb[i][t_]), err.(lb[i][t_]), fmt="d", mfc="none", label=label[2], color="limegreen", capsize=2)
    errorbar(tcut_fm, value.(averb[t_]), err.(averb[t_]), fmt="s", mfc="none", label=label[3], color="purple", capsize=2)
        
    axis("tight")
    PyPlot.title(ens.id*": Bounding Method (impr. set $impr_set; discr. $discr)")
    xlabel(L"$t_{cut}$ [fm]")
    ylabel(latexstring("a_\\mu^{\\rm{HVP}}[\\rm{NLO}_c^{$(comp[1:2]),$(comp[3:4])}]"))
    legend(label, loc="best")
    display(gcf())      #display the figure
    close("all")
end


#-- Plot 2.5: Short distance cut-off (cc_disc and c8 components)

@info("Plot 2.5: Short distance cut-off (cc_disc and c8 components)")

OBS = "c8"
DIAG = ["b"]

println("- Computing for corr. '$OBS' and diag. '$DIAG'")


function windowfunc(t::Union{Int64,Float64,uwreal},tstar::Union{Int64,Float64},Delta::Union{Int64,Float64})
    return 0.5*(1 + tanh((t-tstar)/Delta))
end

myobsC = Dict{String, Vector{uwreal}}("cc" => obsCCdisc, "c8" => obsC8)[OBS][t]

aens = t0_ph / sqrt(t0)
tcut_fm =  value.(aens .* t)

tstar = findfirst(abs.(err.(myobsC[2:end])./value.(myobsC[2:end])) .> 0.5) + 1
Delta = 1.5

for diag in DIAG
    println("   - Starting diagram $diag")

    println("      - Computing some things...")
    if diag in ["a", "b"]
        windvec = 1 .- windowfunc.(t,tstar,Delta)

        TMR = diag == "a" ? TMRa : TMRb 
        int = myobsC .* TMR
        intW = myobsC .* windvec .* TMR

    elseif diag == "c"
        int = (myobsC .* hcat(myobsC...)) .* TMRc

    end
    uwerr.(int); uwerr.(intW)

    println("      - Plotting...")

    errorbar(tcut_fm, value.(int), err.(int), fmt="^", mfc="none", label=latexstring("K(t)G^{$OBS}(t)"), color="gray", capsize=2)
    errorbar(tcut_fm, value.(intW), err.(intW), fmt="o", mfc="none", label=latexstring("K(t)W(t)G^{$OBS}(t)"), color="blue", capsize=2)
        
    axis("tight")
    OBS == "cc" ? PyPlot.title(ens.id*": Charm-charm disc. ("*diag*")") : PyPlot.title(ens.id*": Charm-8 disc. ("*diag*")")
    xlabel(L"$t_{cut}$ [fm]")
    # OBS == "cc" ?  ylim((-5e-6, 15e-6)) : ylim((-15e-7, 5e-7))
    max_=maximum(value.(intW)); min_=minimum(value.(intW))
    index = argmax(abs.([max_,min_]))
    ylim_ = index == 1 ? [-0.5*max_,1.5*max_] : [1.5*min_,-0.5*min_]
    ylim(ylim_)
    legend([latexstring("K(t)G^{$OBS}(t)"),latexstring("K(t)W(t)G^{$OBS}(t)")], loc="best")
    display(gcf())      #display the figure
    close("all")

    amu = (alpha/pi)^3 * sum(intW) * 1e10; uwerr(amu)
    println("      ⇒ Result = $(value(amu)) ± $(err(amu)).")
end



##===========================================================================##
##==========================> (3) HVP & FVC Plots <==========================##
##===========================================================================##

#-- Plot 3.1: kappa_C correction

@info("Plot 3.1: kappa C correction")

SYST = true

obs = [HVP["HVP"]["gcc_lc_conn_sim"],HVP["HVP"]["gcc_lc_conn_sim+"]]
kappaC = [kcd_in[ensid]["kappaC_sim"],kcd_in[ensid]["kappaC_sim_plus"]]
# kappaC_tar = uwreal([kcd_in[ens.id]["kappaC"],kcd_in[ens.id]["kappaC_err"]], "kappaC target")
fb = BDIO_open(joinpath(path_bdio,"kappaC_tar",ensid,"$(ensid)_kappaC"),"r")
kappaC_tar = uwreal(0.0)
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    kappaC_tar = ALPHAdobs_read_next(fb)
end


@. lin_model(x,p) = p[1] + p[2] * x

fit  = fit_routine(lin_model, kappaC, obs, 2)
par = fit.param
obs_tar = lin_model(kappaC_tar, par)[1]

uwerr(kappaC_tar)
uwerr(obs_tar)
obs_val = value.(obs)
obs_err = err.(obs)

minimum(vcat(value(kappaC_tar),kappaC...))

figure()
PyPlot.title("$(ensid)")
errorbar(kappaC, obs_val, obs_err, fmt="d", capsize=2, mfc="none", color="black")
errorbar(value(kappaC_tar), value(obs_tar), xerr = err(kappaC_tar), yerr=err(obs_tar), capsize=2, color="red", fmt="o")
kappaC_arr = vcat(value(kappaC_tar)+err(kappaC_tar),value(kappaC_tar)-err(kappaC_tar),kappaC...)
xarr = Float64.(range(minimum(kappaC_arr), maximum(kappaC_arr), length=100))
yarr = lin_model(xarr, par); uwerr.(yarr)
fill_between(xarr, value.(yarr) .- err.(yarr), value.(yarr) .+ err.(yarr), color="royalblue", alpha=0.2)
if SYST
    if !(kappaC[1] <= value(kappaC_tar) <= kappaC[2])
        if value(kappaC_tar) < kappaC[1]
            extra = 1.5 .* (abs.(xarr.-kappaC[1]) ./ (kappaC[2]-kappaC[1])).^2 .* err(obs[1]) .* (xarr.<kappaC[1])
            syst  = 1.5 * (abs(value(kappaC_tar)-kappaC[1]) / (kappaC[2]-kappaC[1]))^2 * err(obs[1])
        elseif value(kappaC_tar) > kappaC[2]
            extra = 1.5 .* (abs.(xarr.-kappaC[2]) ./ (kappaC[2]-kappaC[1])).^2 .* err(obs[2]) .* (xarr.>kappaC[2])
            syst  = 1.5 * (abs(value(kappaC_tar)-kappaC[2]) / (kappaC[2]-kappaC[1]))^2 * err(obs[2])
        end
        fill_between(xarr, value.(yarr) .- err.(yarr) .- extra, value.(yarr) .- err.(yarr), color="red", alpha=0.05)
        fill_between(xarr, value.(yarr) .+ err.(yarr), value.(yarr) .+ err.(yarr) .+ extra, color="red", alpha=0.05)
        errorbar(value(kappaC_tar), value(obs_tar), xerr = err(kappaC_tar), yerr=sqrt(err(obs_tar)^2+syst^2), capsize=2, color="red", fmt="o")
    end
end
ylabel(L"$a_\mu^{\rm{hvp}}[\rm{NLO}^{\rm{cc-conn}}]$")
xlabel(L"$\kappa_c$")
tight_layout()
display(gcf())
PyPlot.savefig(joinpath(julia_script_directory,"..","Slides & Plots","Plots kappaC corr",ens.id))
close()


##===============================================================================##
##==========================> (4) Extrapolation plots <==========================##
##===============================================================================##

#-- Plot 4.1: Continiuum projection plots

fitres_len = length(fitres[IMPR_SET[1]][mykeys[1]])
key_len    = length(mykeys)

@info("Computing 'y' points for continuum plots...")

xarr = [Float64.(range(0.0, 1.5*maximum(value.(xydata["xdata"][:,1])), length=100)) fill(value(phi2_ph), 100) fill(value(phi4_ph), 100)]


yarr = Dict{String,Any}()
for impr_set in IMPR_SET
    yarr[impr_set] = Dict{String,Any}()
    for key in mykeys
        yarr[impr_set][key] = Vector{Vector{uwreal}}()
        println("- impr. set $impr_set; comp = $key:")
        for i in ProgressBar(1:fitres_len)
            my = f_tot_isov[i](xarr, param[impr_set][key][i]); uwerr.(my)
            push!(yarr[impr_set][key], my)
        end
    end
end

#- 

@info("Plot 4.1.1: Continiuum 'magic' plot")

PROJpoints = true

# color_list = ["blue","purple","red","orange","brown","gray"]
color_list = ["blue","red","green","brown"]
# color_list = ["blue","green","firebrick","orange","brown","gray"]
# color_list = ["red","green"]

fmt_list = ["^","o","s","d"]


label_list = []

fig = figure(figsize=(10,7.5))

println("- Plotting...")
weights = MA["weight_tot"]
for (j,impr_set) in enumerate(IMPR_SET)
    for (k,key) in enumerate(mykeys)
        for i in collect(1:fitres_len)
            fill_between(xarr[:,1], value.(yarr[impr_set][key][i]).-err.(yarr[impr_set][key][i]), value.(yarr[impr_set][key][i]).+err.(yarr[impr_set][key][i]), alpha=weights[i+fitres_len*(key_len*(j-1)+k-1)], color=color_list[key_len*(j-1)+k])
        end
        if PROJpoints
            println("   - [Computing projection for $key and Set $impr_set ...]")
            xproj = xydata["xdata"]
            xproj[:,2] = fill(phi2_ph, length(xydata["xdata"][:,1]))::Vector{uwreal}
            xproj[:,3] = fill(phi4_ph, length(xydata["xdata"][:,1]))::Vector{uwreal}
            yproj_list = []
            for i in collect(1:fitres_len)
                push!(yproj_list, f_tot_isov[i](xproj,param[impr_set][key][i]))
            end

            w_ = weights[1+fitres_len*(key_len*(j-1)+k-1):fitres_len*(key_len*(j-1)+k)]
            w = w_./sum(w_)
            
            yproj = []
            yproj_syst = []
            for arg in collect(1:length(xproj[:,1]))
                y_ = [element[arg] for element in yproj_list]
                val, syst = model_average(y_,w)
                push!(yproj,val)
                push!(yproj_syst,syst)
            end
            uwerr.(yproj)

            errorbar(value.(xproj[:,1]), value.(yproj), err.(yproj), fmt=fmt_list[key_len*(j-1)+k], capsize=2, color=color_list[key_len*(j-1)+k], mfc="none", label="discr. $(key[end-1:end]); set $impr_set")
            errorbar(value.(xproj[:,1]), value.(yproj), sqrt.(err.(yproj).^2 .+ yproj_syst.^2), fmt=fmt_list[key_len*(j-1)+k], capsize=2, color=color_list[key_len*(j-1)+k], mfc="none")
            push!(label_list,"discr. $(key[end-1:end]); set $impr_set")
        end
    end
end
uwerr(MA["amu"])
errorbar([0.0],[value(MA["amu"])],[err(MA["amu"])],fmt="o",mfc="none",color="black", ms=5, capsize=3)
errorbar([0.0],[value(MA["amu"])],[sqrt(err(MA["amu"])^2+(MA["syst"])^2)],fmt="o",color="black", ms=5, capsize=3)
PyPlot.title("Projection to continuum extrapolation")
axvline(ls="dashed", color="black", lw=0.2, alpha=0.7) 
xlabel(L"$a^2/8t_0$")
xlim(-0.002,0.06)
ylim(-9,-7)
if diag == "NLOa&b"
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{NLOa}\\&\\rm{b}^{\\mathrm{$comp}}]"))
else
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{$diag}^{\\mathrm{$comp}}]"))
end
PROJpoints ? legend(label_list, loc="best") : nothing
tight_layout()
display(gcf())
close()

##


label_list = []

println("- Plotting...")
weights = MA["weight_tot"]

for (j,impr_set) in enumerate(["1","2"])
    for (k,key) in enumerate(mykeys)
        if PROJpoints
            println("   - [Computing projection for $key and Set $impr_set ...]")
            xproj = xydata["xdata"]
            xproj[:,2] = fill(phi2_ph, 12)::Vector{uwreal}
            xproj[:,3] = fill(phi4_ph, 12)::Vector{uwreal}
            yproj_list = []
            for i in collect(1:fitres_len)
                push!(yproj_list, f_tot_isov[i](xproj,param[impr_set][key][i]))
            end

            w_ = weights[1+fitres_len*(key_len*(j-1)+k-1):fitres_len*(key_len*(j-1)+k)]
            w = w_./sum(w_)
            
            yproj = []
            yproj_syst = []
            for arg in collect(1:length(xproj[:,1]))
                y_ = [element[arg] for element in yproj_list]
                val, syst = model_average(y_,w)
                push!(yproj,val)
                push!(yproj_syst,syst)
            end
            uwerr.(yproj)

            errorbar(value.(xproj[:,1]), value.(yproj), err.(yproj), fmt=fmt_list[key_len*(j-1)+k], capsize=2, color=color_list[key_len*(j-1)+k], mfc="none", label="discr. $(key[end-1:end]); set $impr_set")
            push!(label_list,"discr. $(key[end-1:end]); set $impr_set")
        end
    end
end
uwerr(MA["amu"])
errorbar([0.0],[value(MA["amu"])],[err(MA["amu"])],fmt="o",mfc="none",color="black", ms=5, capsize=3)
errorbar([0.0],[value(MA["amu"])],[sqrt(err(MA["amu"])^2+(MA["syst"])^2)],fmt="o",color="black", ms=5, capsize=3)
PyPlot.title("Projection to continuum extrapolation")

xlabel(L"$a^2/8t_0$")
xlim(-0.002,0.06)
ylim(-9,-7)
if diag == "NLOa&b"
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{NLOa}\\&\\rm{b}^{\\mathrm{$comp}}]"))
else
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{$diag}^{\\mathrm{$comp}}]"))
end
PROJpoints ? legend(label_list, loc="best") : nothing

# legend(handles=handles)
tight_layout()
display(gcf())
close()


#-

@info("Plot 4.1,2: Best fits continuum plot")

println("- Searching for best fits...")

toprange = 5

maxWeight = sortperm(MA["weight_tot"], rev=true)

toparg = Dict{String,Any}()
toparg = Dict{String,Any}()
for (j,impr_set) in enumerate(IMPR_SET)
    toparg[impr_set] = Dict{String,Any}()
    for (k,key) in enumerate(mykeys)
        toparg[impr_set][key] = filter(x -> (k+2*j-3)*fitres_len < x < (k+2*j-2)*fitres_len, maxWeight)[1:toprange] .- (k+2*j-3)*fitres_len
    end
end


color_list = ["blue","purple","red","orange","brown","gray"]

println("- Plotting...")

label_list = []
weights = MA["weight_tot"]
for (j,impr_set) in enumerate(IMPR_SET)
    for (k,key) in enumerate(mykeys)
        for (i,arg) in enumerate(toparg[impr_set][key])
            if i == 1
                PyPlot.plot(xarr[:,1], value.(yarr[impr_set][key][arg]), color=color_list[2*(j-1)+k], alpha=0.4,label = "Set $impr_set; $key")
                push!(label_list,"Set $impr_set; $key")
            else
                PyPlot.plot(xarr[:,1], value.(yarr[impr_set][key][arg]), color=color_list[2*(j-1)+k], alpha=0.4)
            end
        end
    end
end
errorbar([0.0],[value(MA["amu"])],[sqrt(err(MA["amu"])^2+(MA["syst"])^2)],fmt="o",mfc="none",color="black")
PyPlot.title("Projection to continuum extrapolation")
axvline(ls="dashed", color="black", lw=0.2, alpha=0.7) 
xlabel(L"$a^2/8t_0$")
if diag == "NLOa&b"
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{NLOa}\\&\\rm{b}^{\\mathrm{$comp}}]"))
else
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{$diag}^{\\mathrm{$comp}}]"))
end
# legend(label_list, loc="best")
tight_layout()
display(gcf())
close()


#-- Plot 4.2: Chiral projection plot

@info("Plot 4.2: Chiral projection plot")

fmt = ["o","^","s","d","*"]
color = ["orange","red","purple","blue","green"]

length_model = myinfo["length"]

for (i,impr_set) in enumerate(IMPR_SET) 
    for (j,key) in enumerate(mykeys)
        label = []
        firstarg, lastarg = [(key_len*(i-1)+j-1)*length_model+1,(key_len*(i-1)+j)*length_model]
        arg_max = argmax(MA["weight_tot"][firstarg:lastarg])
        discr = key[5:6]
        fit = fitres[impr_set][key][arg_max]
        model = label_tot_isov[arg_max]
        myres   = res[impr_set][key][arg_max]; uwerr(myres)
        myparam = param[impr_set][key][arg_max]; uwerr.(myparam)

        # ydata = xydata["ydata"][impr_set][key]; uwerr.(ydata) # ¿is the extrapolation for phi4 still required here? YES
        xproj = xydata["xdata"]
        xproj[:,3] = fill(phi4_ph, length(xydata["xdata"][:,1]))::Vector{uwreal}
        ydata = f_tot_isov[arg_max](xproj,myparam); uwerr.(ydata)

        println("Chosen model (diag=$diag, set=$impr_set, discr=$discr): $(label_tot_isov[arg_max]) (n: $arg_max)")
        println("- chi2/dof     = $(fit.chi2/fit.dof)")
        println("- chi2/chi2exp = $(fit.chi2/fit.chi2exp)")

        for (k,b) in  enumerate(sort(unique(getfield.(ensInfo, :beta))))
            push!(label,L"$\beta=$"*"$b")
            n_ = findall(x->x.beta == b, ensInfo)
            a2_aux = mean(value.(xydata["xdata"][:,1][n_]))
            errorbar(value.(xydata["xdata"][n_,2]), value.(ydata[n_]), err.(ydata[n_]), fmt=fmt[k], capsize=2, color=color[k], mfc="none", label=label[k])
            xxx = [fill(a2_aux,100) Float64.(range(0.02, 0.8, length=100)) fill(value(phi4_ph), 100)]
            yyy = f_tot_isov[arg_max](xxx, myparam)
            PyPlot.plot(xxx[:,2],value.(yyy),ls="--",color=color[k],lw=0.5)
        end
        xxx_ph = [fill(0.0,100) Float64.(range(0.02, 0.8, length=100)) fill(value(phi4_ph), 100)]
        yyy_ph = f_tot_isov[arg_max](xxx_ph, myparam); uwerr.(yyy_ph)
        errorbar(value(phi2_ph), value(myres), err(myres), fmt="^", capsize=2, color="black",label="ph.")
        push!(label,"ph.")
        PyPlot.plot(xxx_ph[:,2],value.(yyy_ph),ls="--",color="gray",lw=0.5)
        fill_between(xxx_ph[:,2], value.(yyy_ph)-err.(yyy_ph), value.(yyy_ph)+err.(yyy_ph), alpha=0.2, color="gray")
        PyPlot.title("Chiral and continuum extrapolation (Discr. $(key[end-1:end]); impr. set $impr_set)")
        axvline(x=value(phi2_ph), ls="dashed", color="black", lw=0.2, alpha=0.7) 
        xlabel(L"$\Phi_2$")
        if diag == "NLOa&b"
            ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{NLOa}\\&\\rm{b}^{\\mathrm{$comp}}]"))
        else
            ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{$diag}^{\\mathrm{$comp}}]"))
        end
        legend(label,loc="best")
        tight_layout()
        display(gcf())
        close()
    end
end


#-- Plot 4.3: Model Average analysis

DISCARD = true


@info("Plot 4.3: Model Average analysis")

println("- Starting diag. $diag")

p1 = MA["res_tot"]; uwerr.(p1)
uwerr.(MA["amu"]); v =  value(MA["amu"])
e = err(MA["amu"]); esyst = sqrt(err(MA["amu"])^2+MA["syst"]^2)
w = MA["weight_tot"]

chi2chi2exp=[]; chi2dof=[]; pval=[]; mods=[]
for impr_set in ["1","2"]
    for key in mykeys
        push!(chi2chi2exp, getfield.(fitres[impr_set][key],:chi2)./getfield.(fitres[impr_set][key],:chi2exp))
        push!(chi2dof    , getfield.(fitres[impr_set][key],:chi2)./getfield.(fitres[impr_set][key],:dof    ))
        push!(pval       , getfield.(fitres[impr_set][key],:pval))
        push!(mods, myinfo["label_tot_isov"])
    end
end

chi2chi2exp = vcat(chi2chi2exp...)
chi2dof     = vcat(chi2dof...)
pval        = vcat(pval...)
mods        = vcat(mods...)



fig = figure(figsize=(15,10))

subplots_adjust(hspace=0.1) 
subplot(411)    
ax1 = gca()                
x = 1:length(p1)
y = value.(p1)
dy = err.(p1)

# PyPlot.title("Model Average (diag. $diag; comp. $comp)")
PyPlot.title("Model Average")

length_model = length(f_tot_isov)
length_keys  = length(mykeys)
if comp in ["cc disc","c8 disc"]
    fill_between(x[1:Int64(length(x)/2)], v-e, v+e, color="limegreen", alpha=0.2)
    fill_between(x[1:Int64(length(x)/2)], v-esyst, v+esyst, color="limegreen", alpha=0.2)
    color_list = ["dodgerblue","firebrick"]
    IMPR_SET = ["1"]
else
    fill_between(x, v-e, v+e, color="limegreen", alpha=0.2)
    fill_between(x, v-esyst, v+esyst, color="limegreen", alpha=0.2)
    # color_list = ["blue","purple","red","orange","brown","gray"]
    color_list = ["blue","red","green","brown"]
    IMPR_SET = ["1","2"]
end
for (i,impr_set) in enumerate(IMPR_SET)
    for (j,key) in enumerate(mykeys)
        o,f=[(length_keys*(i-1)+j-1)*length_model+1,(length_keys*(i-1)+j)*length_model]
        errorbar(x[o:f], y[o:f], dy[o:f], fmt="o", mfc="none", color=color_list[length_keys*(i-1)+j], ms=10, capsize=2)
    end
end
setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels
if diag == "NLOa&b"
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{NLOa}\\&\\rm{b}^{\\mathrm{$comp}}]"))
else
    ylabel(latexstring("a_{\\mu}^{\\rm{HVP}}[\\rm{$diag}^{\\mathrm{$comp}}]"))
end
ylim([v-6*e,v+6*e])

subplot(412)
ax2=gca()
for (i,impr_set) in enumerate(IMPR_SET)
    for (j,key) in enumerate(mykeys)
        o,f=[(length_keys*(i-1)+j-1)*length_model+1,(length_keys*(i-1)+j)*length_model]
        PyPlot.bar(x[o:f], w[o:f], alpha=0.4, color=color_list[length_keys*(i-1)+j], edgecolor="black", linewidth=1.0)
    end
end
# PyPlot.bar(x, w, alpha=0.4, color="royalblue", edgecolor="blue", linewidth=1.5)
setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels
# errorbar(mods, weight_model, 0*dy, color="green")
ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))

subplot(413)
ax3=gca()
for (i,impr_set) in enumerate(IMPR_SET)
    for (j,key) in enumerate(mykeys)
        o,f=[(length_keys*(i-1)+j-1)*length_model+1,(length_keys*(i-1)+j)*length_model]
        PyPlot.bar(x[o:f], chi2chi2exp[o:f], alpha=0.4, color=color_list[length_keys*(i-1)+j], edgecolor="black", linewidth=1.0)
    end
end
setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels
ylabel(L"$\chi^2/\chi^2_{\mathrm{exp}}$")
ylim([0.0,2.0])

# subplot(414)
# for (i,impr_set) in enumerate(IMPR_SET)
#     for (j,key) in enumerate(mykeys)
#         o,f=[(length_keys*(i-1)+j-1)*length_model+1,(length_keys*(i-1)+j)*length_model]
#         PyPlot.bar(x[o:f], chi2dof[o:f], alpha=0.4, color=color_list[length_keys*(i-1)+j], edgecolor="black", linewidth=1.0)
#     end
# end
# ylabel(L"$\chi^2/\rm{dof}$")
# ylim([0.0,2.0])

subplot(414)
for (i,impr_set) in enumerate(IMPR_SET)
    for (j,key) in enumerate(mykeys)
        o,f=[(length_keys*(i-1)+j-1)*length_model+1,(length_keys*(i-1)+j)*length_model]
        PyPlot.bar(x[o:f], pval[o:f], alpha=0.4, color=color_list[length_keys*(i-1)+j], edgecolor="black", linewidth=1.0)
    end
end
ylabel(L"$\rm{p-values}$")
ylim([0.0,1.0])

# PyPlot.xticks(x, mods, rotation=90)
xlabel("Models [n = $length_model]")
tight_layout()
display(fig)
close("all")

#--

top = 12
diag = "a"

allmods = false; n = 0 

if allmods
    mods = 1:length(MA[diag]["res_tot"])
    label_tot_isov_ = []
    for _=1:Int64(length(MA[diag]["res_tot"])/length(label_tot_isov))
        push!(label_tot_isov_,label_tot_isov)
    end
    label_tot_isov_ = hcat(label_tot_isov_...)
elseif allmods == false
    mods = n*length(f_tot_isov)+1:(n+1)*length(f_tot_isov)
    label_tot_isov_ = label_tot_isov
    if n==0
        impr_set = "1"; discr = "ll"
    elseif n==1 
        impr_set = "1"; discr = "lc"
    elseif n==2
        impr_set = "2"; discr = "ll"
    elseif n==3
        impr_set = "2"; discr = "lc"
    end
end

@info("Print 3.3: Best & worst performing models for (diag. $diag; comp. $comp; impr. set $impr_set; discr. $discr)")
println("\n")


Vector{Any}(["hola","adios"])

println("- The top $top best working models are:")

bestModels   = label_tot_isov_[reverse(sortperm(MA[diag]["weight_tot"][mods]))]
bestModels_w = reverse(sort(MA[diag]["weight_tot"][mods]))

for i=1:top
    len_str = 0
    for j=1:length(bestModels[i]); len_str += length(bestModels[i][j])+4; end
    println("$(bestModels[i]) $(repeat(" ", 45-len_str))    --->   w = $(round(bestModels_w[i],digits=5))")
end
println("\n")

println("- The top $top worst working models are:")

worstModels   = reverse(bestModels)
worstModels_w = reverse(bestModels_w)

for i=1:top
    len_str = 0
    for j=1:length(worstModels[i]); len_str += length(worstModels[i][j])+4; end
    println("$(worstModels[i]) $(repeat(" ", 45-len_str))    --->   w = $(worstModels_w[i])")
end
println("\n")

println("- The top $top models most distant to the average ($(round(value(MA[diag]["amu"]),digits=2))) are:")

D       = abs.(value.(MA[diag]["res_tot"][mods]).-value(MA[diag]["amu"]))
argsD   = reverse(sortperm(D))
valueD  = value.(MA[diag]["res_tot"][mods][argsD])
wD      = MA[diag]["weight_tot"][mods][argsD]
modelsD = label_tot_isov_[argsD]

for i=1:top
    len_str = 0
    for j=1:length(modelsD[i]); len_str += length(modelsD[i][j])+4; end
    println("$(modelsD[i]) $(repeat(" ", 45-len_str))    --->   val = $(round(valueD[i],digits=2)) ;  w = $(wD[i])")
end
println("\n")

println("- The top $top models most contributing to the systematics ($(round(MA[diag]["syst"],digits=2))) are:")

systcontr  = (D./abs(value(MA[diag]["amu"]))) .* MA[diag]["weight_tot"][mods]
argssyst   = reverse(sortperm(systcontr))
valuesyst  = value.(MA[diag]["res_tot"][mods][argssyst])
wsyst      = MA[diag]["weight_tot"][argssyst]
modelssyst = label_tot_isov_[argssyst]

for i=1:top
    len_str = 0
    for j=1:length(modelssyst[i]); len_str += length(modelssyst[i][j])+4; end
    println("$(modelssyst[i]) $(repeat(" ", 45-len_str))    --->   val = $(round(valuesyst[i],digits=2)) ;  w = $(round(wsyst[i],digits=5))")
end
println("\n")

##=======================================================================##
##==========================> (5) Other plots <==========================##
##=======================================================================##

#-- Plot 5.1.1: Results vs data driven & other lattice

@info("Plot 5.1.1: Results vs data driven & other lattice")

data_res = [
[uwreal([-20.623,0.130],"res a (Jegerlehner)"),uwreal([10.349,0.063],"res b (Jegerlehner)"),uwreal([0.337,0.005],"res c (Jegerlehner)"),uwreal([-9.927,0.067],"res (Jegerlehner)")],
[uwreal([-20.77,0.08],"res a (Keshavarzi et al.)"),uwreal([10.62,0.04],"res b (Keshavarzi et al.)"),uwreal([0.34,0.01],"res c (Keshavarzi et al.)"),uwreal([-9.82,0.04],"res c (Keshavarzi et al.)")]
]

latt_res = []

# latt_res = [
# [uwreal([0.0,0.0],"nothing"),uwreal([0.0,0.0],"nothing"),uwreal([0.0,0.0],"nothing"),uwreal([-9.3,1.2],"res Fermilab Lattice")]
# ]

# latt_res = [
# [uwreal([-20.03,0.82],"res a (a=0.15fm)"),uwreal([10.37,0.11],"res b (a=0.15fm)"),uwreal([0.329,0.012],"res c (a=0.15fm)")],
# [uwreal([-19.82,0.79],"res a (a=0.12fm)"),uwreal([10.204,0.091],"res b (a=0.12fm)"),uwreal([0.321,0.011],"res c (a=0.12fm)")]
# ]

my_res   = [uwreal([-20.92,0.42],"my res a"),uwreal([10.49,0.28],"my res b"),uwreal([0.338,0.013],"my res c"),uwreal([-10.07,0.14],"my res c")]; uwerr.(my_res)
my_syst  = [0.57,0.35,0.014,0.21]

data_label = ["Jegerlehner","Keshavarzi et al."]
# latt_label = []
# latt_label = ["Fermilab 2018"]
latt_label = []
# latt_label = [L"Fermilab 2018 ($a\sim0.15$fm)",L"Fermilab 2018 ($a\sim0.12$fm)"]
my_label   = ["this work"]

yticks = vcat("",data_label,latt_label,my_label,"",""...)
yticks

data_y = collect(1:1+length(data_res))
latt_y = collect(length(data_res)+1+1:length(data_res)+length(latt_res)+1)
my_y   = collect(length(data_res)+length(latt_res)+1+1:length(data_res)+length(latt_res)+1+2)

fig = figure(figsize=(20,10))

subplots_adjust(hspace=0.1) 
subplot(541)    
ax1 = gca()
PyPlot.title(L"$a_\mu[\rm{NLO}_a]$")
data_a = [res[1] for res in data_res]; uwerr.(data_a)
latt_a = [res[1] for res in latt_res]; uwerr.(latt_a)
errorbar(value.(data_a), data_y[2:end], 0.0, xerr=err.(data_a), fmt="o", color="black", ms=10, capsize=4)
errorbar(value.(latt_a), latt_y, 0.0, xerr=err.(latt_a), fmt="o", color="gray", ms=10, capsize=4)
errorbar(value(my_res[1]), my_y[1:end-1], 0.0, xerr=err(my_res[1]), fmt="o", color="red", ms=10, capsize=4)
errorbar(value(my_res[1]), my_y[1:end-1], 0.0, xerr=sqrt(err(my_res[1])^2+my_syst[1]^2), fmt="o", color="red", ms=10, capsize=4)
y = vcat(data_y,latt_y,my_y...)
PyPlot.yticks(y, yticks, rotation=0)


vec = vcat(data_a,latt_a,my_res[1]+uwreal([0.0,my_syst[1]],"syst a")); uwerr.(vec)
trueres = value.(vec) .!= 0.0
xlim([minimum(value.(vec[trueres]).-err.(vec[trueres]))-0.1,maximum(value.(vec[trueres]).+err.(vec[trueres]))+0.1])


subplot(542)    
ax2 = gca()
PyPlot.title(L"$a_\mu[\rm{NLO}_b]$")
data_b = [res[2] for res in data_res]; uwerr.(data_b)
latt_b = [res[2] for res in latt_res]; uwerr.(latt_b)
errorbar(value.(data_b), data_y[2:end], 0.0, xerr=err.(data_b), fmt="o", color="black", ms=10, capsize=4)
errorbar(value.(latt_b), latt_y, 0.0, xerr=err.(latt_b), fmt="o", color="gray", ms=10, capsize=4)
errorbar(value(my_res[2]), my_y[1:end-1], 0.0, xerr=err(my_res[2]), fmt="o", color="red", ms=10, capsize=4)
errorbar(value(my_res[2]), my_y[1:end-1], 0.0, xerr=sqrt(err(my_res[2])^2+my_syst[2]^2), fmt="o", color="red", ms=10, capsize=4)
# setp(ax2.get_yticklabels(),visible=false) # Disable y tick labels
PyPlot.yticks(y, fill("", length(y)), rotation=0)

vec = vcat(data_b,latt_b,my_res[2]+uwreal([0.0,my_syst[2]],"syst b")); uwerr.(vec)
trueres = value.(vec) .!= 0.0
xlim([minimum(value.(vec[trueres]).-err.(vec[trueres]))-0.1,maximum(value.(vec[trueres]).+err.(vec[trueres]))+0.1])

subplot(543)    
ax3 = gca()
PyPlot.title(L"$a_\mu[\rm{NLO}_c]$")
data_c = [res[3] for res in data_res]; uwerr.(data_c)
latt_c = [res[3] for res in latt_res]; uwerr.(latt_c)
errorbar(value.(data_c), data_y[2:end], 0.0, xerr=err.(data_c), fmt="o", color="black", ms=10, capsize=2)
errorbar(value.(latt_c), latt_y, 0.0, xerr=err.(latt_c), fmt="o", color="gray", ms=10, capsize=2)
errorbar(value(my_res[3]), my_y[1:end-1], 0.0, xerr=err(my_res[3]), fmt="o", color="red", ms=10, capsize=2)
errorbar(value(my_res[3]), my_y[1:end-1], 0.0, xerr=sqrt(err(my_res[3])^2+my_syst[3]^2), fmt="o", color="red", ms=10, capsize=2)
# setp(ax2.get_yticklabels(),visible=false) # Disable y tick labels
PyPlot.yticks(y, fill("", length(y)), rotation=0)

vec = vcat(data_c,latt_c,my_res[3]+uwreal([0.0,my_syst[3]],"syst c")); uwerr.(vec)
trueres = value.(vec) .!= 0.0
xlim([minimum(value.(vec[trueres]).-err.(vec[trueres]))-0.02,maximum(value.(vec[trueres]).+err.(vec[trueres]))+0.02])

subplot(544)    
ax3 = gca()
PyPlot.title(L"$a_\mu[\rm{NLO}]$")
data_res = [res[4] for res in data_res]; uwerr.(data_res)
latt_res = [res[4] for res in latt_res]; uwerr.(latt_res)
errorbar(value.(data_res), data_y[2:end], 0.0, xerr=err.(data_res), fmt="o", color="black", ms=10, capsize=2)
errorbar(value.(latt_res), latt_y, 0.0, xerr=err.(latt_res), fmt="o", color="gray", ms=10, capsize=2)
errorbar(value(my_res[4]), my_y[1:end-1], 0.0, xerr=err(my_res[4]), fmt="o", color="red", ms=10, capsize=2)
errorbar(value(my_res[4]), my_y[1:end-1], 0.0, xerr=sqrt(err(my_res[4])^2+my_syst[4]^2), fmt="o", color="red", ms=10, capsize=2)
# setp(ax2.get_yticklabels(),visible=false) # Disable y tick labels
PyPlot.yticks(y, fill("", length(y)), rotation=0)

vec = vcat(data_res,latt_res,my_res[4]+uwreal([0.0,my_syst[4]],"syst res")); uwerr.(vec)
trueres = value.(vec) .!= 0.0
xlim([minimum(value.(vec[trueres]).-err.(vec[trueres]))-0.1,maximum(value.(vec[trueres]).+err.(vec[trueres]))+0.1])

tight_layout()
display(fig)
close("all")

#-- 5.1.2: Final result vs data driven & other lattice

@info("Plot 4.1.2: Final result vs data driven & other lattice")

data_res = [
[uwreal([-20.623,0.130],"res a (Jegerlehner)"),uwreal([10.349,0.063],"res b (Jegerlehner)"),uwreal([0.337,0.005],"res c (Jegerlehner)"),uwreal([-9.927,0.067],"res (Jegerlehner)")],
[uwreal([-20.77,0.08],"res a (Keshavarzi et al.)"),uwreal([10.62,0.04],"res b (Keshavarzi et al.)"),uwreal([0.34,0.01],"res c (Keshavarzi et al.)"),uwreal([-9.82,0.04],"res c (Keshavarzi et al.)")]
]
my_res   = [uwreal([-20.92,0.42],"my res a"),uwreal([10.49,0.28],"my res b"),uwreal([0.338,0.013],"my res c"),uwreal([-10.07,0.14],"my res c")]; uwerr.(my_res)
my_syst  = [0.57,0.35,0.014,0.21]





#-- 5.1.3: Other result comparison

amuLO   = [uwreal([977.4642,88.4803],"LO amu"),uwreal([717.7441,18.609],"LO amu"),uwreal([719.9554,19.4075],"LO amu")]
amusyst = [uwreal([0.0,335.4771],"LO syst"),uwreal([0.0,40.0547],"LO syst"),uwreal([0.0,32.8617],"LO syst")]

res     = [uwreal([720.0,12.5],"LO res")]
ressyst = [uwreal([0.0,9.9],"LO res syst")]


fig = figure(figsize=(8,4))

label = ["All","All-(a4)","All-(a4,a2phi4)","1904.03120v1"]
for i=collect(1:length(amuLO))
    uwerr(amuLO[i])
    amuLOtot = amuLO[i]+amusyst[i]; uwerr(amuLOtot)
    errorbar(value(amuLO[i]), i, 0.0, xerr=err(amuLO[i]), fmt="o", color="black", ms=10, capsize=2)
    errorbar(value(amuLOtot), i, 0.0, xerr=err(amuLOtot), fmt="o", color="black", ms=10, capsize=2)
end

for i=collect(1:length(res))
    uwerr(res[i])
    restot = res[i]+ressyst[i]; uwerr(restot)
    errorbar(value(res[i]), length(amuLO)+i, 0.0, xerr=err(res[i]), fmt="x", color="gray", ms=10, capsize=2)
    errorbar(value(restot), length(amuLO)+i, 0.0, xerr=err(restot), fmt="x", color="gray", ms=10, capsize=2)
end


PyPlot.yticks(collect(1:length(amuLO)+length(res)), label, rotation=0)


tight_layout()
display(fig)
close("all")


#-- 5.1.2: Final result vs data driven & other lattice

@info("Plot 4.1.2: Final result vs data driven & other lattice")

data_res = [
[uwreal([-20.623,0.130],"res a (Jegerlehner)"),uwreal([10.349,0.063],"res b (Jegerlehner)"),uwreal([0.337,0.005],"res c (Jegerlehner)")],
[uwreal([-20.77,0.08],"res a (Keshavarzi et al.)"),uwreal([10.62,0.04],"res b (Keshavarzi et al.)"),uwreal([0.34,0.01],"res c (Keshavarzi et al.)")]
]
my_res   = [uwreal([-21.47,0.40],"my res a"),uwreal([10.95,0.28],"my res b"),uwreal([0.363,0.013],"my res c")]; uwerr.(my_res)
my_syst  = [0.93,0.59,0.0251]

data_tot = sum.(data_res); uwerr.(data_tot)
latt_tot = uwreal([],"")
my_tot = sum(my_res); uwerr.(my_tot)

#-- 5.1.3: Other result comparison

amuLO   = [uwreal([977.4642,88.4803],"LO amu"),uwreal([717.7441,18.609],"LO amu"),uwreal([719.9554,19.4075],"LO amu")]
amusyst = [uwreal([0.0,335.4771],"LO syst"),uwreal([0.0,40.0547],"LO syst"),uwreal([0.0,32.8617],"LO syst")]

res     = [uwreal([720.0,12.5],"LO res")]
ressyst = [uwreal([0.0,9.9],"LO res syst")]


fig = figure(figsize=(8,4))

# label = ["All","All-(a4)","All-(a4,a2phi4)","1904.03120v1"]
label = ["All","All-(a4)","All-(a4,a2phi4)","[1]"]
for i=collect(1:length(amuLO))
    uwerr(amuLO[i])
    amuLOtot = amuLO[i]+amusyst[i]; uwerr(amuLOtot)
    errorbar(value(amuLO[i]), i, 0.0, xerr=err(amuLO[i]), fmt="o", color="black", ms=10, capsize=2)
    errorbar(value(amuLOtot), i, 0.0, xerr=err(amuLOtot), fmt="o", color="black", ms=10, capsize=2)
end

for i=collect(1:length(res))
    uwerr(res[i])
    restot = res[i]+ressyst[i]; uwerr(restot)
    errorbar(value(res[i]), length(amuLO)+i, 0.0, xerr=err(res[i]), fmt="x", color="gray", ms=10, capsize=2)
    errorbar(value(restot), length(amuLO)+i, 0.0, xerr=err(restot), fmt="x", color="gray", ms=10, capsize=2)
end


PyPlot.yticks(collect(1:length(amuLO)+length(res)), label, rotation=0)


tight_layout()
display(fig)
close("all")


#-- 5.1.4: Window paper plots

amu     = [uwreal([187.29,0.663],"LO amu"),uwreal([717.7441,18.609],"LO amu"),uwreal([719.9554,19.4075],"LO amu")]
amusyst = [uwreal([0.0,1.146],"LO syst"),uwreal([0.0,40.0547],"LO syst"),uwreal([0.0,32.8617],"LO syst")]

res     = [uwreal([720.0,12.5],"LO res")]
ressyst = [uwreal([0.0,9.9],"LO res syst")]


fig = figure(figsize=(8,4))

# label = ["All","All-(a4)","All-(a4,a2phi4)","1904.03120v1"]
label = ["All","All-(a4)","All-(a4,a2phi4)","[1]"]
for i=collect(1:length(amuLO))
    uwerr(amuLO[i])
    amuLOtot = amuLO[i]+amusyst[i]; uwerr(amuLOtot)
    errorbar(value(amuLO[i]), i, 0.0, xerr=err(amuLO[i]), fmt="o", color="black", ms=10, capsize=2)
    errorbar(value(amuLOtot), i, 0.0, xerr=err(amuLOtot), fmt="o", color="black", ms=10, capsize=2)
end

for i=collect(1:length(res))
    uwerr(res[i])
    restot = res[i]+ressyst[i]; uwerr(restot)
    errorbar(value(res[i]), length(amuLO)+i, 0.0, xerr=err(res[i]), fmt="x", color="gray", ms=10, capsize=2)
    errorbar(value(restot), length(amuLO)+i, 0.0, xerr=err(restot), fmt="x", color="gray", ms=10, capsize=2)
end


PyPlot.yticks(collect(1:length(amuLO)+length(res)), label, rotation=0)


tight_layout()
display(fig)
close("all")

#---5 something, comparison plots


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

## EXTRA PLOTS

@info("Result comparison") # no need to read fits and MA  (just mykeys dict)

# IMPR_SET_COMB = [["1"],["1old"],["2"],["1old","2"],["1","2"]]
IMPR_SET_COMB = [["1"],["2"],["1","2"]]

STD_DERIV  = false
SimpleBase = false

sqrtt0_ph_TAR = nothing  # 0.1443  sqrtt0_ph_CLS  nothing

if comp in ["gCCconn","∆lc_b"]
    mDs_SHIFT = true
    MD_ph_prime = BDIOread_mDs(path_bdio)["mDs ph"]
else
    mDs_SHIFT = false
end

factor = Dict(
    "g33" => 1., "g88" => 1/3., "gCCconn" => 4/9., "∆ls_amu" => 1/3., "∆lc_b" => 4/9.,
    "g33s" => "", "g88s" => "(1/3)", "gCCconns" => "(4/9)", "∆ls_amus" => "(1/3)", "∆lc_bs" => "(4/9)"
)

mykeys = DictComptoKey[comp]

# paper_res = 35.03
# paper_err = [0.04,sqrt(0.04^2+0.21^2)]

# paper_res = -0.495
# paper_err = [0.007,sqrt(0.007^2+0.034^2)]

# paper_res = 186.30
# paper_err = [0.75,sqrt(0.75^2+1.08^2)]

# paper_res = −2.42
# paper_err = [0.04,sqrt(0.04^2+0.1^2)]

# paper_res = 6.81
# paper_err = [0.09,sqrt(0.09^2+0.21^2)]

# paper_res = 378.7
# paper_err = [3.7,sqrt(3.7^2+3.1^2)]

# paper_res = 44.5
# paper_err = [1.2,sqrt(1.2^2+1.1^2)]

paper_res = 0.01409
paper_err = [0.00035,sqrt(0.00035^2+0.00060^2)]

fig = figure(figsize=(6,6))
y = 0.0
for (i,IMPR_SET) in enumerate(IMPR_SET_COMB[length.(IMPR_SET_COMB) .== 1])
    impr_set = IMPR_SET[1]
    for (k,key) in enumerate(mykeys)
        y += 1.0

        myres = factor[comp] * MA[impr_set]["res"][key]; uwerr(myres)
        syst  = factor[comp] * info[impr_set]["syst"][key]

        if !isnothing(sqrtt0_ph_TAR)
            der = mchist(myres, "sqrtt0 [fm]")[1] / artificial_err
            res_t0shift = myres + value(sqrtt0_ph_TAR - sqrtt0_ph) * der; uwerr(res_t0shift)
        else
            res_t0shift = myres
        end

        errorbar(value(myres), xerr=err(myres), y-0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
        errorbar(value(myres), xerr=sqrt(err(myres)^2+syst^2), y-0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
        if !isnothing(sqrtt0_ph_TAR)
            errorbar(value(res_t0shift), xerr=err(res_t0shift), y+0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
            errorbar(value(res_t0shift), xerr=sqrt(err(res_t0shift)^2+syst^2), y+0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
        end
        if mDs_SHIFT
            der_mDs = mchist(myres, "MD_ph [GeV]")[1] / artificial_err
            res_t0mDsshift = res_t0shift + value(MD_ph - MD_ph_prime) * der_mDs; uwerr(res_t0mDsshift)
            errorbar(value(res_t0mDsshift), xerr=err(res_t0mDsshift), y+0.3, 0.0, fmt="^", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
            errorbar(value(res_t0mDsshift), xerr=sqrt(err(res_t0mDsshift)^2+syst^2), y+0.3, 0.0, fmt="^", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
        end
    end
    if length(mykeys) > 1
        y += 1.0

        weight = vcat([info[impr_set]["weight"][key] for key in mykeys]...)./(length(mykeys))
        res_tot = vcat([res[impr_set][key] for key in mykeys]...)

        amu, syst_ = model_average(res_tot, weight); uwerr.(amu)

        myres  = factor[comp] * amu[1]; uwerr(myres)
        syst = factor[comp] * syst_
        
        if !isnothing(sqrtt0_ph_TAR)
            der = mchist(myres, "sqrtt0 [fm]")[1] / artificial_err
            res_t0shift = myres + value(sqrtt0_ph_TAR - sqrtt0_ph) * der; uwerr(res_t0shift)
        else
            res_t0shift = myres
        end

        errorbar(value(myres), xerr=err(myres), y-0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
        errorbar(value(myres), xerr=sqrt(err(myres)^2+syst^2), y-0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
        if !isnothing(sqrtt0_ph_TAR)
            errorbar(value(res_t0shift), xerr=err(res_t0shift), y+0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
            errorbar(value(res_t0shift), xerr=sqrt(err(res_t0shift)^2+syst^2), y+0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
        end
        if mDs_SHIFT
            der_mDs = mchist(myres, "MD_ph [GeV]")[1] / artificial_err
            res_t0mDsshift = res_t0shift + value(MD_ph - MD_ph_prime) * der_mDs; uwerr(res_t0mDsshift)
            errorbar(value(res_t0mDsshift), xerr=err(res_t0mDsshift), y+0.3, 0.0, fmt="^", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
            errorbar(value(res_t0mDsshift), xerr=sqrt(err(res_t0mDsshift)^2+syst^2), y+0.3, 0.0, fmt="^", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
        end
    end
end
for (i,IMPR_SET) in enumerate(IMPR_SET_COMB[length.(IMPR_SET_COMB) .!= 1])    
    y += 1.0

    weight = vcat([vcat([info[impr_set]["weight"][key] for key in mykeys]...) for impr_set in IMPR_SET]...)./(length(mykeys)*length(IMPR_SET))
    res_tot = vcat([vcat([res[impr_set][key] for key in mykeys]...) for impr_set in IMPR_SET]...)

    amu, syst_ = model_average(res_tot, weight); uwerr.(amu)

    myres  = factor[comp] * amu[1]; uwerr(myres)
    syst = factor[comp] * syst_

    if !isnothing(sqrtt0_ph_TAR)
        der = mchist(myres, "sqrtt0 [fm]")[1] / artificial_err
        res_t0shift = myres + value(sqrtt0_ph_TAR - sqrtt0_ph) * der; uwerr(res_t0shift)
    else
        res_t0shift = myres
    end

    errorbar(value(myres), xerr=err(myres), y-0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
    errorbar(value(myres), xerr=sqrt(err(myres)^2+syst^2), y-0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
    if !isnothing(sqrtt0_ph_TAR)
        errorbar(value(res_t0shift), xerr=err(res_t0shift), y+0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
        errorbar(value(res_t0shift), xerr=sqrt(err(res_t0shift)^2+syst^2), y+0.1, 0.0, fmt="o", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
    end
    if mDs_SHIFT
        der_mDs = mchist(myres, "MD_ph [GeV]")[1] / artificial_err
        res_t0mDsshift = res_t0shift + value(MD_ph - MD_ph_prime) * der_mDs; uwerr(res_t0mDsshift)
        errorbar(value(res_t0mDsshift), xerr=err(res_t0mDsshift), y+0.3, 0.0, fmt="^", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
        errorbar(value(res_t0mDsshift), xerr=sqrt(err(res_t0mDsshift)^2+syst^2), y+0.3, 0.0, fmt="^", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
    end
end
gca().axvspan(paper_res-paper_err[1], paper_res+paper_err[1], color="limegreen", alpha=0.3)
gca().axvspan(paper_res-paper_err[2], paper_res+paper_err[2], color="limegreen", alpha=0.3, label="SD paper result")
diag_str = diag == "NLOa&b" ? "\\rm{NLOa}\\&\\rm{b}" : diag
comp_str = comp[1] == '∆' ? (comp[end] != 'b' ? "\\Delta_{ls}(a_{\\mu})" : "\\Delta_{lc}(b)") : "\\mathrm{$comp}"
wind_str = wind != "NW" ? "_{\\mathrm{$wind}}" : ""
Q_str    = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "($Q\\ GeV)" : ""
xlabel(latexstring("$(factor[comp*"s"])a_{\\mu}^{\\rm{hvp}}[\\rm{$(diag_str)}^{$comp_str}$wind_str]$Q_str"))
PyPlot.yticks(1:7,["VV - set 1","VVc - set 1","set 1","VV - set 2","VVc - set 2","set 2","sets [1,2]"],rotation = 0, fontsize=22)
# PyPlot.yticks(1:11,["VV - set 1","VVc - set 1","set 1","VV - set 1old","VVc - set 1old","set 1old","VV - set 2","VVc - set 2","set 2","sets [1old,2]","sets [1,2]"],rotation = 0, fontsize=22)
# PyPlot.yticks(1:5,["VVc - set 1old","VVc - set 1","VVc - set 2","sets [1old,2]","sets [1,2]"],rotation = 0, fontsize=22)
# PyPlot.yticks([])
# gca()[:tick_params](axis="x", pad=-130)  # Shift ticks closer to the plot
tight_layout()
display(gcf())
close()

#----

comp = "gCCconn"

# MYRES  = [uwreal([379.243,4.735],"a"),uwreal([382.271,4.649],"a"),uwreal([379.272,4.690],"a"),uwreal([381.820,4.651],"a"),uwreal([380.651,4.511],"a")]
# MYSYST = [3.942,4.116,4.63,4.074,4.427]
# MYRES  = [uwreal([133.681,3.785],"a"),uwreal([137.573,3.897],"a"),uwreal([134.456,3.678],"a"),uwreal([138.399,3.859],"a"),uwreal([136.027,3.779],"a")]
# MYSYST = [1.719,2.556,1.623,2.568,2.946]
MYRES  = [uwreal([0.01366,0.00037],"a"),uwreal([0.01481,0.00028],"a"),uwreal([0.01341,0.0003],"a"),uwreal([0.01502,0.00029],"a"),uwreal([0.01423,0.0003],"a")]
MYSYST = [0.00030,0.00047,7.1e-5,0.00044,0.00081]


factor = Dict(
    "g33" => 1., "g88" => 1/3., "gCCconn" => 4/9., "∆ls_amu" => 1/3., "∆lc_b" => 4/9.,
    "g33s" => "", "g88s" => "(1/3)", "gCCconns" => "(4/9)", "∆ls_amus" => "(1/3)", "∆lc_bs" => "(4/9)"
)

# paper_res = 378.7
# paper_err = [3.7,sqrt(3.7^2+3.1^2)]

# paper_res = 44.5
# paper_err = [1.2,sqrt(1.2^2+1.1^2)]

paper_res = 0.01409
paper_err = [0.00035,sqrt(0.00035^2+0.00060^2)]

fig = figure(figsize=(6,4))
y = 0.0
for i=1:length(MYRES)
    y += 1.0

    myres = MYRES[i]; uwerr(myres)
    syst  = MYSYST[i]

    errorbar(value(myres), xerr=err(myres), y, 0.0, fmt="d", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
    errorbar(value(myres), xerr=sqrt(err(myres)^2+syst^2), y, 0.0, fmt="d", mfc="none", color="gray", ms=10, capsize=2, alpha=0.4)
end
gca().axvspan(paper_res-paper_err[1], paper_res+paper_err[1], color="limegreen", alpha=0.3)
gca().axvspan(paper_res-paper_err[2], paper_res+paper_err[2], color="limegreen", alpha=0.3, label="SD paper result")
diag_str = diag == "NLOa&b" ? "\\rm{NLOa}\\&\\rm{b}" : diag
comp_str = comp[1] == '∆' ? (comp[end] != 'b' ? "\\Delta_{ls}(a_{\\mu})" : "\\Delta_{lc}(b)") : "\\mathrm{$comp}"
wind_str = wind != "NW" ? "_{\\mathrm{$wind}}" : ""
Q_str    = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "($Q\\ GeV)" : ""
xlabel(latexstring("$(factor[comp*"s"])a_{\\mu}^{\\rm{hvp}}[\\rm{$(diag_str)}^{$comp_str}$wind_str]$Q_str"))
PyPlot.yticks(1:5,["VV - set 1","VVc - set 1","VV - set 2","VVc - set 2","TOTAL"],rotation = 0, fontsize=22)
tight_layout()
display(gcf())
close()

