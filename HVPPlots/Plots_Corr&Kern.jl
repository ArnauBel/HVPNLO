# Import packages

using Revise

include("../HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs

using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2
using DelimitedFiles

using Plots
using PyPlot
using Colors

using ProgressBars
using Suppressor

# include uwreal constants

include("../HVPtool/uwConst.jl")

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_coef = joinpath(julia_script_directory, "..", "..", "KernelCoeff")

path_spec = joinpath(julia_script_directory, "..", "..", "HVPData", "spectroscopy")

pFVC_MLL  = joinpath(julia_script_directory, "..", "..", "HVPData", "FSE_MLL")


path_plot = joinpath(julia_script_directory,"..","..","Slides & Plots","Plots")

# We do not have charm or disconnected data for some of the ensembles

ensNOcharm = ["C102","D150","D201","D251","D451","F300","H200","H650","J304","J306","J307","J501","N451","N452"]
ensNOdisc  = ["F300","J306"]

ensSPECdata = ["D200","E250"]  # J303

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

# Set plot parameters

ens = "D150"; ens = EnsInfo(ens)

readIMPR_SET = ["1","2"] # ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

RESC       = false
STD_DERIV  = false

BLIND = true

path_bdio = path_bdio_dict["local"]

# Data reading and definitions

@info("Reading obsBDIO data:\n - ens: $(ens.id)")

println("- Reading t0...")

t0 = BDIOread_t0(path_bdio, ens)
t0su3 = t0sym(ens.beta)

println("- Reading fPi...")

fPi = BDIOread_fPS(path_bdio,ens)["fPi"]

println("- Reading TMR...")

diag = BLIND ? "NLO_1D" : "1D"

TMR   = BDIOread_TMR(path_bdio,ens.id,resc=RESC,diag,beta=false,BLIND=BLIND)
TMR_b = BDIOread_TMR(path_bdio,ens.id,resc=RESC,diag,beta=true,BLIND=false)

TMR["TMRa&b"] = TMR["TMRa"] + TMR["TMRb"]
TMR_b["TMRa&b"] = TMR_b["TMRa"] + TMR_b["TMRb"]

corr = Dict(); fvc_hp_dict = Dict()
# HVP  = Dict(); FVC = Dict()
for impr_set in readIMPR_SET
    println("- Reading corr (impr. set $impr_set)...")

    corr[impr_set] = BDIOread_corr(path_bdio,ens,impr_set,STD=STD_DERIV)

    if ens.kappa_l == ens.kappa_s || ens.id in ensNOdisc
        println("      - SU(3) flavour sym point or no disc. data available")
    else
        corr[impr_set]["g88_ll"] = corr[impr_set]["g88conn_ll"] .+ corr[impr_set]["g88disc_ll"] .+ (2).*corr[impr_set]["g08conn_ll"] .+ corr[impr_set]["g08disc_ll"] .+ corr[impr_set]["g80disc_ll"]
        corr[impr_set]["g88_lc"] = corr[impr_set]["g88conn_lc"] .+ corr[impr_set]["g88disc_lc"] .+ corr[impr_set]["g08conn_lc"] .+ corr[impr_set]["g08disc_lc"]
    end
end

println("- Reading fvc...")

fvc_hp_dict    = BDIOread_FVCcorr(path_bdio,ens)
fvc_hpRef_dict = BDIOread_FVCcorr(path_bdio,ens,Vref=true)

if ens.id ∉ ["A653","H650"]
    # GS_data_matrix, headers = readdlm("../gsFVC/$(ens.id)_gs_fvc.txt", '\t', header=true)

    # t_gs   = GS_data_matrix[:, 1]
    # fvc_gs = [uwreal([GS_data_matrix[:, 2][k],GS_data_matrix[:, 3][k]],ens.id*"_gs") for k=1:length(t_gs)]

    t_gs, fvc_gs    = TXTread_FVCcorr_GS(pFVC_MLL,ens)
    t_gs, fvc_gsRef = TXTread_FVCcorr_GS(pFVC_MLL,ens,Vref=true)
else
    println("   (No gs fvc available for ens = $(ens.id))")
end

if ens.id in ["D200","E250","J303"]
    println("- Reading spectral data...")
    E, Z, Z_impr = get_spectr_data(path_spec,ens)
end

@info("All data ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<----------------------------------------------------- KERNEL PLOTS ----------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("SD-Substracted Kernel plot")

diag = ""  # [LO]  [NLOa,NLOb]  [LO,NLOa,NLOc]

SAVE     = false
OVERSAVE = false

QLIST = Qlist  # Required for wind = SDsub

t = collect(range(0,0.9,1000))
TMRDict = Dict{String,Vector{Float64}}(
    "LO" => hbarc^2 .* Tildef2((massmu/hbarc) .* t, path_coef),
    "NLOa" => hbarc^2 .* Tildef4a((massmu/hbarc) .* t, path_coef),
    "NLOb" => hbarc^2 .* Tildef4b((massmu/hbarc) .* t, path_coef),
    "NLOa&b" => hbarc^2 .* (Tildef4a((massmu/hbarc) .* t, path_coef) .+ Tildef4b((massmu/hbarc) .* t, path_coef))
)

fig = figure(figsize=(8,6))

myTMR  = TMRDict[diag]
myTMRw = myTMR .* Window("SD")(t)

TMRb(Q::Float64) = ((16/(Q/hbarc)^2)^2 * π^2 * (massmu/hbarc)^2) .* C4[diag].((massmu/hbarc).*t) .* sin.((Q/hbarc/4) .* t).^4
TMRsub(Q::Float64) = myTMRw .- (Window("SD")(0) .* TMRb(Q))

PyPlot.plot(t, (hbarc/massmu)^3*myTMRw./(t.^3), color = "black", label=L"Q=\infty")
for Q in sort(QLIST,rev=true)
    PyPlot.plot(t, (hbarc/massmu)^3*TMRsub(Q)./(t.^3), label="Q=$Q")
end

title("SD-substracted Kernel")
xlabel("t [fm]")
diag_str = diag != "NLOa&b" ? diag : "NLOa\\&b"
ylabel(latexstring("\\frac{1}{\\left(m_\\mu t\\right)^3}\\tilde{K}^{(\\rm{$diag_str})}_{\\rm{sub}}(m_\\mu t;Q)"))
legend()
tight_layout()
display(gcf())
if SAVE
    p = create_path(path_plot,["Other","SDsubKernel_$diag.pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p)
end
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Kernel comparison plot")

DIAG = ["LO","NLOa","NLOb","NLOa&b"]  # [LO]  [NLOa,NLOb]  [LO,NLOa,NLOc]   ["LO","NLOa","NLOb","NLOa&b"]

NORMtoLO = false

SAVE     = false
OVERSAVE = false

t = collect(range(0.0,6.0,1000))
TMRDict = Dict{String,Vector{Float64}}(
    "LO"     => hbarc^2 .* Tildef2((massmu/hbarc) .* t, path_coef),
    "NLOa"   => hbarc^2 .* Tildef4a((massmu/hbarc) .* t, path_coef),
    "NLOb"   => hbarc^2 .* Tildef4b((massmu/hbarc) .* t, path_coef),
    "NLOa&b" => hbarc^2 .* (Tildef4a((massmu/hbarc) .* t, path_coef) .+ Tildef4b((massmu/hbarc) .* t, path_coef))
)

color = ["orange","lightblue","blue","darkblue"]
linestyle = [":","--","--","-"]

fig = figure(figsize=(8,6))
if !NORMtoLO
    PyPlot.plot(t, 0.0.*t, color="gray", linestyle="--")
end
for (d,diag) in enumerate(DIAG)
    diag_str = diag != "NLOa&b" ? diag : "NLOa\\&b"
    if !NORMtoLO
        PyPlot.plot(t, TMRDict[diag], label=diag_str , color=color[d], linestyle=linestyle[d])
    else
        PyPlot.plot(t, TMRDict[diag]./TMRDict["LO"], label=diag_str , color=color[d], linestyle=linestyle[d])
    end
end

# title("Kernels")
xlabel(latexstring("t\\ [\\textrm{fm}]"))
if !NORMtoLO
    ylabel(latexstring("\\tilde{f}^{(i)}(m_\\mu t)"))
else
    ylabel(latexstring("\\tilde{f}^{(i)}(m_\\mu t)/\\tilde{f}^{(\\rm{LO})}(m_\\mu t)"))
end
xlim(0.0,5.0)
ylim(-100,100)
legend()
tight_layout()
display(gcf())
if SAVE
    p = create_path(path_plot,["Other","Kernel_$(DIAG).pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p)
end
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Kernel ratio plot")

diag_ref = "LO"
diag     = "NLOa"

t = collect(range(0.0,10.0,1000))
TMRDict = Dict{String,Vector{Float64}}(
    "LO" => hbarc^2 .* Tildef2((massmu/hbarc) .* t, path_coef),
    "NLOa" => hbarc^2 .* Tildef4a((massmu/hbarc) .* t, path_coef),
    # "NLOb" => hbarc^2 .* Tildef4b((massmu/hbarc) .* t, path_coef),
    # "NLOa&b" => hbarc^2 .* (Tildef4a((massmu/hbarc) .* t, path_coef) .+ Tildef4b((massmu/hbarc) .* t, path_coef))
)

# color = ["orange","lightblue","darkblue","blue"]
# linestyle = ["-","-","-","--"]

fig = figure(figsize=(8,6))
PyPlot.plot(t, TMRDict[diag]./TMRDict[diag_ref]) # , color=color[d], linestyle=linestyle[d])


title("Kernels")
xlabel("t [fm]")
diag_str = diag != "NLOa&b" ? diag : "NLOa\\&b"
diag_ref_str = diag_ref != "NLOa&b" ? diag : "NLOa\\&b"
ylabel(latexstring("\\tilde{f}^{(\\rm{$diag_str})}(m_\\mu t)\\ /\\ \\tilde{f}^{(\\rm{$diag_ref_str})}(m_\\mu t)"))
legend()
tight_layout()
ylim(bottom = -18)

# Remove y-axis ticks and labels
# PyPlot.gca().axes.get_yaxis().set_visible(false)

display(gcf())
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Kernel error plot")

# rtol = 1e-8
tmax = 15.0 # fm

t = collect(range(0.0,tmax,10))

TMRDict = Dict{String,Vector{Float64}}(
    "LO" => hbarc^2 .* Tildef2((massmu/hbarc) .* t, path_coef),
    "NLOa" => hbarc^2 .* Tildef4a((massmu/hbarc) .* t, path_coef),
    "NLOb" => hbarc^2 .* Tildef4b((massmu/hbarc) .* t, path_coef, OverErr=true),
    "NLOa&b" => hbarc^2 .* (Tildef4a((massmu/hbarc) .* t, path_coef) .+ Tildef4b((massmu/hbarc) .* t, path_coef, OverErr=true))
)


TMRnumDict = Dict{String,Vector{Float64}}(
    "LO" => hbarc^2 .* Tildef2_num.((massmu/hbarc) .* t),
    "NLOa" => hbarc^2 .* Tildef4a_num.((massmu/hbarc) .* t),
    "NLOb" => hbarc^2 .* Tildef4b_num.((massmu/hbarc) .* t),
    "NLOa&b" => hbarc^2 .* (Tildef4a_num.((massmu/hbarc) .* t) .+ Tildef4b_num.((massmu/hbarc) .* t))
)

##

diag = "NLOa"

fig = figure(figsize=(8,6))
title("Error (expanded vs. numeric)")

PyPlot.plot(t, abs.(TMRDict[diag].-TMRnumDict[diag]))
yscale("log")
xlabel("t [fm]")
diag_str = diag != "NLOa&b" ? diag : "NLOa\\&b"
ylabel(latexstring("|\\tilde{f}^{(\\rm{$diag_str})}(m_\\mu t)\\ -\\ \\tilde{f}^{(\\rm{$diag_str})}(m_\\mu t)|"))
legend()
tight_layout()

display(gcf())
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<--------------------------------------------------- INTEGRAND PLOTS ---------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("1D integrand plot")

# comp = "g33_ll"
impr_set = "2"

DIAG = ["NLOb"]  # ["LO"]  ["NLOa","NLOb"]  ["LO","NLOa","NLOb"]
COMP = ["g33_ll"]  # ["g33_ll"]  ["g33_ll","g88_ll"]  ["g33_ll","g88_ll","gCCconn_ll"]
WIND = ["NW"]  # ["SD","SDsub"]  ["SD","ILD"]  ["SD","ID","LD"]  ["NW","SD","ILD"]  ["NW","SD","ID","LD"]

SAVE     = false
OVERSAVE = false

QLIST = Qlist  # Required for wind = SDsub

T = HVPobs.Data.get_T(ens.id)
t = collect(1:Int64(T/2+1))

fig = figure(figsize=(10,6))

sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1); t = collect(1:sym_points)

# label = ["Isovector (3,3)","Isoscalar (8,8)","Charm connected (c,c)"]
label  = ["No window","W='SD'","W='ID'","W='LD'"]
colour = ["black","lightgreen","brown","purple"]

i = 0
for diag in DIAG
    for wind in sort(WIND, by = x -> Dict(s => i for (i, s) in enumerate(["NW","SD","SDsub","ID","LD","ILD"]))[x])
        for comp in COMP
            i += 1
            sym = comp in ["gCCconn_ll_sim","gCCconn_lc_sim","gCCconn_ll_sim+","gCCconn_lc_sim"]

            factor = !sym ? hbarc * sqrt(t0)/sqrtt0_ph : hbarc * sqrt(t0su3)/sqrtt0_ph
            aens = !sym ? (sqrtt0_ph / sqrt(t0)) : (sqrtt0_ph / sqrt(t0sym(ens.beta)))
            tfm = aens.*(t.-1)

            for Q in (wind == "SDsub" ? sort(QLIST,rev=true) : 1.0)
                tmr = diag == "LO" ? TMR["TMR"] : (diag != "NLOa&b" ? TMR["TMR"*diag[end]] : TMR["TMRa"] .+ TMR["TMRb"])
                tmr = tmr[t]
                if sym
                    tmr_b = diag == "LO" ? TMR_b["TMR"] : (diag != "NLOa&b" ? TMR_b["TMR"*diag[end]] : TMR_b["TMRa"] .+ TMR_b["TMRb"])
                    tmr_b = tmr_b[t]
                end
                if wind == "NW"
                    TMRw = !sym ? tmr : tmr_b
                elseif wind == "SDsub"
                    TMRw_ = !sym ? (tmr .* Window("SD")(tfm)) : (tmr .* Window("SD")(tfm_SU3))
                    TMRb(Q::Float64) = ((16/(Q/factor)^2)^2 * π^2 * (massmu/factor)^2) .* C4[diag].((massmu/factor).*(t .- 1)) .* sin.((Q/factor/4) .* (t .- 1)).^4
                    TMRsub(Q::Float64) = TMRw_ .- (Window("SD")(0) .* TMRb(Q))
                    TMRw = TMRsub(Q)
                else
                    TMRw = !sym ? (tmr .* Window(wind)(tfm)) : (tmr .* Window(wind)(tfm_SU3))
                end

                int = ((alpha/π)^2*1e10) .* corr[impr_set][comp][t] .* TMRw; uwerr.(int)

                diag_str = diag != "NLOa&b" ? diag : "NLOa\\&b"
                comp_str = comp[2]*","*comp[3]
                wind_str = wind != "NW" ? "; $wind" : ""
                Qstr = wind == "SDsub" ? "; Q=$Q GeV" : ""
                errorbar(value.(tfm), value.(int), err.(int),  c=colour[i], capsize=2, fmt="-o", mfc="none", label=label[i])  #  label[i])  "$diag_str; $comp_str$wind_str$Qstr"
            end
        end
    end
end
# title("Integrands for $(ens.id) impr. set $(impr_set)")
title(ens.id)
xlabel(latexstring("t\\ [\\rm{fm}]"))
ylabel(latexstring("\\tilde{f}^{(4b)}(\\hat{t})\\times G^{(3,3)}(t)\\times\\Theta_{\\rm{W}}(t)"))
xMin, xMax = xlim()
xMin, xMax = any(x -> x in ["LD","ILD","NW"],WIND) ? [xMin,xMax] : ("ID" in WIND ? [-0.1,2.0] : [-0.05,1.0])
# xlim(xMin,xMax)
xlim(0.0,6.0)
ylim(bottom=-40.)
legend()
tight_layout()
display(gcf())
if SAVE
    RESCstr = !RESC ? "" : "_resc"
    # p = create_path(path_plot,["Other","IntegradIsospin.pdf"],OVERWRITE=OVERSAVE)
    p = create_path(path_plot,["Other","IntegradWindows_$(ens.id).pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p)
end
close()


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Integrand reconstruction plot")

discr    = "ll"
impr_set = "2"

diag = "LO" # [LO]  [NLOa,NLOb]  [LO,NLOa,NLOc]
wind = "NW" # NW  LD  ILD

wind ∉ ["NW","LD","ILD"] ? error("Window choice not possible for corr reconstruction\n - Please choose between 'NW', 'LD' or 'ILD'") : nothing
ens.id ∉ ["D200","E250","J303"] ? error("No spectroscopy data for ens $(ens.id)\n - Please choose between 'D200', 'E250' or 'J303'") : nothing

NMAX = 4 # nothing  E250: 4 ,  D200: 2


comp = "g33_"*discr

fig = figure(figsize=(12,6))

sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1); t = collect(1:sym_points)
aens = (sqrtt0_ph / sqrt(t0))
tfm = aens.*(t.-1)


corrVec, corrVec_PiPi, corrVec_Impr = reconstr_corr(ens,E,Z,Z_impr,nmax=NMAX,impr_set=impr_set,IMPR=true,RENORM=true,total=true); [uwerr.(corr) for corr in corrVec]

tmr  = diag == "LO" ? TMR["TMR"] : (diag != "NLOa&b" ? TMR["TMR"*diag[end]] : TMR["TMRa"] .+ TMR["TMRb"])
TMRw = wind == "NW" ? tmr[t] : tmr[t] .* Window(wind)(tfm)

int  = corr[impr_set][comp][t] .* TMRw; uwerr.(int)
intR = [corr[t] .* TMRw for corr in corrVec]; [uwerr.(intr) for intr in intR]

errorbar(value.(tfm), value.(int), err.(int), capsize=2, fmt="o", color="black", label="Full inetgrand")
for (n,intr) in enumerate(intR)
    errorbar(value.(tfm), value.(intr), err.(intr), capsize=2, fmt="o", color="black", mfc="none", alpha=0.1*n, label="Reconstr. n = $n")
end

title("Integrands for $(ens.id) [impr. set $(impr_set), discr. $(comp[end-1:end])]")
xlabel("t [fm]")
diag_str = diag != "NLOa&b" ? diag : "NLOa\\&b"
wind_str = wind != "NW" ? ";$wind" : ""
ylabel(latexstring("K^{(\\rm{$diag_str$wind_str})}(m_\\mu t)\\times G(t)"))
xlim(0,4.3)
legend()
tight_layout()
display(gcf())
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------ FVC PLOTS ------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Gunaris-Sakurai vs Hansen-Patella FSE")

diag = "NLOa&b"

VREF = true

# Tover2 = Int64(HVPobs.Data.get_T(ens.id)/2+1)
# t = collect(1:Tover2)

T = HVPobs.Data.get_T(ens.id)+1
t = collect(1:T)

# tmr = uwreal.([0.0, 0.0012675632782413627, 0.02023464304487228, 0.10210555939582061, 0.3214066384157293]) 
tmr = (diag == "LO") ? TMR["TMR"] : TMR["TMR"*diag[4:end]]; tmr = tmr[t]

fvcPi_hp   = VREF ? vcat(tmr[1],-fvc_hpRef_dict["FVCPi6"][1:T-1]) : vcat(tmr[1],-fvc_hp_dict["FVCPi6"][1:T-1])
fvcPi_hpgs = VREF ? vcat(vcat(tmr[1],-fvc_hpRef_dict["FVCPi6"][1:Int64(t_gs[1])-1]),-fvc_gsRef...)[1:T] : vcat(vcat(tmr[1],-fvc_hp_dict["FVCPi6"][1:Int64(t_gs[1])-1]),-fvc_gs...)[1:T]
fvcK_hp = VREF ? vcat(tmr[1],-fvc_hpRef_dict["FVCK6"][1:T-1]) : vcat(tmr[1],-fvc_hp_dict["FVCK6"][1:T-1])


if ens.kappa_l != ens.kappa_s
    int_hp   = tmr .* (fvcPi_hp   .+ fvcK_hp)
    int_hpgs = tmr .* (fvcPi_hpgs .+ fvcK_hp)
else
    int_hp   = 3/2 .* (tmr .* fvcPi_hp)  
    int_hpgs = 3/2 .* (tmr .* fvcPi_hpgs)
end

aens = value(sqrtt0_ph / sqrt(t0))

if ens.id in ["H105","N200","N300","N302"] && !VREF
    enstoinf = Dict("H105"=>"N101","N200"=>"D251","N300"=>"J307","N302"=>"J306")

    corr_large = BDIOread_corr(path_bdio,enstoinf[ens.id],readIMPR_SET[1],STD=STD_DERIV)
    fvc_large_hp  = BDIOread_FVCcorr(path_bdio,enstoinf[ens.id])

    fvc_hp_dict_large = BDIOread_FVCcorr(path_bdio,enstoinf[ens.id])
    fvc_hpRef_dict_large = BDIOread_FVCcorr(path_bdio,enstoinf[ens.id],Vref=true)

    t_gs_large, fvc_large_mll = TXTread_FVCcorr_GS(pFVC_MLL,enstoinf[ens.id])
    fvcPi_hp_large   =  vcat(tmr[1],-fvc_hp_dict_large["FVCPi6"][1:Tover2-1])
    fvcPi_hpgs_large = vcat(vcat(tmr[1],-fvc_hp_dict_large["FVCPi6"][1:Int64(t_gs_large[1])-1]),-fvc_large_mll...)[1:Tover2]
    fvcK_hp_large  = vcat(tmr[1],-fvc_hp_dict_large["FVCK6"][1:Tover2-1])

    if ens.kappa_l != ens.kappa_s
        int_hp_large   = tmr .* (fvcPi_hp_large   .+ fvcK_hp_large)
        int_hpgs_large = tmr .* (fvcPi_hpgs_large .+ fvcK_hp_large)
    else
        int_hp_large   = 3/2 .* (tmr .* fvcPi_hp_large)  
        int_hpgs_large = 3/2 .* (tmr .* fvcPi_hpgs_large)
    end

    int_hp .-= int_hp_large
    int_hpgs .-= int_hpgs_large

    int_data = tmr .* (corr_large["g33_ll"][1:length(corr[readIMPR_SET[1]]["g33_ll"])] .- corr[readIMPR_SET[1]]["g33_ll"])[1:T]; uwerr.(int_data)
    title("$(enstoinf[ens.id]) vs. $(ens.id)")
else
    VREF ? title("FSE for $(ens.id) (to Lref)") : title("FSE for $(ens.id)")
end

uwerr.(int_hp)
uwerr.(int_hpgs)

errorbar(aens .* (collect(0:T-1).+0.3), value.(int_hp)  , err.(int_hp)  , fmt="^", mfc="none", color="blue", capsize=2, label="FSE, HP")
errorbar(aens .* (collect(0:T-1).+0.6), value.(int_hpgs), err.(int_hpgs), fmt="v", mfc="none", color="red" , capsize=2, label=L"FSE, HP$\&$MLL")
if ens.id in ["H105","N200","N300","N302"] && !VREF
    errorbar(aens .* collect(0:T-1), value.(int_data)  , err.(int_data)  , fmt="d", color="black", capsize=2, label="FSE, data")
end
axvline(aens * (t_gs[1]-1), color="gray", ls="--", alpha=0.6)
# xlabel(latexstring("t/a"))
xlabel(latexstring("t\\ [\\rm{fm}]"))
# ylabel(latexstring("\\Delta G(t)^{(3,3)}\\tilde{K}^{(\\rm{$diag})}(m_\\mu t)"))
# yscale("log")
# xlim(right=50)
xlim(right=3.3) # fm
legend()
tight_layout()
display(gcf())
close()


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------ BOUNDING METHOD PLOTS ------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Corr boundings & reconstruction plot")

discr    = "ll"
impr_set = "1"

SPEC_REC = true

NMAX = 4 # nothing  E250: 4 ,  D200: 2


tcut = 10
tEeffFix = 12


mpi  = m_ens[ens.id]["mPi"] 
mrho = m_ens[ens.id]["mRho"]; uwerr(mrho)
L = ens.L
E2pi = 2*sqrt(mpi^2 + (2*pi/L)^2); uwerr(E2pi)

sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1); t = collect(1:sym_points)
aens = sqrtt0_ph / sqrt(t0)
tfm = aens.*(t.-1)

E0_33 = mrho < E2pi ? mrho : E2pi
E0_88 = mrho

if ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s
    corr[impr_set]["g88_ll"] = corr[impr_set]["g88conn_ll"] .+ corr[impr_set]["g88disc_ll"] .+ (2).*corr[impr_set]["g08conn_ll"] .+ corr[impr_set]["g08disc_ll"] .+ corr[impr_set]["g80disc_ll"]
    corr[impr_set]["g88_lc"] = corr[impr_set]["g88conn_lc"] .+ corr[impr_set]["g88disc_lc"] .+ corr[impr_set]["g08conn_lc"] .+ corr[impr_set]["g08disc_lc"]
end

fig = figure(figsize=(8,6))

corr33 = corr[impr_set]["g33_$discr"][t]; uwerr.(corr33)
errorbar(value.(tfm), value.(corr33), err.(corr33), color="blue", fmt="-o", mfc="none", label="33")
if ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc
    corr88 = corr[impr_set]["g88_$discr"][t]; uwerr.(corr88)
    errorbar(value.(tfm), value.(corr88), err.(corr88), color="red",  fmt="-s", mfc="none", label="88")
end

UB33 = corr_bound(t, tcut, corr33, E0_33); uwerr.(UB33)
LB33 = corr_bound(t, tcut, corr33, Eeff(tEeffFix, corr33)); uwerr.(LB33)
if ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc
    UB88 = corr_bound(t, tcut, corr88, E0_88); uwerr.(UB88)
    LB88 = corr_bound(t, tcut, corr88, Eeff(tEeffFix, corr88)); uwerr.(LB88)
end

axvline((tEeffFix-1)*value(aens), color="black", alpha=0.5, linestyle="--", linewidth=2)

errorbar(value.(tfm[tcut+1:end]), value.(UB33), err.(UB33), color="lightblue", fmt="-o", mfc="none", capsize=2, label="33 UB")
errorbar(value.(tfm[tcut+1:end]), value.(LB33), err.(LB33), color="darkblue",  fmt="-o", mfc="none", capsize=2, label="33 LB")

if ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc
    errorbar(value.(tfm[tcut+1:end]), value.(UB88), err.(UB88), color="indianred", fmt="-s", mfc="none", capsize=2, label="88 UB")
    errorbar(value.(tfm[tcut+1:end]), value.(LB88), err.(LB88), color="darkred",   fmt="-s", mfc="none", capsize=2, label="88 LB")
end

ncol = 1
if ens.id in ["D200","E250","J303"] && SPEC_REC
    corrVec, corrVec_PiPi, corrVec_Impr = reconstr_corr(ens,E,Z,Z_impr,nmax=nothing,impr_set=impr_set,IMPR=true,RENORM=true,total=true)
    # corrVec = [corr[t] for corr in corrVec]
    if discr == "lc"
        corrVec = [corrRec .* (corr[impr_set]["g$(comp)_lc"][t]./corr[impr_set]["g$(comp)_ll"][t]) for corrRec in corrVec]
    end
    [uwerr.(corr) for corr in corrVec]

    for n=1:min(length(E),length(Z))
        errorbar(value.(tfm), value.(corrVec[n]), err.(corrVec[n]), color="gray", alpha=0.1*n, fmt="d", mfc="none", capsize=2, label="33 recons. (n=$n)")
    end
    ncol += 1
end

title("Correlator boundings for $(ens.id) [impr. set $impr_set, discr. $discr]")
yscale("log")
xlabel("t [fm]")
ylabel("G(t)")
legend(ncol=ncol)
tight_layout()
display(gcf())
close()


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Corr. saturation plot")

impr_set = "2"

NMAX = nothing # nothing  E250: 4 ,  D200: 2


comp = "g33_ll"

fig = figure(figsize=(8,6))

sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1); t = collect(1:sym_points)
aens = (sqrtt0_ph / sqrt(t0))
tfm = aens.*(t.-1)

corrVec, corrVec_PiPi, corrVec_Impr = reconstr_corr(ens,E,Z,Z_impr,nmax=NMAX,impr_set=impr_set,IMPR=true,RENORM=true,total=true); [uwerr.(corr) for corr in corrVec]
mycorr = corr[impr_set][comp][t]; uwerr.(mycorr)

fill_between(value.(tfm), 1.0 .+ err.(mycorr)./value.(mycorr), 1.0 .- err.(mycorr)./value.(mycorr), color="gray", alpha=0.3, label="LMA relative error")
for (n,corr) in enumerate(corrVec)
    errorbar(value.(tfm), value.(corr)./value.(mycorr), err.(corr)./abs.(value.(mycorr)), capsize=2, fmt="o", color="black", mfc="none", alpha=0.1*n, label="Reconstr. n = $n")
end
axvline(value(tfm[findfirst_uninterrupted(err.(corrVec[end]) .< err.(mycorr))]-aens/2), color="black", alpha=0.2, linestyle="dotted", linewidth=2)

title("Saturation for $(ens.id) [impr. set $(impr_set), discr. $(comp[end-1:end])]")
xlabel("t [fm]")
ylabel(latexstring("G_{\\rm{rec}}(t)/G_{\\rm{LMA}}(t)"))
xlim(-0.1,3.4)
ylim(-0.1,1.3)
legend()
tight_layout()
display(gcf())
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Bounding Method plot")

diag  = "NLOa&b"  # LO  NLOa  NLOb  NLOa&b
wind  = "LD"  # NW  LD  ILD
discr = "ll"  # ll  lc

impr_set = "1"

NMAX = 4  # nothing  D200: 2 ,  E250: 4

BMRec  = false
# BMImpr = false


tcut0 = 10
tstep = 1


mpi  = m_ens[ens.id]["mPi"]
mrho = m_ens[ens.id]["mRho"]
if ens.id == "H650"    
    uwerr(mrho); mrho = mrho.mean-mrho.err
    # mrho = 0.3646
end
E2pi = 2*sqrt(mpi^2 + (2π/ens.L)^2)
# uwerr(E0_ens[ens.id]["E0"]); E0 = E0_ens[ens.id]["E0"].mean + E0_ens[ens.id]["E0"].err

tEeffFix = 1.5 # round(4.5/E0)  1.2  1.5

tmr  = diag == "LO" ? TMR["TMR"] : (diag != "NLOa&b" ? TMR["TMR"*diag[end]] : TMR["TMRa"] .+ TMR["TMRb"])

sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1)
t   = collect(1:sym_points)
tBM = collect(1:length(tmr))

aens = !RESC ? sqrtt0_ph / sqrt(t0) : hbarc * fPi / fPi_ph
# aens = asym(ens.beta)
tfm = aens.*(tBM.-1)

exp_diag = diag == "LO" ? 2 : 3

tmrw = wind == "NW" ? tmr : tmr .* Window(wind)(tfm)

res =  Dict()
HVP = Dict(); HVPsyst = Dict(); plateau = Dict()
ub = Dict(); lb = Dict(); lb0 = Dict();  averb = Dict()
HVP33rec = 0.0
trec = 0.0
for comp in ((ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc) ? ["33","88conn","88"] : ["33"])
    
    ub[comp]  = Vector{uwreal}()
    lb[comp]  = Vector{uwreal}()
    lb0[comp] = Vector{uwreal}()

    E0 = 0.0
    if comp == "33"
        E0 = mrho < E2pi ? mrho : E2pi
    elseif comp in ["88conn","88"]
        E0 = mrho
    end

    # if ens.id in ["D200","J303"] && comp == "33" && BMImpr
    #     corrVec = reconstr_corr(ens,E,Z,Z_impr,nmax=2,impr_set=impr_set,IMPR=true,RENORM=true,total=false)
    #     # corrVec = [value.(corr[t]) for corr in corrVec]
    #     corrVec = [corr[t] for corr in corrVec]
    #     obs = corr[impr_set]["g$(comp)_$(discr)"][t] .- corrVec[end]
    #     intRec = corrVec[end] .* tmrw
    #     int = obs .* tmrw
    #     uwerr(E[length(corrVec)+1]); E0 = value(E[length(corrVec)+1])-err(E[length(corrVec)+1])
    #     Impr = true
    # else
    #     obs = corr[impr_set]["g$(comp)_$(discr)"]
    #     int = obs[t] .* tmrw
    #     E0 = min(value(E2pi)-err(E2pi),value(mrho)-err(mrho))
    #     Impr = false
    # end

    obs = corr[impr_set]["g$(comp)_$(discr)"]
    int = obs[t] .* tmrw[t]
    # E0 = min(value(E2pi)-err(E2pi),value(mrho)-err(mrho))

    Eeff_ = uwreal(0.0)

    for tcut in tcut0:tstep:t[end-1]  # compute the upper and lower corr bounds 

        if aens.mean.*tcut < tEeffFix  # we fix the eff energy at some point
            Eeff_=Eeff(tcut, obs)
        end
        UB = corr_bound(tBM, tcut, obs, E0)
        LB = corr_bound(tBM, tcut, obs, Eeff_)

        UBInt = UB .* tmrw[tcut+1:end]
        LBint = LB .* tmrw[tcut+1:end]

        ub_  = (alpha/pi)^exp_diag * (sum(int[1:tcut])+sum(UBInt)) * 1e10
        lb_  = (alpha/pi)^exp_diag * (sum(int[1:tcut])+sum(LBint)) * 1e10
        lb0_ = (alpha/pi)^exp_diag * (sum(int[1:tcut])) * 1e10

        push!(ub[comp], ub_)
        push!(lb[comp], lb_)
        push!(lb0[comp], lb0_)
    end
    # HVP[comp], HVPsyst[comp], plateau[comp], averb[comp] = bounding_method(ub[comp],lb[comp],aens,AVER=true,PLAT=true,tcut0=tcut0)
    HVP[comp], HVPsyst[comp], plateau[comp], averb[comp] = bounding_method(ub[comp],lb[comp],aens,AVER=true,PLAT=true,tcut0=tcut0,correlations=true)

    # if Impr
    #     hvpRec = (alpha/pi)^exp_diag * sum(intRec) * 1e10
    #     HVP[comp] += hvpRec
    #     ub[comp] .+= hvpRec; lb[comp] .+= hvpRec; lb0[comp] .+= hvpRec;averb[comp] .+= hvpRec
    # end
    uwerr(HVP[comp])
    uwerr.(ub[comp]); uwerr.(lb[comp]); uwerr.(lb0[comp]) # ; uwerr.(averb[comp])

    if ens.id in ensSPECdata && comp == "33" && BMRec
        corrVec = reconstr_corr(ens,E,Z,Z_impr,impr_set=impr_set,IMPR=true,RENORM=true,total=false)
        # corrVec = [corr[t] for corr in corrVec]
        uwerr.(corrVec[end]); uwerr.(obs[t])
        trec = findfirst_uninterrupted(err.(corrVec[end]) .< err.(obs[t]))
        if discr == "lc"
            corrRec = corrVec[end] .* (corr[impr_set]["g$(comp)_lc"][t]./corr[impr_set]["g$(comp)_ll"][t])
        elseif discr == "ll"
            corrRec = corrVec[end]
        end
        recInt = corrRec[trec:end] .* tmrw[trec:end]

        HVP33rec = (alpha/pi)^exp_diag * (sum(int[1:trec-1])+sum(recInt)) * 1e10
    end
end

# Plot:

tcut0_fm = value(aens).*(collect(tcut0:tstep:t[end-1]).-1)

for comp in ((ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc) ? ["33","88conn","88"] : ["33"])
    fig = figure(figsize=(8,6))

    # axvline(aens.mean*(tEeffFix-1), color="black", alpha=0.2, linestyle="--", linewidth=2)
    axvline(tEeffFix, color="black", alpha=0.2, linestyle="--", linewidth=2)

    errorbar(tcut0_fm, value.(ub[comp]), err.(ub[comp]), color="green", fmt="o", capsize=2, label="UB from E0")
    errorbar(tcut0_fm, value.(lb[comp]), err.(lb[comp]), color="darkblue",  fmt="s", capsize=2, label="LB from Eeff")
    # errorbar(tcut0_fm, value.(averb[comp]), err.(averb[comp]), color="orange",  fmt="^", capsize=2, mfc="none", label="Average")
    # yMin, yMax = ylim()
    errorbar(tcut0_fm, value.(lb0[comp]), err.(lb0[comp]), color="lightblue",  fmt="d", capsize=2, label="LB from 0")

    # [tcut0_fm[1],tcut0_fm[end]]
    # plateau[comp]
    fill_between(plateau[comp], value(HVP[comp])+err(HVP[comp]), value(HVP[comp])-err(HVP[comp]), color="orange", alpha=0.2, label="Bounding Method")
    fill_between(plateau[comp], value(HVP[comp])+sqrt(err(HVP[comp])^2+HVPsyst[comp]^2), value(HVP[comp])-sqrt(err(HVP[comp])^2+HVPsyst[comp]^2), color="orange", alpha=0.2)
    if ens.id in ensSPECdata && comp == "33"
        uwerr(HVP33rec)
        axvline(value(tfm[trec]-aens/2), color="black", alpha=0.2, linestyle="dotted", linewidth=2)
        fill_between([tcut0_fm[1], tcut0_fm[end]], value(HVP33rec)+err(HVP33rec), value(HVP33rec)-err(HVP33rec), color="red", alpha=0.2, label="Reconstructed")
    end

    title("Bounding Method for $(ens.id) [impr. set $impr_set, discr. $discr]  ($(print_uwreal(HVP[comp])))")
    # ylim(yMin,yMax)
    ylim(value(HVP[comp])-4*err(HVP[comp]),value(HVP[comp])+6*err(HVP[comp]))
    comp == "33" ? xlim(0.7,min(tcut0_fm[end],5)) : xlim(1,min(tcut0_fm[end],3))
    xlabel(latexstring("t_{\\rm{cut}}\\ [\\rm{fm}]"))
    diag_str = diag == "LO" ? "\\rm{LO}" : (diag == "NLOa&b" ? "\\rm{NLO}_{\\rm{a}\\&\\rm{b}}" : "\\rm{NLO}_{\\rm{$(diag[end])}}")
    comp_str = comp != "88conn" ? "\\mathrm{$(comp[1]),$(comp[2])}" : "\\mathrm{$(comp[1]),$(comp[2])(conn.)}"
    fact_str = "\\times10^{10}" # diag == "LO" ? "\\times10^{10}" : "\\times10^{11}"
    if wind == "NW"
        ylabel(latexstring("a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}](t_c)$fact_str"))
    else
        ylabel(latexstring("\\left(a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]\\right)^{\\rm{$wind}}(t_c)$fact_str"))
    end

    legend()
    tight_layout()
    display(gcf())
    close()
end


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Tail reconstruction: Small periodic boxes")

ensid = "A653"  # A653  A654  B450

diag     = "NLOa&b"  # LO  NLOa  NLOb  NLOa&b
wind     = "LD"  # NW  LD  ILD
key      = "g33_ll"  # g33_ll  g33_lc  g88_ll  g88_lc
impr_set = "2"

path_bdio = path_bdio_dict["local"]


p0_dict = Dict(
    "A653_set1_g33_ll" => 15:20,
    "A653_set1_g33_lc" => 15:20,
    # "A653_set2_g33_ll" => 22:22,
    "A653_set2_g33_ll" => 14:22,
    "A653_set2_g33_lc" => 15:16,

    "A654_set1_g33_ll" => 16:20,
    "A654_set1_g33_lc" => 15:17,
    "A654_set2_g33_ll" => 13:15,
    "A654_set2_g33_lc" => 14:16,

    "B450_set1_g33_ll" => 20:24,
    "B450_set1_g33_lc" => 20:24,
    "B450_set2_g33_ll" => 20:24,
    "B450_set2_g33_lc" => 20:24,

    "N452_set1_g33_ll" => 20:21,
)

ens = EnsInfo(ensid)

corr = BDIOread_corr(path_bdio,ens,impr_set,STD=false)
tmr  = BDIOread_TMR(path_bdio,ens,diag,beta=false,BLIND=BLIND)

t0 = BDIOread_t0(path_bdio,ens)

T = HVPobs.Data.get_T(ens.id)
sym_points = Int64(T/2+1)
t = collect(1:length(tmr))
aens = sqrtt0_ph / sqrt(t0); tfm = aens.*(t.-1)

# if ens.kappa_l != ens.kappa_s
#     corr["g88_ll"] = corr["g88conn_ll"] .+ corr["g88disc_ll"] .+ (2).*corr["g08conn_ll"] .+ corr["g08disc_ll"] .+ corr["g80disc_ll"]
#     corr["g88_lc"] = corr["g88conn_lc"] .+ corr["g88disc_lc"] .+ corr["g08conn_lc"] .+ corr["g08disc_lc"]
# end

tmrw  = (wind == "NW") ? (tmr) : (tmr .* Window(wind)(tfm))

obs = corr[key][1:sym_points]

@. exp_model(x0,p)   = p[2] * exp(-p[1] * x0)
@. cosh_model(x0,p)  = p[2] * (exp(-p[1] * x0) + exp(-p[1] * (T-x0))) # + p[4] * (exp(-p[3] * x0) + exp(-p[3] * (T-x0))) 

p0_tuple = p0_dict["$(ens.id)_set$(impr_set)_$(key)"]
fit_vec = []
for p0 in p0_tuple
    data = obs[p0:end] 
    fit  = fit_routine(cosh_model,collect(p0:sym_points).-1, data, 2, pval=true, info=false, lineprint=false)
    push!(fit_vec,fit)
end

w = get_w_from_fitres(vcat(fit_vec...), AIC=true)

param = getfield.(vcat(fit_vec...),:param)
p1_vec = [par[1] for par in param]
p2_vec = [par[2] for par in param]

p1_res, p1_sys = model_average(p1_vec, w)
p2_res, p2_sys = model_average(p2_vec, w)

p1 = p1_res[1] + uwreal([0.0,p1_sys],"p1 fit $(ens.id)"); uwerr(p1)
p2 = p2_res[1] + uwreal([0.0,p2_sys],"p2 fit $(ens.id)"); uwerr(p2)

obs_rec = vcat(obs[1:(p0_tuple[1]-1)+(argmax(w)-1)],exp_model(t[p0_tuple[1]+(argmax(w)-1):end].-1,[p1,p2])...)

# MA plot

len = length(obs)

fig = figure(figsize=(16,12))
gs = fig.add_gridspec(5, 1, height_ratios=[4, 1, 1, 1, 1])  # Adjust the height_ratios as needed

ax1 = fig.add_subplot(gs[1, 1])
title("$(ens.id)  [$(key[2:3])-$(key[end-1:end]) impr $impr_set] corr")

x0 = collect(Int64,p0_tuple[1]-5:len) .- 1
vec = obs[x0[1]+1:end]; uwerr.(vec)

errorbar(x0, value.(vec), err.(vec), fmt="o", capsize=2, color="black")

x_fit = collect(p0_tuple[1]-0.4:0.1:len).-1.
y_fit = cosh_model(x_fit,[p1,p2]); uwerr.(y_fit)


fill_between(x_fit, value.(y_fit)+err.(y_fit), value.(y_fit)-err.(y_fit), alpha=0.3, color="orange")

axvline(x=p0_tuple[1]-1.0-0.4, color="red", linestyle="--")
axvline(x=p0_tuple[end]-1.0+0.4, color="red", linestyle="--")

yscale("log")
axis("tight")
setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels

ax21 = fig.add_subplot(gs[2, 1])

uwerr.(p1_vec)
errorbar(collect(p0_tuple) .- 1.0, value.(p1_vec), err.(p1_vec), fmt="d", mfc="none", color="blue")
fill_between(x0, value.(p1)+err.(p1), value.(p1)-err.(p1), alpha=0.4, color="green")

ylabel("p[1]")
setp(ax21.get_xticklabels(),visible=false) # Disable x tick labels

ax22 = fig.add_subplot(gs[3, 1])

uwerr.(p2_vec)
errorbar(collect(p0_tuple) .- 1.0, value.(p2_vec), err.(p2_vec), fmt="d", mfc="none", color="blue")
fill_between(x0, value.(p2)+err.(p2), value.(p2)-err.(p2), alpha=0.4, color="green")

ylabel("p[2]")
setp(ax22.get_xticklabels(),visible=false) # Disable x tick labels


ax3 = fig.add_subplot(gs[4, 1])

fill_between(x0, maximum(w)/2, maximum(w)/2, alpha=0.0, color="white")
PyPlot.plot(collect(p0_tuple) .- 1., w, linestyle="none", marker="o", mfc="none", color="blue")

ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))

setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels

ax4 = fig.add_subplot(gs[5, 1])

pval_vec = getfield.(fit_vec,:pval)
fill_between(x0, maximum(vcat(pval_vec...))/2, maximum(vcat(pval_vec...))/2, alpha=0.0, color="white")
PyPlot.plot(collect(p0_tuple) .- 1.0, pval_vec, linestyle="none", marker="o", mfc="none", color="blue")

ylabel(L"$\rm{p-values}$")

xlabel(L"$t/a$")


tight_layout()
display(fig)
close("all")

# corr and int plot

int = obs .* tmrw[1:sym_points]
int_rec = obs_rec .* tmrw
uwerr.(obs); uwerr.(int)
uwerr.(obs_rec); uwerr.(int_rec)

fig = figure(figsize=(8,12))
gs = fig.add_gridspec(2, 1)  # Adjust the height_ratios as needed

ax1 = fig.add_subplot(gs[1, 1])

expdiag = (diag == "LO") ? 2 : 3
res     = (alpha/π)^expdiag * sum(int) * 1e10; uwerr(res)
res_rec = (alpha/π)^expdiag * sum(int_rec) * 1e10; uwerr(res_rec)
dif = abs(res - res_rec); uwerr(dif)

title("$(ens.id) [$(key[2:3])-$(key[end-1:end]) impr $impr_set] : $(print_uwreal(res)) vs. $(print_uwreal(res_rec))  [$(print_uwreal(dif))]")

errorbar(collect(1:sym_points).-1, value.(obs), err.(obs), color="black", fmt="o", mfc="none" ,capsize=2 ,label="corr")
errorbar(t.-1, value.(obs_rec), err.(obs_rec), color="blue", fmt="s", mfc="none" ,capsize=2 ,label="rec corr")
yscale("log")
ylabel(L"G(t)")
legend()

setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels

ax2 = fig.add_subplot(gs[2, 1])

errorbar(collect(1:sym_points).-1, value.(int), err.(int), color="black", fmt="o", mfc="none" ,capsize=2 ,label="int")
errorbar(t.-1, value.(int_rec), err.(int_rec), color="blue", fmt="s", mfc="none" ,capsize=2 ,label="rec int")
legend()
xlabel("t/a")
ylabel(L"\tilde{K}(\hat{t})\times G(t)")

display(gcf())
close()

##

impr_set = "1"

path_bdio = path_bdio_dict["local"]

TMRD200 = BDIOread_TMR(path_bdio,"D200",resc=false,"1D",beta=false)
TMRD201 = BDIOread_TMR(path_bdio,"D201",resc=false,"1D",beta=false)

TMRD450 = BDIOread_TMR(path_bdio,"D450",resc=false,"1D",beta=false)
TMRD451 = BDIOread_TMR(path_bdio,"D451",resc=false,"1D",beta=false)

corrC101 = BDIOread_corr(path_bdio,"C101",impr_set,STD=false)
corrC102 = BDIOread_corr(path_bdio,"C102",impr_set,STD=false)

corrD200 = BDIOread_corr(path_bdio,"D200",impr_set,STD=false)
corrD201 = BDIOread_corr(path_bdio,"D201",impr_set,STD=false)

corrD450 = BDIOread_corr(path_bdio,"D450",impr_set,STD=false)
corrD451 = BDIOread_corr(path_bdio,"D451",impr_set,STD=false)

corrJ303 = BDIOread_corr(path_bdio,"J303",impr_set,STD=false)
corrJ304 = BDIOread_corr(path_bdio,"J304",impr_set,STD=false)

##

RC1 = (corrC101["g33_ll"] ./ corrC102["g33_ll"])[1:49]; uwerr.(RC1)
RD2 = (corrD200["g33_ll"] ./ corrD201["g33_ll"])[1:65]; uwerr.(RD2)
RD4 = (corrD450["g33_ll"] ./ corrD451["g33_ll"])[1:65]; uwerr.(RD4)
RJ3 = (corrJ303["g33_ll"] ./ corrJ304["g33_ll"])[1:97]; uwerr.(RJ3)

# errorbar(collect(1:length(RC1)), value.(RC1), err.(RC1), fmt="d", mfc="none", color="red" , label = "C101/C102")
errorbar(collect(1:length(RD2)).*asym(3.55).mean, value.(RD2), err.(RD2), fmt="o", mfc="none", color="black", label = "D200/D201")
errorbar(collect(1:length(RD4)).*asym(3.46).mean, value.(RD4), err.(RD4), fmt="o", mfc="none", color="blue" , label = "D450/D451")
# errorbar(collect(1:length(RJ3)), value.(RJ3), err.(RJ3), fmt="s", mfc="none", color="green", label = "J303/J304")
PyPlot.plot([1,length(RD2)], [1.0,1.0], ls="--", color="gray", alpha=0.5)
ylabel(L"G(t)/G^{m_s^{\rm{ph}}}(t)")
# xlabel(L"t/a")
xlabel(L"t/\rm{fm}")
xlim(0,3.0)
ylim(0.9,1.1)
legend()
display(gcf())
close()
