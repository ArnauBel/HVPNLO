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

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "..", "ObsBDIO"),
    # "local" => joinpath(julia_script_directory, "..", "..", "ObsCrosschecks", "Obs4LD(t0=Regensburg)", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "..", "ObsExternal", "PortableSSD", "ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_bPert   = joinpath(julia_script_directory, "..", "..", "PertSD")
path_FVCcont = joinpath(julia_script_directory, "..", "..", "FVCcont")

path_plot    = joinpath(julia_script_directory, "..", "..", "Slides & Plots", "Plots")

charge_factor = Dict(
    "33" => 1., "88" => 1/3., "88conn" => 1/3., "SS" => 1/9., "CC" => 4/9., "CCdisc" => 4/9., "C8disc" => 2/(3*sqrt(3)), "BB" => 1/9., "∆ls_amu" => 1/3., "∆ls_amuconn" => 1., "∆lc_b" => 4/9., "disc" => 1.,
    "33s" => "     ", "88s" => "(1/3)", "88conns" => "(1/3)", "SSs" => "(1/9)", "CCs" => "(4/9)", "CCdisc" => "(4/9)", "C8disc" => "2/(3*√3)", "BBs" => "(1/9)", "∆ls_amus" => "(1/3)", "∆ls_amuconns" => "", "∆lc_bs" => "(4/9)", "discs" => "     ",
    "3333" => 1.,  "3388" => 2/3., "33CC" => 8/9., "8888" => 1/9., "88CC" => 8/27., "CCCC" => 16/81., 
    "3333s" => "",  "3388s" => "(2/3)", "33CCs" => "(8/9)", "8888s" => "(1/9)", "88CCs" => "(8/27)", "CCCCs" => "(16/81)",
)

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

BL_factor = 1.45

BLIND_LD  = Any[true,BL_factor]

path_bdio = path_bdio_dict["local"]


Q33 = 5.0
QCC = 4.0

STD_DERIV = false
tl_IMPR   = true
VREF      = true
RESC      = false


amu      = Dict()
amusyst  = Dict()
AMU      = Dict()
AMUSYST  = Dict()
AMUt0ERR = Dict()
AMUERR   = Dict()
if VREF
    AMUFVC = Dict()
end

amu_bb = Dict(
    "NLOa"   => -0.233484,
    "NLOb"   => 0.0480212,
    "NLOa&b" => -0.185463
)

AMUgamma = Dict(
    "NLOa"   => uwreal([0.143 ,0.072],"HVPg"),
    "NLOb"   => uwreal([-0.121,0.060],"HVPg"),
    "NLOc"   => uwreal([-5.49 ,2.75].*1e-3,"HVPg"),
    "NLOa&b" => uwreal([2.22  ,1.11].*1e-2,"HVPg"),
    # "NLO"    => uwreal([1.67  ,0.84].*1e-2,"HVPg"),
)

AMU38 = Dict(
    "NLOa"   => uwreal([-7.45,3.73].*1e-2,"HVP38"),
    "NLOb"   => uwreal([5.19 ,2.60].*1e-2,"HVP38"),
    "NLOc"   => uwreal([2.80 ,1.40].*1e-3,"HVP38"),
    "NLOa&b" => uwreal([-2.26,1.13].*1e-2,"HVP38"),
    # "NLO"    => uwreal([-1.98,0.99].*1e-2,"HVP38"),
)

AMUIB = Dict()
[AMUIB[diag] = AMUgamma[diag] + AMU38[diag] for diag in keys(AMUgamma)]
AMUIB["NLO"] = AMUIB["NLOa&b"] + AMUIB["NLOc"]

MD_ph_prime, mDs_SU3 = BDIOread_mDs(path_bdio)

for diag in ["NLOa","NLOb","NLOa&b"]

    # For the SD window

    amu[diag]      = Dict("SD" => Dict(), "SDsub" => Dict(), "ID" => Dict(), "LD" => Dict())
    amusyst[diag]  = Dict("SD" => Dict(), "SDsub" => Dict(), "ID" => Dict(), "LD" => Dict())
    AMU[diag]      = Dict()
    AMUSYST[diag]  = Dict()
    AMUt0ERR[diag] = Dict()
    AMUERR[diag]   = Dict()
    if VREF
        AMUFVC[diag] = Dict()
    end

    amu33sub, info33sub = BDIOread_MAtot(path_bdio,diag,"SDsub","g33",StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,Q=Q33)
    b33Pert_Q33 = TXTread_bQ(path_bPert,diag)[Qlist .== Q33][1]; uwerr(b33Pert_Q33)
    b33Pert_QCC = TXTread_bQ(path_bPert,diag)[Qlist .== QCC][1]; uwerr(b33Pert_QCC)
    if VREF
        FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,"SDsub",Q=Q33)
        amu33sub += FVC_ChPT
        AMUFVC[diag]["SD"] = abs(0.1*FVC_ChPT)
    end
    amu[diag]["SDsub"]["33"] = amu33sub
    amu[diag]["SD"]["33"]    = amu33sub + Window("SD")(0) * b33Pert_Q33.mean
    amusyst[diag]["SDsub"]["33"] = info33sub["syst"]
    amusyst[diag]["SD"]["33"]    = sqrt(info33sub["syst"]^2 + Window("SD")(0)^2 * b33Pert_Q33.err^2)

    ∆ls_amu, info∆ls_amu      = BDIOread_MAtot(path_bdio,diag,"SDsub","∆ls_amu",StdDer=STD_DERIV)
    amu[diag]["SDsub"]["∆ls"] = ∆ls_amu
    amu[diag]["SD"]["88"]     = amu[diag]["SD"]["33"] + ∆ls_amu
    amusyst[diag]["SDsub"]["∆ls"] = info∆ls_amu["syst"]
    amusyst[diag]["SD"]["88"]     = sqrt(amusyst[diag]["SD"]["33"]^2 + info∆ls_amu["syst"]^2)

    ∆ls_amu_conn, info∆ls_amu_conn = BDIOread_MAtot(path_bdio,diag,"SDsub","∆ls_amuconn",StdDer=STD_DERIV)
    amu[diag]["SDsub"]["∆ls_conn"] = ∆ls_amu_conn
    amu[diag]["SD"]["SS"]          = 2*amu[diag]["SD"]["33"] + 3*∆ls_amu_conn
    amusyst[diag]["SDsub"]["∆ls_conn"] = info∆ls_amu_conn["syst"]
    amusyst[diag]["SD"]["SS"]          = sqrt(4*amusyst[diag]["SD"]["33"]^2 + 9*info∆ls_amu_conn["syst"]^2)

    amu[diag]["SD"]["disc"]     = 1/3. * amu[diag]["SD"]["88"] - 1/9. * (amu[diag]["SD"]["33"]+amu[diag]["SD"]["SS"])
    amusyst[diag]["SD"]["disc"] = 1/3. * sqrt(amusyst[diag]["SD"]["88"]^2 + 1/9. * (amusyst[diag]["SD"]["33"]^2 + amusyst[diag]["SD"]["SS"]^2))

    amuCCsub, infoCCsub = BDIOread_MAtot(path_bdio,diag,"SDsub","gCCconn",StdDer=STD_DERIV,Q=QCC)
    ∆lc_b, info∆lc_b    = BDIOread_MAtot(path_bdio,diag,"SDsub","∆lc_b",StdDer=STD_DERIV,Q=QCC)
    amu[diag]["SDsub"]["CC"]  = amuCCsub
    amu[diag]["SDsub"]["∆lc"] = ∆lc_b
    amu[diag]["SD"]["CC"] = amuCCsub + Window("SD")(0) * (2*b33Pert_QCC.mean + ∆lc_b); uwerr(amu[diag]["SD"]["CC"])
    der_mDs = mchist(amu[diag]["SD"]["CC"], "MD_ph [GeV]")[1] / artificial_err
    amu[diag]["SD"]["CC"]    += value(MD_ph - MD_ph_prime) * der_mDs
    amusyst[diag]["SDsub"]["CC"]  = infoCCsub["syst"]
    amusyst[diag]["SDsub"]["∆lc"] = info∆lc_b["syst"]
    amusyst[diag]["SD"]["CC"] = sqrt(infoCCsub["syst"]^2 + Window("SD")(0)^2 * (4*b33Pert_QCC.err^2 + info∆lc_b["syst"]^2))

    AMU[diag]["SD"]      = amu[diag]["SD"]["33"] + (1/3) * amu[diag]["SD"]["88"] + (4/9) * amu[diag]["SD"]["CC"] + (1/9) * amu_bb[diag]
    AMUSYST[diag]["SD"]  = sqrt(amusyst[diag]["SD"]["33"]^2 + 1/9 * amusyst[diag]["SD"]["88"]^2 + 16/81 * amusyst[diag]["SD"]["CC"]^2)
    AMUt0ERR[diag]["SD"] = get_t0err([AMU[diag]["SD"]],sqrtt0_ph_Madrid)[1]
    AMUERR[diag]["SD"]   = !VREF ? sqrt(AMU[diag]["SD"].err^2 + AMUSYST[diag]["SD"]^2 + AMUt0ERR[diag]["SD"]^2) : sqrt(AMU[diag]["SD"].err^2 + AMUSYST[diag]["SD"]^2 + AMUt0ERR[diag]["SD"]^2 + AMUFVC[diag]["SD"]^2)

    # Add botttom effects on it (only considered to affect the SD piece)
    amu[diag]["SD"]["BB"]     = AMU[diag]["BB"]     = uwreal(amu_bb[diag])
    amusyst[diag]["SD"]["BB"] = AMUSYST[diag]["BB"] = 0.5*amu_bb[diag]

    # Add charm-disconnected piece
    amuCCdisc, infoCCdisc = BDIOread_MAtot(path_bdio,diag,"SD","gCCdisc",StdDer=STD_DERIV)
    amuC8disc, infoC8disc = BDIOread_MAtot(path_bdio,diag,"SD","gC8disc",StdDer=STD_DERIV)
    amu[diag]["SD"]["CCdisc"] = AMU[diag]["CCdisc"] = amuCCdisc
    amu[diag]["SD"]["C8disc"] = AMU[diag]["C8disc"] = amuC8disc
    amusyst[diag]["SD"]["CCdisc"] = AMUSYST[diag]["CCdisc"] = infoCCdisc["syst"]
    amusyst[diag]["SD"]["C8disc"] = AMUSYST[diag]["C8disc"] = infoC8disc["syst"]

    # For the ID & LD windows

    for wind in ["ID","LD"]
        BLIND        = (wind=="LD") ? BLIND_LD[1] : false
        BLIND_factor = BLIND ? BLIND_LD[2] : 1.0


        amu33, info33 = BDIOread_MAtot(path_bdio,diag,wind,"g33",StdDer=STD_DERIV,BLIND=BLIND,Vref=VREF)
        if VREF
            FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,wind)
            AMUFVC[diag][wind] = abs(0.1*FVC_ChPT)
            amu33 += FVC_ChPT*BLIND_factor
        end
        amu[diag][wind]["33"]     = amu33/BLIND_factor
        amusyst[diag][wind]["33"] = info33["syst"]/BLIND_factor

        amu88, info88 = BDIOread_MAtot(path_bdio,diag,wind,"g88",StdDer=STD_DERIV,BLIND=BLIND)
        amu[diag][wind]["88"]     = amu88/BLIND_factor
        amusyst[diag][wind]["88"] = info88["syst"]/BLIND_factor

        amuSS, infoSS = BDIOread_MAtot(path_bdio,diag,wind,"gSS",StdDer=STD_DERIV,BLIND=BLIND)
        amu[diag][wind]["SS"]     = amuSS/BLIND_factor
        amusyst[diag][wind]["SS"] = infoSS["syst"]/BLIND_factor

        amu[diag][wind]["disc"]     = 1/3. * amu[diag][wind]["88"] - 1/9. * (amu[diag][wind]["33"]+amu[diag][wind]["SS"])
        amusyst[diag][wind]["disc"] = 1/3. * sqrt(amusyst[diag][wind]["88"]^2 + 1/9. * (amusyst[diag][wind]["33"]^2 + amusyst[diag][wind]["SS"]^2))

        amuCC, infoCC = BDIOread_MAtot(path_bdio,diag,wind,"gCCconn",StdDer=STD_DERIV,BLIND=false); uwerr(amuCC)
        der_mDs = mchist(amuCC, "MD_ph [GeV]")[1] / artificial_err
        amu[diag][wind]["CC"]     = amuCC + value(MD_ph - MD_ph_prime) * der_mDs
        amusyst[diag][wind]["CC"] = infoCC["syst"]


        AMU[diag][wind]      = (amu[diag][wind]["33"] + (1/3) * amu[diag][wind]["88"]) + (4/9) * amu[diag][wind]["CC"]
        AMUSYST[diag][wind]  = sqrt(amusyst[diag][wind]["33"]^2 + 1/9 * amusyst[diag][wind]["88"]^2 + 16/81 * amusyst[diag][wind]["CC"]^2)
        AMUt0ERR[diag][wind] = get_t0err([AMU[diag][wind]],sqrtt0_ph_Madrid)[1]
        AMUERR[diag][wind]   = !VREF ? sqrt(AMU[diag][wind].err^2 + AMUSYST[diag][wind]^2 + AMUt0ERR[diag][wind]^2) : sqrt(AMU[diag][wind].err^2 + AMUSYST[diag][wind]^2 + AMUt0ERR[diag][wind]^2 + AMUFVC[diag][wind]^2)
    end

    if VREF
        FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,"NW")
        AMUFVC[diag]["NW"] = abs(0.1*FVC_ChPT)
    end

    for comp in ["33","88","CC","SS","disc"]
        # BLIND_factor = (BLIND_LD[1] && comp != "CC") ? BLIND_LD[2] : 1.0

        AMU[diag][comp]      = amu[diag]["SD"][comp] + amu[diag]["ID"][comp] + amu[diag]["LD"][comp]
        AMUSYST[diag][comp]  = sqrt(amusyst[diag]["SD"][comp]^2 + amusyst[diag]["ID"][comp]^2 + amusyst[diag]["LD"][comp]^2)
        AMUt0ERR[diag][comp] = get_t0err([AMU[diag][comp]],sqrtt0_ph_Madrid)[1]
        AMUERR[diag][comp]   = (!VREF || comp != "33") ? sqrt(AMU[diag][comp].err^2 + AMUSYST[diag][comp]^2 + AMUt0ERR[diag][comp]^2) : sqrt(AMU[diag][comp].err^2 + AMUSYST[diag][comp]^2 + AMUt0ERR[diag][comp]^2 + AMUFVC[diag]["NW"]^2)
    end
end

AMU["NLOc"]      = Dict()
AMUSYST["NLOc"]  = Dict()
AMUt0ERR["NLOc"] = Dict()
AMUERR["NLOc"]   = Dict()

for comp in ["3333","3388","33CC","8888","88CC","CCCC"]
    AMU["NLOc"][comp], info = BDIOread_MAtot(path_bdio,"NLOc","NW","g"*comp,StdDer=STD_DERIV,BLIND=false,Vref=false); uwerr(AMU["NLOc"][comp])
    AMUSYST["NLOc"][comp]   = info["syst"]
    AMUt0ERR["NLOc"][comp]  = get_t0err([AMU["NLOc"][comp]],sqrtt0_ph_Madrid)[1]
    AMUERR["NLOc"][comp]    = sqrt(AMU["NLOc"][comp].err^2 + AMUSYST["NLOc"][comp]^2 + AMUt0ERR["NLOc"][comp]^2)
end

for comp in ["33CC","88CC","CCCC"]
    der_mDs = mchist(AMU["NLOc"][comp], "MD_ph [GeV]")[1] / artificial_err
    AMU["NLOc"][comp] += value(MD_ph - MD_ph_prime) * der_mDs
end


@info("Window-results :")
for diag in ["NLOa","NLOb","NLOa&b"]
    for wind in ["SD","ID","LD"]
        println("   (aµ[$diag])($wind) = $(print_uwreal(10*AMU[diag][wind],10*[AMUSYST[diag][wind],AMUt0ERR[diag][wind],AMUFVC[diag][wind]],total=true))")
    end
    println("--------------------------------------------------------")
end
println("\n")


@info("Channel-results :")
for diag in ["NLOa","NLOb","NLOa&b"]
    for ch in ["33","88","CC","SS","disc"]
        println("   $(charge_factor[ch*"s"])(aµ[$diag])($ch) = $(print_uwreal(10*charge_factor[ch]*AMU[diag][ch],10*charge_factor[ch]*[AMUSYST[diag][ch]],total=true))")
    end
    println("--------------------------------------------------------")
end
println("\n")

@info("Results in isoQCD :")
for diag in ["NLOa","NLOb","NLOa&b"]
    # AMU[diag]["tot"] = AMU[diag]["33"] + (1/3)*AMU[diag]["88"] + (4/9)*AMU[diag]["CC"]; print_uwreal(AMU[diag]["tot"])
    AMU[diag]["tot"]      = AMU[diag]["SD"] + AMU[diag]["ID"] + AMU[diag]["LD"]; uwerr(AMU[diag]["tot"])
    AMUSYST[diag]["tot"]  = sqrt(AMUSYST[diag]["SD"]^2 + AMUSYST[diag]["ID"]^2 + AMUSYST[diag]["LD"]^2)
    AMUt0ERR[diag]["tot"] = get_t0err([AMU[diag]["tot"]],sqrtt0_ph_Madrid)[1]
    AMUFVC[diag]["tot"]   = AMUFVC[diag]["NW"]
    AMUERR[diag]["tot"]   = sqrt(AMU[diag]["tot"].err^2 + AMUSYST[diag]["tot"]^2 + AMUt0ERR[diag]["tot"]^2 + AMUFVC[diag]["tot"]^2)

    println("   aµ[$diag] = $(print_uwreal(10*AMU[diag]["tot"],10*[AMUSYST[diag]["tot"],AMUt0ERR[diag]["tot"],AMUFVC[diag]["tot"]],total=true))")
end


AMU["NLOc"]["tot"] = uwreal(0.0); AMUSYSTNLOc2 = 0.0

for comp in ["3333","3388","33CC","8888","88CC","CCCC"]
    AMU["NLOc"]["tot"] += charge_factor[comp]*AMU["NLOc"][comp]
    AMUSYSTNLOc2       += (charge_factor[comp]*AMUSYST["NLOc"][comp])^2
end

AMUSYST["NLOc"]["tot"]  = sqrt(AMUSYSTNLOc2)
AMUt0ERR["NLOc"]["tot"] = get_t0err([AMU["NLOc"]["tot"]],sqrtt0_ph_Madrid)[1]
AMUERR["NLOc"]["tot"]   = sqrt(AMU["NLOc"]["tot"].err^2 + AMUSYST["NLOc"]["tot"]^2 + AMUt0ERR["NLOc"]["tot"]^2)

println("   aµ[NLOc] = $(print_uwreal(10*AMU["NLOc"]["tot"],10*[AMUSYST["NLOc"]["tot"],AMUt0ERR["NLOc"]["tot"]],total=true))")

println("\n")
println("--------------------------------------------------------")
println("--------------------------------------------------------")

AMU["NLO"]      = Dict()
AMUSYST["NLO"]  = Dict()
AMUt0ERR["NLO"] = Dict()
AMUFVC["NLO"]   = Dict()
AMUERR["NLO"]   = Dict()

AMU["NLO"]["tot"]      = AMU["NLOa&b"]["tot"] + AMU["NLOc"]["tot"]
AMUSYST["NLO"]["tot"]  = sqrt(AMUSYST["NLOa&b"]["tot"]^2 + AMUSYST["NLOc"]["tot"]^2)
AMUt0ERR["NLO"]["tot"] = get_t0err([AMU["NLO"]["tot"]],sqrtt0_ph_Madrid)[1]
AMUFVC["NLO"]["tot"]   = AMUFVC["NLOa&b"]["tot"]
AMUERR["NLO"]["tot"]   = sqrt(AMU["NLO"]["tot"].err^2 + AMUSYST["NLO"]["tot"]^2 + AMUt0ERR["NLO"]["tot"]^2 + AMUFVC["NLO"]["tot"]^2)

println("   aµ[NLO] = $(print_uwreal(10*AMU["NLO"]["tot"],10*[AMUSYST["NLO"]["tot"],AMUt0ERR["NLO"]["tot"],AMUFVC["NLO"]["tot"]],total=true))")

println("\n")
@info("Final results :")

for diag in ["NLOa","NLOb","NLOa&b","NLOc","NLO"]
    if diag == "NLO"
            println("--------------------------------------------------------")
    end
    uwerr(AMUIB[diag])
    if diag != "NLOc"
        println(
            "   aµ[$diag] = $(print_uwreal(
            10*(AMU[diag]["tot"]+AMUIB[diag].mean),
            10*[AMUSYST[diag]["tot"],AMUt0ERR[diag]["tot"],AMUFVC[diag]["tot"],AMUIB[diag].err],
            total=true))"
            )
    else
        println(
            "   aµ[$diag] = $(print_uwreal(
            10*(AMU[diag]["tot"]+AMUIB[diag].mean),
            10*[AMUSYST[diag]["tot"],AMUt0ERR[diag]["tot"],AMUIB[diag].err],
            total=true))"
            )
    end
end

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

# Comparison plots

SAVE     = false
OVERSAVE = false

RES = Dict(
    "NLOa" => [
        [uwreal([-20.73 ,0.19 ],"HLMNT04 NLOa"),"HLMNT04"],
        [uwreal([-20.90 ,0.21 ],"KLMS14"),"KLMS14"],
        [uwreal([-20.613,0.130],"Jeger NLOa"),"Jegerlehner"],
        [uwreal([-20.77 ,0.08 ],"kesh NLOa"),"KNT18"],
        ],
    "NLOb" => [
        [uwreal([10.60 ,0.10 ],"HLMNT04 NLOb"),"HLMNT04"],
        [uwreal([10.68 ,0.11 ],"KLMS14 NLOb"),"KLMS14"],
        [uwreal([10.349,0.063],"Jeger NLOb"),"Jegerlehner"],
        [uwreal([10.62 ,0.04 ],"kesh NLOb"),"KNT18"],
        ],
    "NLOc" => [
        [uwreal([0.34 ,0.01 ],"HLMNT04 NLOc"),"HLMNT04"],
        # [uwreal([0.35 ,0.014],"KLMS14 NLOc"),"KLMS14 *"],
        [uwreal([0.337,0.005],"Jeger NLOc"),"Jegerlehner"],
        [uwreal([0.34 ,0.01 ],"kesh NLOc"),"KNT18"],
        ],
    "NLO"  => [
        [uwreal([-9.79 ,0.10 ],"HLMNT04 NLO"),"HLMNT04"],
        [uwreal([-9.84 ,0.07 ],"HLMNT11 NLO"),"HLMNT11"],
        [uwreal([-9.87 ,0.09 ],"KLMS14 NLO"),"KLMS14"],
        [uwreal([-9.927,0.067],"Jeger NLO" ),"Jegerlehner"],
        [uwreal([-9.82 ,0.04 ],"kesh NLO" ),"KNT18"],
        ]
)

list_old = Dict(
    "NLOa"   => ["KNT18"],
    "NLOb"   => ["KNT18"],
    "NLOc"   => ["KNT18"],
    "NLO"    => ["HLMNT04","KNT18"]
)

RES_sl = Dict(
    "NLOa"   => [[uwreal([-21.5,0.15],"Spacelike"),"Spacelike - Mainz 2025"]],
    "NLOb"   => [[uwreal([11,0.1],"Spacelike"),"Spacelike - Mainz 2025"]],
    "NLOc"   => [
        # [uwreal([0.39710,0.00966],"Spacelike"),"Spacelike - Mainz 2022"],
        # [uwreal([0.37685,0.00837],"Spacelike"),"Spacelike - Mainz 2025"],
        [uwreal([0.3735,0.0081],"Spacelike"),"Spacelike - Mainz 2025"],
        ],
    "NLO"    => [[uwreal([-10.5,0.1],"Spacelike"),"Spacelike - Mainz 2025"]],
)

# RES_WP = Dict("NLO" => [uwreal([-9.83,0.07],"WP25"),"WP 2025"])
RES_WP = Dict(
    "NLOa" => Dict(
        "part" => [
            [uwreal([-20.75 ,0.07 ],"KNT19 NLO"),"KNT19"],
            [uwreal([-21.34 ,0.13 ],"KNT19/CMD3 NLO"),"KNT19/CMD-3"]
        ],
        "aver" => [uwreal([-21.045,0.295],"WP25"),"WP 2025"]
    ),
    "NLOb" => Dict(
        "part" => [
            [uwreal([10.59 ,0.04 ],"KNT19 NLO"),"KNT19"],
            [uwreal([10.92 ,0.07 ],"KNT19/CMD3 NLO"),"KNT19/CMD-3"]
        ],
        "aver" => [uwreal([10.755,0.165],"WP25"),"WP 2025"]
    ),
    "NLOc" => Dict(
        "part" => [
            [uwreal([0.34 ,0.01 ],"KNT19 NLO"),"KNT19"],
            [uwreal([0.36 ,0.01 ],"KNT19/CMD3 NLO"),"KNT19/CMD-3"]
        ],
        "aver" => [uwreal([0.35,0.01],"WP25"),"WP 2025"]
    ),
    "NLO" => Dict(
        "part" => [
            [uwreal([-9.83 ,0.04 ],"KNT19 NLO"),"KNT19"],
            [uwreal([-10.08,0.06 ],"KNT19/CMD3 NLO"),"KNT19/CMD-3"]
        ],
        "aver" => [uwreal([-9.96,0.13],"WP25"),"WP 2025"]
    )
)

for diag in ["NLOa","NLOb","NLOc","NLO"]
    PyPlot.title(diag)
    y_ticks = ["this work"]
    errorbar(10*(AMU[diag]["tot"].mean+AMUIB[diag].mean), xerr=10*AMU[diag]["tot"].err, 0.25, 0.0, fmt="o", color="black", ms=10, capsize=2)
    errorbar(10*(AMU[diag]["tot"].mean+AMUIB[diag].mean), xerr=10*sqrt(AMUERR[diag]["tot"]^2+AMUIB[diag].err^2), 0.25, 0.0, fmt="o", color="black", ms=10, capsize=2)
    errorbar(10*AMU[diag]["tot"].mean, xerr=10*AMU[diag]["tot"].err, -0.25, 0.0, mfc="none", fmt="o", color="black", ms=10, capsize=2)
    errorbar(10*AMU[diag]["tot"].mean, xerr=10*AMUERR[diag]["tot"] , -0.25, 0.0, mfc="none", fmt="o", color="black", ms=10, capsize=2)
    i = 0
    for res in reverse(RES_sl[diag])
        i -= 1
        uwerr(res[1])
        errorbar(10*(res[1].mean+AMUIB[diag].mean), xerr=10*sqrt(res[1].err^2+AMUIB[diag].err^2), i, 0.0, fmt="d", color="green", ms=10, capsize=2)
        push!(y_ticks,res[2])
    end
    i -= 1
    res = RES_WP[diag]["aver"]
    uwerr(res[1])
    errorbar(10*res[1].mean, xerr=10*res[1].err, i, 0.0, fmt="s", color="red", ms=10, capsize=2)
    push!(y_ticks,res[2])
    for res in reverse(RES_WP[diag]["part"])
        i -= 1
        uwerr(res[1])
        errorbar(10*res[1].mean, xerr=10*res[1].err, i, 0.0, fmt="s", color="purple", ms=10, capsize=2)
        push!(y_ticks,res[2])
    end
    for res in reverse(RES[diag])
        i -= 1
        uwerr(res[1])
        mfc = (res[2] in list_old[diag]) ? "none" : "blue"
        errorbar(10*res[1].mean, xerr=10*res[1].err, i, 0.0, fmt="o", mfc=mfc, color="blue", ms=10, capsize=2)
        push!(y_ticks,res[2])
    end
    fill_betweenx([1,i-1], 
        10*((AMU[diag]["tot"].mean+AMUIB[diag].mean)-sqrt(AMUERR[diag]["tot"]^2+AMUIB[diag].err^2)), 
        10*((AMU[diag]["tot"].mean+AMUIB[diag].mean)+sqrt(AMUERR[diag]["tot"]^2+AMUIB[diag].err^2)), 
        color="gray", alpha=0.4
    )
    
    xlabel(latexstring("a_{\\mu}^{\\rm{hvp}}[\\rm{$(diag)}]\\times10^{11}"))
    PyPlot.yticks(reverse(collect(i:0)), y_ticks, rotation = 30, fontsize=15)

    tight_layout()
    display(gcf())
    if SAVE
        p = create_path(path_plot,["Results","Res_$diag.pdf"],OVERWRITE=OVERSAVE)
        PyPlot.savefig(p)
    end
    close()
end

# All together

# Create 1×4 axes, shared y axis, minimal spacing
fig, axs = subplots(
    1, 4,
    figsize=(14, 6),
    sharey=true,
    gridspec_kw=Dict("wspace" => 0.05)   # small gap between plots
)

# Example of filling each subplot
for (j, diag) in enumerate(["NLOa", "NLOb", "NLOc", "NLO"])
    ax = axs[j]

    # ax.set_title(diag, fontsize=16)

    y_ticks = ["this work"]
    ax.errorbar(10*(AMU[diag]["tot"].mean+AMUIB[diag].mean), xerr=10*AMU[diag]["tot"].err, 0.25, 0.0, fmt="o", color="black", ms=10, capsize=2)
    ax.errorbar(10*(AMU[diag]["tot"].mean+AMUIB[diag].mean), xerr=10*sqrt(AMUERR[diag]["tot"]^2+AMUIB[diag].err^2), 0.25, 0.0, fmt="o", color="black", ms=10, capsize=2)
    ax.errorbar(10*AMU[diag]["tot"].mean, xerr=10*AMU[diag]["tot"].err, -0.25, 0.0, mfc="none", fmt="o", color="black", ms=10, capsize=2)
    ax.errorbar(10*AMU[diag]["tot"].mean, xerr=10*AMUERR[diag]["tot"] , -0.25, 0.0, mfc="none", fmt="o", color="black", ms=10, capsize=2)
    i = 0
    # for res in reverse(RES_sl[diag])
    #     i -= 1
    #     uwerr(res[1])
    #     ax.errorbar(10*(res[1].mean+AMUIB[diag].mean), xerr=10*sqrt(res[1].err^2+AMUIB[diag].err^2), i, 0.0, fmt="d", color="green", ms=10, capsize=2)
    #     push!(y_ticks,res[2])
    # end
    i -= 1
    res = RES_WP[diag]["aver"]
    uwerr(res[1])
    ax.errorbar(10*res[1].mean, xerr=10*res[1].err, i, 0.0, fmt="s", color="red", ms=10, capsize=2)
    push!(y_ticks,res[2])
    for res in reverse(RES_WP[diag]["part"])
        i -= 1
        uwerr(res[1])
        ax.errorbar(10*res[1].mean, xerr=10*res[1].err, i, 0.0, fmt="s", color="purple", ms=10, capsize=2)
        push!(y_ticks,res[2])
    end
    RES_NLO_str = getindex.(reverse(RES["NLO"]),2)
    for (l,res) in enumerate(reverse(RES[diag]))
        k=1
        while res[2] != RES_NLO_str[l+k-1]
            k+=1
        end
        i -= k
        uwerr(res[1])
        mfc = (res[2] in list_old[diag]) ? "none" : "blue"
        ax.errorbar(10*res[1].mean, xerr=10*res[1].err, i, 0.0, fmt="o", mfc=mfc, color="blue", ms=10, capsize=2)
        # push!(y_ticks,res[2])
    end
    ax.fill_betweenx([1,i-1], 
        10*((AMU[diag]["tot"].mean+AMUIB[diag].mean)-sqrt(AMUERR[diag]["tot"]^2+AMUIB[diag].err^2)), 
        10*((AMU[diag]["tot"].mean+AMUIB[diag].mean)+sqrt(AMUERR[diag]["tot"]^2+AMUIB[diag].err^2)), 
        color="gray", alpha=0.4
    )
    
    ax.set_xlabel(latexstring("a_{\\mu}^{\\rm{hvp}}[\\rm{$(diag)}]\\times10^{11}"), fontsize=16)
    ax.set_yticks(reverse(collect(i:0)), vcat(y_ticks,RES_NLO_str), rotation = 30, fontsize=14)

    if j == 1
        ax.set_yticklabels(vcat(y_ticks,RES_NLO_str), rotation=30, fontsize=14)
    else
        # ax.set_yticklabels([])   # remove y labels for right panels
    end

    ax.grid(true)
end

# tight_layout()
display(fig)
if SAVE
    p = create_path(path_plot,["Results","Res_joint.pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p)
end
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

# Pie-chart for the central value

SAVE     = false
OVERSAVE = false

clamp01(x) = max(0.0, min(1.0, float(x)))

"""
    scale_color(base::Tuple, factor::Real)

Returns an RGB tuple scaled by `factor` and clamped to [0,1].
`base` should be a 3-tuple (r,g,b) with each in 0..1.
This corresponds to Mathematica's `scaleColor[base, factor]`.
"""
function scale_color(base::Tuple{<:Real, <:Real, <:Real}, factor::Real)
    return (clamp01(base[1]*factor), clamp01(base[2]*factor), clamp01(base[3]*factor))
end

for diag in ["NLOa","NLOb","NLOa&b"]
    res = Dict()
    res["SD"] = [charge_factor[key]*amu[diag]["SD"][key].mean for key in ["33","88","CC","BB"]]
    res["ID"] = [charge_factor[key]*amu[diag]["ID"][key].mean for key in ["33","88","CC"]]
    res["LD"] = [(charge_factor[key]/BLIND_LD[2])*amu[diag]["LD"][key].mean for key in ["33","88","CC"]]

    baseColors = Dict(
        "SD" => (0.2,0.6,0.9), 
        "ID" => (0.9,0.5,0.2),
        "LD" => (0.4,0.8,0.3)
        )

    scales = (1/1.2, 1/1.4, 1/1.6, 1/1.8)
    innerColors = reduce(vcat, [baseColors[key] for key in ["SD","ID","LD"]])
    outerColors = reduce(vcat, [[baseColors[key].*scales[i] for i=1:length(res[key])] for key in ["SD","ID","LD"]])

    innerValues = [sum(res[key]) for key in ["SD","ID","LD"]]
    outerValues = vcat([res[key] for key in ["SD","ID","LD"]]...)

    fig = figure(figsize=(6,4))
    ax = fig.add_subplot(1,1,1)

    # Outer ring: radius 1.0, width 0.3
    ax.pie(abs.(reverse(outerValues)),
            radius=1.0,
            colors=reverse(outerColors),
            startangle=180,
            wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"))

    # Inner ring: radius 0.65, width 0.35
    ax.pie(abs.(reverse(innerValues)),
            radius=0.68,
            colors=reverse(innerColors),
            startangle=180,
            wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"),
            labels=nothing)

    diag_str = diag == "NLOa&b" ? "NLOa\\&b" : diag
    ax.text(0.0, 0.0, latexstring("\\rm{$(diag_str)}"), ha="center", va="center", fontsize=12, fontweight="bold", color="black")
    ax.set(aspect="equal")

    # legend_colors = vcat(innerColors, [(0.9,0.9,0.9), (0.8,0.8,0.8), (0.7,0.7,0.7), (0.6,0.6,0.6)])
    # legend_labels = vcat(["SD","ID","LD"], ["3,3","8,8","c,c","b,b"])
    legend_colors = vcat(innerColors, [(0.9,0.9,0.9), (0.8,0.8,0.8), (0.7,0.7,0.7)])
    legend_labels = vcat(["SD","ID","LD"], ["3,3","8,8","c,c"])

    mpatches = PyPlot.matplotlib[:patches]
    handles = [mpatches.Patch(facecolor=c) for c in legend_colors]
    ax.legend(handles, legend_labels, loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=2, frameon=false)

    display(gcf())
    if SAVE
        p = create_path(path_plot,["Results","PieChart_CV_$diag.pdf"],OVERWRITE=OVERSAVE)
        PyPlot.savefig(p, bbox_inches="tight")
    end
    close()  # avoid showing/accumulating figures in loops
end

# For NLOc

fig = figure(figsize=(6,4))
ax = fig.add_subplot(1,1,1)

res_vec = value.([charge_factor[comp]*AMU["NLOc"][comp] for comp in ["3333","3388","33CC","8888","88CC","CCCC"]])

ax.pie( res_vec,
        radius=1.0,
        startangle=180,
        wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"))

ax.text(0.0, 0.0, latexstring("\\rm{NLOc}"), ha="center", va="center", fontsize=12, fontweight="bold", color="black")
ax.set(aspect="equal")

ax.legend(["3,3-3,3","3,3-8,8","3,3-c,c","8,8-8,8","8,8-c,c","c,c-c,c"], loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=2, frameon=false)

display(gcf())
if SAVE
    p = create_path(path_plot,["Results","PieChart_CV_NLOc.pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p, bbox_inches="tight")
end
close()

# All together

fig, axs = subplots(
    1, 3,
    figsize=(14, 6),
    sharey=true,
    gridspec_kw=Dict("wspace" => 0.05)   # small gap between plots
)

for (j, diag) in enumerate(["NLOa", "NLOb", "NLOa&b"])
    ax = axs[j]

    res = Dict()
    res["SD"] = [charge_factor[key]*amu[diag]["SD"][key].mean for key in ["33","88","CC","BB"]]
    res["ID"] = [charge_factor[key]*amu[diag]["ID"][key].mean for key in ["33","88","CC"]]
    res["LD"] = [(charge_factor[key]/BLIND_LD[2])*amu[diag]["LD"][key].mean for key in ["33","88","CC"]]

    baseColors = Dict(
        "SD" => (0.2,0.6,0.9), 
        "ID" => (0.9,0.5,0.2),
        "LD" => (0.4,0.8,0.3)
        )

    scales = (1/1.2, 1/1.4, 1/1.6, 1/1.8)
    innerColors = reduce(vcat, [baseColors[key] for key in ["SD","ID","LD"]])
    outerColors = reduce(vcat, [[baseColors[key].*scales[i] for i=1:length(res[key])] for key in ["SD","ID","LD"]])

    innerValues = [sum(res[key]) for key in ["SD","ID","LD"]]
    outerValues = vcat([res[key] for key in ["SD","ID","LD"]]...)

    # Outer ring: radius 1.0, width 0.3
    ax.pie(abs.(reverse(outerValues)),
            radius=1.0,
            colors=reverse(outerColors),
            startangle=180,
            wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"))

    # Inner ring: radius 0.65, width 0.35
    ax.pie(abs.(reverse(innerValues)),
            radius=0.68,
            colors=reverse(innerColors),
            startangle=180,
            wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"),
            labels=nothing)

    diag_str = diag == "NLOa&b" ? "NLOa\\&b" : diag
    ax.text(0.0, 0.0, latexstring("\\rm{$(diag_str)}"), ha="center", va="center", fontsize=12, fontweight="bold", color="black")
    ax.set(aspect="equal")

    # legend_colors = vcat(innerColors, [(0.9,0.9,0.9), (0.8,0.8,0.8), (0.7,0.7,0.7), (0.6,0.6,0.6)])
    # legend_labels = vcat(["SD","ID","LD"], ["3,3","8,8","c,c","b,b"])
    legend_colors = vcat(innerColors, [(0.9,0.9,0.9), (0.8,0.8,0.8), (0.7,0.7,0.7)])
    legend_labels = vcat(["SD","ID","LD"], ["3,3","8,8","c,c"])

    mpatches = PyPlot.matplotlib[:patches]
    handles = [mpatches.Patch(facecolor=c) for c in legend_colors]
    if j == 3
        ax.legend(handles, legend_labels, loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=2, frameon=false)
    end
end
display(gcf())
if SAVE
    p = create_path(path_plot,["Results","PieChart_CV_joint.pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p, bbox_inches="tight")
end
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

# Pie-chart for the variance

SAVE     = false
OVERSAVE = false

for diag in ["NLOa","NLOb","NLOa&b"]
    var = Dict()
    for wind in ["SD","ID","LD"]
        var[wind] = [
            AMU[diag][wind].err^2,
            AMUSYST[diag][wind]^2,
            AMUt0ERR[diag][wind]^2,
            AMUFVC[diag][wind]^2
            ]
    end

    baseColors = Dict(
        "SD" => (0.2,0.6,0.9), 
        "ID" => (0.9,0.5,0.2),
        "LD" => (0.4,0.8,0.3)
        )

    scales = (1/1.2, 1/1.4, 1/1.6, 1/1.8)
    innerColors = reduce(vcat, [baseColors[key] for key in ["SD","ID","LD"]])
    outerColors = reduce(vcat, [[baseColors[key].*scales[i] for i=1:length(var[key])] for key in ["SD","ID","LD"]])

    innerValues = [sum(var[key]) for key in ["SD","ID","LD"]]
    outerValues = vcat([var[key] for key in ["SD","ID","LD"]]...)

    push!(innerColors, (0.7,0.0,0.0))
    push!(outerColors, (0.7,0.0,0.0))

    push!(innerValues, AMUIB[diag].err^2)
    push!(outerValues, AMUIB[diag].err^2)

    fig = figure(figsize=(6,4))
    ax = fig.add_subplot(1,1,1)

    # Outer ring: radius 1.0, width 0.3
    ax.pie(abs.(reverse(outerValues)),
            radius=1.0,
            colors=reverse(outerColors),
            startangle=180,
            wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"))

    # Inner ring: radius 0.65, width 0.35
    ax.pie(abs.(reverse(innerValues)),
            radius=0.68,
            colors=reverse(innerColors),
            startangle=180,
            wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"),
            labels=nothing)

    diag_str = diag == "NLOa&b" ? "NLOa\\&b" : diag
    ax.text(0.0, 0.0, latexstring("\\rm{$(diag_str)}"), ha="center", va="center", fontsize=12, fontweight="bold", color="black")
    ax.set(aspect="equal")

    legend_colors = vcat(innerColors, [(0.9,0.9,0.9),(0.8,0.8,0.8),(0.7,0.7,0.7),(0.6,0.6,0.6)])
    legend_labels = vcat(["SD","ID","LD","IB"], ["stat.","syst.","scale","fvc"])
    # legend_labels = vcat(["SD","ID","LD"], ["stat.","syst.","scale","fvc"])

    mpatches = PyPlot.matplotlib[:patches]
    handles = [mpatches.Patch(facecolor=c) for c in legend_colors]

    ax.legend(handles, legend_labels, loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=2, frameon=false)

    display(gcf())
    if SAVE
        p = create_path(path_plot,["Results","PieChart_err_$diag.pdf"],OVERWRITE=OVERSAVE)
        PyPlot.savefig(p, bbox_inches="tight")
    end
    close()  # avoid showing/accumulating figures in loops
end

# For NLOc

fig = figure(figsize=(6,4))
ax = fig.add_subplot(1,1,1)

var_vec = [
    AMU["NLOc"]["tot"].err^2,
    AMUSYST["NLOc"]["tot"]^2,
    AMUt0ERR["NLOc"]["tot"]^2,
    AMUIB["NLOc"].err^2
    ]

ax.pie( var_vec,
        radius=1.0,
        startangle=180,
        wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"))

ax.text(0.0, 0.0, latexstring("\\rm{NLOc}"), ha="center", va="center", fontsize=12, fontweight="bold", color="black")
ax.set(aspect="equal")

ax.legend(["stat.","syst.","scale","IB"], loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=1, frameon=false)

display(gcf())
if SAVE
    p = create_path(path_plot,["Results","PieChart_err_NLOc.pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p, bbox_inches="tight")
end
close()

# All together

fig, axs = subplots(
    1, 3,
    figsize=(14, 6),
    sharey=true,
    gridspec_kw=Dict("wspace" => 0.05)   # small gap between plots
)

for (j, diag) in enumerate(["NLOa", "NLOb", "NLOa&b"])
    ax = axs[j]

    var = Dict()
    for wind in ["SD","ID","LD"]
        var[wind] = [
            AMU[diag][wind].err^2,
            AMUSYST[diag][wind]^2,
            AMUt0ERR[diag][wind]^2,
            AMUFVC[diag][wind]^2
            ]
    end

    baseColors = Dict(
        "SD" => (0.2,0.6,0.9), 
        "ID" => (0.9,0.5,0.2),
        "LD" => (0.4,0.8,0.3)
        )

    scales = (1/1.2, 1/1.4, 1/1.6, 1/1.8)
    innerColors = reduce(vcat, [baseColors[key] for key in ["SD","ID","LD"]])
    outerColors = reduce(vcat, [[baseColors[key].*scales[i] for i=1:length(var[key])] for key in ["SD","ID","LD"]])

    innerValues = [sum(var[key]) for key in ["SD","ID","LD"]]
    outerValues = vcat([var[key] for key in ["SD","ID","LD"]]...)

    push!(innerColors, (0.7,0.0,0.0))
    push!(outerColors, (0.7,0.0,0.0))

    push!(innerValues, AMUIB[diag].err^2)
    push!(outerValues, AMUIB[diag].err^2)

    # Outer ring: radius 1.0, width 0.3
    ax.pie(abs.(reverse(outerValues)),
            radius=1.0,
            colors=reverse(outerColors),
            startangle=180,
            wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"))

    # Inner ring: radius 0.65, width 0.35
    ax.pie(abs.(reverse(innerValues)),
            radius=0.68,
            colors=reverse(innerColors),
            startangle=180,
            wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"),
            labels=nothing)

    diag_str = diag == "NLOa&b" ? "NLOa\\&b" : diag
    ax.text(0.0, 0.0, latexstring("\\rm{$(diag_str)}"), ha="center", va="center", fontsize=12, fontweight="bold", color="black")
    ax.set(aspect="equal")

    legend_colors = vcat(innerColors, [(0.9,0.9,0.9),(0.8,0.8,0.8),(0.7,0.7,0.7),(0.6,0.6,0.6)])
    legend_labels = vcat(["SD","ID","LD","IB"], ["stat.","syst.","scale","fvc"])
    # legend_labels = vcat(["SD","ID","LD"], ["stat.","syst.","scale","fvc"])

    mpatches = PyPlot.matplotlib[:patches]
    handles = [mpatches.Patch(facecolor=c) for c in legend_colors]

    if j == 3
        ax.legend(handles, legend_labels, loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=2, frameon=false)
    end
end
display(gcf())
if SAVE
    p = create_path(path_plot,["Results","PieChart_err_joint.pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p, bbox_inches="tight")
end
close()  # avoid showing/accumulating figures in loops

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

# LO pie-chart

SAVE     = false
OVERSAVE = false

clamp01(x) = max(0.0, min(1.0, float(x)))

"""
    scale_color(base::Tuple, factor::Real)

Returns an RGB tuple scaled by `factor` and clamped to [0,1].
`base` should be a 3-tuple (r,g,b) with each in 0..1.
This corresponds to Mathematica's `scaleColor[base, factor]`.
"""
function scale_color(base::Tuple{<:Real, <:Real, <:Real}, factor::Real)
    return (clamp01(base[1]*factor), clamp01(base[2]*factor), clamp01(base[3]*factor))
end

amu_LO = Dict(
    "SD" => Dict(
        "33" => uwreal([43.06,0.22],"SD-33"),
        "88" => uwreal(3*[13.857,0.081],"SD-88"),
        "CC" => uwreal(9/4*[11.53,0.30],"SD-CC"),
        "BB" => uwreal(9*[0.29,0.03],"SD-BB")
    ),
    "ID" => Dict(
        "33" => uwreal([186.3,1.1],"ID-33"),
        "88" => uwreal(3*[47.41,0.29],"ID-88"),
        "CC" => uwreal(9/4*[2.89,0.14],"ID-CC"),
    ),
    "LD" => Dict(
        "33" => uwreal([378.7,4.8],"LD-33"),
        "88" => uwreal(3*[44.5,1.6],"LD-88"),
        "CC" => uwreal(9/4*[0.01409,0.00069],"LD-CC"),
    ),
)

res = Dict()
res["SD"] = [charge_factor[key]*amu_LO["SD"][key].mean for key in ["33","88","CC","BB"]]
res["ID"] = [charge_factor[key]*amu_LO["ID"][key].mean for key in ["33","88","CC"]]
res["LD"] = [charge_factor[key]*amu_LO["LD"][key].mean for key in ["33","88","CC"]]

baseColors = Dict(
    "SD" => (0.2,0.6,0.9), 
    "ID" => (0.9,0.5,0.2),
    "LD" => (0.4,0.8,0.3)
    )

scales = (1/1.2, 1/1.4, 1/1.6, 1/1.8)
innerColors = reduce(vcat, [baseColors[key] for key in ["SD","ID","LD"]])
outerColors = reduce(vcat, [[baseColors[key].*scales[i] for i=1:length(res[key])] for key in ["SD","ID","LD"]])

innerValues = [sum(res[key]) for key in ["SD","ID","LD"]]
outerValues = vcat([res[key] for key in ["SD","ID","LD"]]...)

fig = figure(figsize=(6,4))
ax = fig.add_subplot(1,1,1)

# Outer ring: radius 1.0, width 0.3
ax.pie(abs.(reverse(outerValues)),
        radius=1.0,
        colors=reverse(outerColors),
        startangle=180,
        wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"))

# Inner ring: radius 0.65, width 0.35
ax.pie(abs.(reverse(innerValues)),
        radius=0.68,
        colors=reverse(innerColors),
        startangle=180,
        wedgeprops=Dict("width"=>0.3, "edgecolor"=>"white"),
        labels=nothing)

ax.text(0.0, 0.0, latexstring("\\rm{LO}"), ha="center", va="center", fontsize=12, fontweight="bold", color="black")
ax.set(aspect="equal")

# legend_colors = vcat(innerColors, [(0.9,0.9,0.9), (0.8,0.8,0.8), (0.7,0.7,0.7), (0.6,0.6,0.6)])
# legend_labels = vcat(["SD","ID","LD"], ["3,3","8,8","c,c","b,b"])
legend_colors = vcat(innerColors, [(0.9,0.9,0.9), (0.8,0.8,0.8), (0.7,0.7,0.7)])
legend_labels = vcat(["SD","ID","LD"], ["3,3","8,8","c,c"])

mpatches = PyPlot.matplotlib[:patches]
handles = [mpatches.Patch(facecolor=c) for c in legend_colors]
ax.legend(handles, legend_labels, loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=2, frameon=false)

display(gcf())
if SAVE
    p = create_path(path_plot,["Results","PieChart_CV_LO.pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p, bbox_inches="tight")
end
close() 

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

# Charm sector: lattice vs. pQCD

AMU_pQCD = Dict(
    "NLOa" => uwreal([-8.034,0.127],"NLOa CC-pQCD"),
    "NLOb" => uwreal([2.247,0.038],"NLOb CC-pQCD")
)
AMU_pQCD["NLOa&b"] = uwreal([-5.787,0.089],"NLOab CC-pQCD")
# AMU_pQCD["NLOa&b"] = AMU_pQCD["NLOa"] + AMU_pQCD["NLOb"]

# [uwerr(AMU_pQCD[key]) for key in keys(AMU_pQCD)]

fig, (ax1, ax2, ax3) = subplots(1, 3,
    gridspec_kw = Dict("width_ratios" => [1, 1, 1], "wspace" => 0.1),
    figsize = (14, 2)
)

factor = 10*4/9

@info("pQCD vs. lattice in the charm-quark sector :")
println("-------------------------------------------------")
for diag in ["NLOa","NLOb","NLOa&b"]
    println(
        "   aµ[$diag](CC-lat)  = $(print_uwreal(
        factor*AMU[diag]["CC"],
        factor*[AMUSYST[diag]["CC"],AMUt0ERR[diag]["CC"]],
        total=true))"
        )
    println("   aµ[$diag](CC-pQCD) = $(print_uwreal(AMU_pQCD[diag]))")
    println("-------------------------------------------------")
end

for (k,diag) in enumerate(["NLOa","NLOb","NLOa&b"])
    ax = [ax1,ax2,ax3][k]

    ax.errorbar(factor*AMU[diag]["CC"].mean, xerr=factor*AMU[diag]["CC"].err, 0.0, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2)
    ax.errorbar(factor*AMU[diag]["CC"].mean, xerr=factor*AMUERR[diag]["CC"], 0.0, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2)
    ax.errorbar(AMU_pQCD[diag].mean, xerr=AMU_pQCD[diag].err, -1.0, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2)

    ax.set_ylim(-1.5,0.5)

    if k==1
        ax.set_yticks([0.0,-1.0],["CLS","pQCD"])
    else
        ax.set_yticks([])  # Remove redundant y-ticks on the second plot
    end
    diag_str = diag != "NLOa&b" ? diag : "NLOa\\&b"
    ax.set_xlabel(latexstring("(4/9)\\ a_\\mu^{c,c}[\\mathrm{$(diag_str)}]\\ \\times\\ 10^{-11}"))
end

tight_layout()
display(gcf())

close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

# Contribution tales

@info("Table for SDsub window partial results")
println("
\\begin{table}[H]
    \\centering
    \\small
    \\begin{tabular}{c|c c}
        \$(i)\$ & \$(a_\\mu^{3,3}[(i)])^{\\mathrm{SD}}_{\\mathrm{sub}}(Q^{(3,3)})\$ & \$\\frac{4}{9}(a_\\mu^{c,c}[(i)])^{\\mathrm{SD}}_{\\mathrm{sub}}(Q^{(c,c)})\$ \\\\
        \\hline")
for diag in ["NLOa","NLOb","NLOa&b"]
    i_str = diag != "NLOa&b" ? "(4$(diag[4:end]))   " : "(4a\\&b)"
    println("        $(i_str) & $(print_uwreal(10*amu[diag]["SDsub"]["33"],10*amusyst[diag]["SDsub"]["33"])) & $(print_uwreal(10*(4/9)*amu[diag]["SDsub"]["CC"],10*(4/9)*amusyst[diag]["SDsub"]["CC"])) \\\\")
end
println("    \\end{tabular}
    \\break\\break\\break
    \\begin{tabular}{c|c c c}
        \$(i)\$ & \$\\frac{1}{3}\\Delta_{ls}(a_\\mu[(i)])^{\\mathrm{SD}}\$ & \$\\frac{1}{3}\\Delta_{ls}^{\\mathrm{conn}}(a_\\mu[(i)])^{\\mathrm{SD}}\$ & \$\\frac{4}{9}\\Delta_{lc}(b[(i)])(Q^{(c,c)})\$  \\\\
        \\hline")
for diag in ["NLOa","NLOb","NLOa&b"]
    i_str = diag != "NLOa&b" ? "(4$(diag[4:end]))   " : "(4a\\&b)"
    println("        $(i_str) & $(print_uwreal(10/3*amu[diag]["SDsub"]["∆ls"],10/3*amusyst[diag]["SDsub"]["∆ls"])) & $(print_uwreal(10/3*amu[diag]["SDsub"]["∆ls_conn"],10/3*amusyst[diag]["SDsub"]["∆ls_conn"])) & $(print_uwreal(10*(4/9)*amu[diag]["SDsub"]["∆lc"],10*(4/9)*amusyst[diag]["SDsub"]["∆lc"])) \\\\")
end 
println("    \\end{tabular}
    \\caption{Partial results for the SD subtracted quantities. All results are given in units of \$10^{-11}\$. And \$Q^{(3,3)}=5\\, \\mathrm{GeV}\$ and \$Q^{(c,c)}=4\\, \\mathrm{GeV}\$ has been chosen.}
    \\label{tab:isoQCD_SDpre}
\\end{table}
")

#

for wind in ["SD","ID","LD"]
    @info("Table for $wind window results")
    println("
\\begin{table}[H]
    \\centering
    \\small
    \\begin{tabular}{c|c c c|c c}
        \$(i)\$ & \$(a_\\mu^{3,3}[(i)])^{\\mathrm{$(wind)}}\$ & \$\\frac{1}{3}(a_\\mu^{8,8}[(i)])^{\\mathrm{$(wind)}}\$ & \$\\frac{4}{9}(a_\\mu^{c,c}[(i)])^{\\mathrm{$(wind)}}\$ & \$\\frac{1}{9}(a_\\mu^{s,s}[(i)])^{\\mathrm{$(wind)}}\$ & \$(a_\\mu^{\\mathrm{disc}}[(i)])^{\\mathrm{$(wind)}}\$ \\\\
        \\hline")
    for diag in ["NLOa","NLOb","NLOa&b"]
        i_str = diag != "NLOa&b" ? "(4$(diag[4:end]))   " : "(4a\\&b)"
        println("        $(i_str) & $(print_uwreal(10*amu[diag][wind]["33"],10*amusyst[diag][wind]["33"])) & $(print_uwreal(10/3*amu[diag][wind]["88"],10/3*amusyst[diag][wind]["88"])) & $(print_uwreal(10*(4/9)*amu[diag][wind]["CC"],10*(4/9)*amusyst[diag][wind]["CC"]))  & $(print_uwreal(10*(1/9)*amu[diag][wind]["SS"],10*(1/9)*amusyst[diag][wind]["SS"]))  & $(print_uwreal(10*amu[diag][wind]["disc"],10*amusyst[diag][wind]["disc"])) \\\\")
    end
    if wind == "SD"
        cap = "Isopspin and flavor decomposition of the SD window for diagrams NLOa, NLOb and their combination. The light contribution can be obtained from the isovector channel, \$(5/9)\\, a_\\mu^{l,l}[(i)] = 10/9\\, a_\\mu^{3,3}[(i)]\$.  All results are given in units of \$10^{-11}\$."
    else
        cap = "Same as Tab.~\\ref{tab:isoQCD_SD} but for the $wind window."
    end
    println("    \\end{tabular}
    \\caption{$cap}
    \\label{tab:isoQCD_$(wind)}
\\end{table}
")
end

#

@info("isoQCD final results")
println("
\\begin{table}[H]
    \\centering
    \\small
    \\begin{tabular}{c|c c c}
        \$(i)\$ & \$(a_\\mu^{\\mathrm{hvp}}[(i)])^{\\mathrm{SD}}\$ & \$(a_\\mu^{\\mathrm{hvp}}[(i)])^{\\mathrm{ID}}\$ & \$(a_\\mu^{\\mathrm{hvp}}[(i)])^{\\mathrm{LD}}\$ \\\\
        \\hline")
for diag in ["NLOa","NLOb","NLOa&b"]
    i_str = diag != "NLOa&b" ? "(4$(diag[4:end]))   " : "(4a\\&b)"
    println("        $(i_str) & $(print_uwreal(10*AMU[diag]["SD"],10*AMUSYST[diag]["SD"])) & $(print_uwreal(10*AMU[diag]["ID"],10*AMUSYST[diag]["ID"])) & $(print_uwreal(10*AMU[diag]["LD"],10*AMUSYST[diag]["LD"])) \\\\")
end
println("    \\end{tabular}
    \\break\\break\\break
    \\begin{tabular}{c|c c c|c c}
        \$(i)\$ & \$a_\\mu^{3,3}[(i)]\$ & \$\\frac{1}{3}a_\\mu^{8,8}[(i)]\$ & \$\\frac{4}{9}a_\\mu^{c,c}[(i)]\$ & \$\\frac{1}{9}a_\\mu^{s,s}[(i)]\$ & \$a_\\mu^{\\mathrm{disc}}[(i)]\$  \\\\
        \\hline")
for diag in ["NLOa","NLOb","NLOa&b"]
    i_str = diag != "NLOa&b" ? "(4$(diag[4:end]))   " : "(4a\\&b)"
    println("        $(i_str) & $(print_uwreal(10*AMU[diag]["33"],10*AMUSYST[diag]["33"])) & $(print_uwreal(10/3*AMU[diag]["88"],10/3*AMUSYST[diag]["88"])) & $(print_uwreal(10*(4/9)*AMU[diag]["CC"],10*(4/9)*AMUSYST[diag]["CC"])) & $(print_uwreal(10*(1/9)*AMU[diag]["SS"],10*(1/9)*AMUSYST[diag]["SS"]))  & $(print_uwreal(10*AMU[diag]["disc"],10*AMUSYST[diag]["disc"])) \\\\")
end 
println("    \\end{tabular}
    \\caption{Final estimation in isoQCD for diagrams NLOa, NLOb, and their combination. We break down our results in the subsequent time-windows, and the isospin and flavor decomposition.}
    \\label{tab:isoQCD_result}
\\end{table}
")

#

@info("Table for NLOc results")
println("
\\begin{table}[H]
    \\centering
    \\small
    \\begin{tabular}{c c c c c c}
        \$a_\\mu^{3,3-3,3}\$ & \$\\frac{2}{3}a_\\mu^{3,3-8,8}\$ & \$\\frac{8}{9}a_\\mu^{3,3-c,c}\$ & \$\\frac{1}{9}a_\\mu^{8,8-8,8}\$ & \$\\frac{8}{27}a_\\mu^{8,8-c,c}\$ & \$\\frac{16}{91}a_\\mu^{c,c-c,c}\$ \\\\
        \\hline
        $(print_uwreal(10*AMU["NLOc"]["3333"],10*AMUSYST["NLOc"]["3333"])) & $(print_uwreal(10*(2/3)*AMU["NLOc"]["3388"],10*(2/3)*AMUSYST["NLOc"]["3388"])) & $(print_uwreal(10*(8/9)*AMU["NLOc"]["33CC"],10*(8/9)*AMUSYST["NLOc"]["33CC"])) & $(print_uwreal(10*(1/9)*AMU["NLOc"]["8888"],10*(1/9)*AMUSYST["NLOc"]["8888"])) & $(print_uwreal(10*(8/27)*AMU["NLOc"]["88CC"],10*(8/27)*AMUSYST["NLOc"]["88CC"])) & $(print_uwreal(10*(16/91)*AMU["NLOc"]["CCCC"],10*(16/91)*AMUSYST["NLOc"]["CCCC"])) \\\\
    \\end{tabular}
    \\caption{IsoQCD estimation for all NLOc pieces. All values are given in units of \$10^{-11}\$.}
    \\label{tab:isoQCD_NLOc}
\\end{table}
")

##

@info("Table for charm-diisconnected results")
println("
\\begin{table}[H]
    \\centering
    \\small
    \\begin{tabular}{c c c}
        \$(i)\$ & \$\\frac{4}{9}a_\\mu^{c,c(\\mathrm{disc})}[(i)]\$ & \$\\frac{2}{3\\sqrt{3}}a_\\mu^{c,8(\\mathrm{disc})}[(i)]\$ \\\\
        \\hline")
for diag in ["NLOa","NLOb","NLOa&b"]
    i_str = diag != "NLOa&b" ? "(4$(diag[4:end]))   " : "(4a\\&b)"
    println("        $(i_str) & \$$(print_uwreal(1e-10*(4/9)*AMU[diag]["CCdisc"],1e-10*(4/9)*AMUSYST[diag]["CCdisc"],LaTeX=true))\$ & \$$(print_uwreal(1e-10*(2/(3*sqrt(3)))*AMU[diag]["C8disc"],1e-10*(2/(3*sqrt(3)))*AMUSYST[diag]["C8disc"],LaTeX=true))\$ \\\\")
end
println("    \\end{tabular}
    \\caption{Charm quark disconnected contributions to the main diagrams.}
    \\label{tab:charm_disc}
\\end{table}
")


##

# @info("Table for NLOc results")
# println("
# \\begin{table}[H]
#     \\centering
#     \\small
#     \\begin{tabular}{c c c c c c}
#         \$a_\\mu^{3,3-3,3}\$ & \$\\frac{2}{3}a_\\mu^{3,3-8,8}\$ & \$\\frac{8}{9}a_\\mu^{3,3-c,c}\$ & \$\\frac{1}{9}a_\\mu^{8,8-8,8}\$ & \$\\frac{8}{27}a_\\mu^{8,8-c,c}\$ & \$\\frac{16}{91}a_\\mu^{c,c-c,c}\$ \\\\
#         \\hline
#         $(print_uwreal(1000*AMU["NLOc"]["3333"],1000*AMUSYST["NLOc"]["3333"])) & $(print_uwreal(1000*(2/3)*AMU["NLOc"]["3388"],1000*(2/3)*AMUSYST["NLOc"]["3388"])) & $(print_uwreal(1000*(8/9)*AMU["NLOc"]["33CC"],1000*(8/9)*AMUSYST["NLOc"]["33CC"])) & $(print_uwreal(1000*(1/9)*AMU["NLOc"]["8888"],1000*(1/9)*AMUSYST["NLOc"]["8888"])) & $(print_uwreal(1000*(8/27)*AMU["NLOc"]["88CC"],1000*(8/27)*AMUSYST["NLOc"]["88CC"])) & $(print_uwreal(1000*(16/91)*AMU["NLOc"]["CCCC"],1000*(16/91)*AMUSYST["NLOc"]["CCCC"])) \\\\
#     \\end{tabular}
#     \\caption{IsoQCD estimation for all NLOc pieces. All values are given in units of \$10^{-13}\$.}
#     \\label{tab:isoQCD_NLOc}
# \\end{table}
# ")

