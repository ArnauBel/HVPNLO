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

# include("HVPtool/uwConst.jl")

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
    "33" => 1., "88" => 1/3., "CC" => 4/9., "BB" => 1/9., "∆ls_amu" => 1/3., "∆lc_b" => 4/9.,
    "33s" => "", "88s" => "(1/3)", "CCs" => "(4/9)", "BBs" => "(1/9)", "∆ls_amus" => "(1/3)", "∆lc_bs" => "(4/9)",
    "3333" => 1., "3388" => 2/3., "33CC" => 8/9., "8888" => 1/9., "88CC" => 8/27., "CCCC" => 16/81., 
    "3333s" => "", "3388s" => "(2/3)", "33CCs" => "(8/9)", "8888s" => "(1/9)", "88CCs" => "(8/27)", "CCCCs" => "(16/81)",
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

BL_factor = 1.43

BLIND_LD  = Any[true,BL_factor]

path_bdio = path_bdio_dict["local"]


Q33 = 5.0
QCC = 5.0

STD_DERIV  = false
tl_IMPR    = true
VREF       = true
RESC       = false


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
    "NLOa"   => -0.210485,
    "NLOb"   => 0.0432925,
    "NLOa&b" => -0.167192
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

    amu[diag]      = Dict("SD" => Dict(), "ID" => Dict(), "LD" => Dict())
    amusyst[diag]  = Dict("SD" => Dict(), "ID" => Dict(), "LD" => Dict())
    AMU[diag]      = Dict()
    AMUSYST[diag]  = Dict()
    AMUt0ERR[diag] = Dict()
    AMUERR[diag]   = Dict()
    if VREF
        AMUFVC[diag] = Dict()
    end

    amu33sub, info33sub = BDIOread_MAtot(path_bdio,diag,"SDsub","g33",StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,Q=Q33)
    b33Pert = TXTread_bQ(path_bPert,diag)[Qlist .== Q33][1]; uwerr(b33Pert)
    if VREF
        FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,"SDsub",Q=Q33)
        amu33sub += FVC_ChPT
        AMUFVC[diag]["SD"] = abs(0.1*FVC_ChPT)
    end
    amu[diag]["SD"]["33"]     = amu33sub + Window("SD")(0) * b33Pert.mean
    amusyst[diag]["SD"]["33"] = sqrt(info33sub["syst"]^2 + Window("SD")(0)^2 * b33Pert.err^2)

    ∆ls_amu, info∆ls_amu = BDIOread_MAtot(path_bdio,diag,"SDsub","∆ls_amu",StdDer=STD_DERIV)

    amu[diag]["SD"]["88"]     = amu[diag]["SD"]["33"] + ∆ls_amu
    amusyst[diag]["SD"]["88"] = amusyst[diag]["SD"]["33"] + info∆ls_amu["syst"]

    amuCCsub, infoCCsub = BDIOread_MAtot(path_bdio,diag,"SDsub","gCCconn",StdDer=STD_DERIV,Q=QCC)
    ∆lc_b, info∆lc_b    = BDIOread_MAtot(path_bdio,diag,"SDsub","∆lc_b",StdDer=STD_DERIV,Q=QCC)
    amu[diag]["SD"]["CC"] = amuCCsub + Window("SD")(0) * (2*b33Pert.mean + ∆lc_b); uwerr(amu[diag]["SD"]["CC"])
    der_mDs = mchist(amu[diag]["SD"]["CC"], "MD_ph [GeV]")[1] / artificial_err
    amu[diag]["SD"]["CC"]    += value(MD_ph - MD_ph_prime) * der_mDs
    amusyst[diag]["SD"]["CC"] = sqrt(infoCCsub["syst"]^2 + Window("SD")(0)^2 * (2*err(b33Pert)^2 + info∆lc_b["syst"]^2))

    AMU[diag]["SD"]      = amu[diag]["SD"]["33"] + (1/3) * amu[diag]["SD"]["88"] + (4/9) * amu[diag]["SD"]["CC"] + (1/9) * amu_bb[diag]
    AMUSYST[diag]["SD"]  = sqrt(amusyst[diag]["SD"]["33"]^2 + 1/9 * amusyst[diag]["SD"]["88"]^2 + 16/81 * amusyst[diag]["SD"]["CC"]^2)
    AMUt0ERR[diag]["SD"] = get_t0err([AMU[diag]["SD"]],sqrtt0_ph_Madrid)[1]
    AMUERR[diag]["SD"]   = !VREF ? sqrt(AMU[diag]["SD"].err^2 + AMUSYST[diag]["SD"]^2 + AMUt0ERR[diag]["SD"]^2) : sqrt(AMU[diag]["SD"].err^2 + AMUSYST[diag]["SD"]^2 + AMUt0ERR[diag]["SD"]^2 + AMUFVC[diag]["SD"]^2)

    # Add botttom effects on it (only considered to affedt the SD piece)
    amu[diag]["SD"]["BB"] = AMU[diag]["BB"] = uwreal(amu_bb[diag])

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
        amu[diag][wind]["33"]     = amu33
        amusyst[diag][wind]["33"] = info33["syst"]

        amu88, info88 = BDIOread_MAtot(path_bdio,diag,wind,"g88",StdDer=STD_DERIV,BLIND=BLIND)
        amu[diag][wind]["88"]     = amu88
        amusyst[diag][wind]["88"] = info88["syst"]

        amuCC, infoCC = BDIOread_MAtot(path_bdio,diag,wind,"gCCconn",StdDer=STD_DERIV,BLIND=false); uwerr(amuCC)
        der_mDs = mchist(amuCC, "MD_ph [GeV]")[1] / artificial_err
        amu[diag][wind]["CC"]     = amuCC + value(MD_ph - MD_ph_prime) * der_mDs
        amusyst[diag][wind]["CC"] = infoCC["syst"]


        AMU[diag][wind]      = (amu[diag][wind]["33"] + (1/3) * amu[diag][wind]["88"])/BLIND_factor + (4/9) * amu[diag][wind]["CC"]
        AMUSYST[diag][wind]  = sqrt((amusyst[diag][wind]["33"]^2 + 1/9 * amusyst[diag][wind]["88"]^2)/BLIND_factor^2 + 16/81 * amusyst[diag][wind]["CC"]^2)
        AMUt0ERR[diag][wind] = get_t0err([AMU[diag][wind]],sqrtt0_ph_Madrid)[1]
        AMUERR[diag][wind]   = !VREF ? sqrt(AMU[diag][wind].err^2 + AMUSYST[diag][wind]^2 + AMUt0ERR[diag][wind]^2) : sqrt(AMU[diag][wind].err^2 + AMUSYST[diag][wind]^2 + AMUt0ERR[diag][wind]^2 + AMUFVC[diag][wind]^2)
    end

    if VREF
        FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,"NW")
        AMUFVC[diag]["NW"] = abs(0.1*FVC_ChPT)
    end

    for comp in ["33","88","CC"]
        BLIND_factor = BLIND_LD[1] ? BLIND_LD[2] : 1.0

        AMU[diag][comp]      = amu[diag]["SD"][comp] + amu[diag]["ID"][comp] + amu[diag]["LD"][comp]/BLIND_factor
        AMUSYST[diag][comp]  = sqrt((amusyst[diag]["SD"][comp])^2 + (amusyst[diag]["ID"][comp])^2 + (amusyst[diag]["LD"][comp]/BLIND_factor)^2)
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

@info("Results in isoQCD :")
for diag in ["NLOa","NLOb","NLOa&b"]
    # AMU[diag]["tot"] = AMU[diag]["33"] + (1/3)*AMU[diag]["88"] + (4/9)*AMU[diag]["CC"]; print_uwreal(AMU[diag]["tot"])
    AMU[diag]["tot"]      = AMU[diag]["SD"] + AMU[diag]["ID"] + AMU[diag]["LD"]; uwerr(AMU[diag]["tot"])
    AMUSYST[diag]["tot"]  = sqrt(AMUSYST[diag]["SD"]^2 + AMUSYST[diag]["ID"]^2 + AMUSYST[diag]["LD"]^2)
    AMUt0ERR[diag]["tot"] = get_t0err([AMU[diag]["tot"]],sqrtt0_ph_Madrid)[1]
    AMUFVC[diag]["tot"]   = AMUFVC[diag]["NW"]
    AMUERR[diag]["tot"] = sqrt(AMU[diag]["tot"].err^2 + AMUSYST[diag]["tot"]^2 + AMUt0ERR[diag]["tot"]^2 + AMUFVC[diag]["tot"]^2)

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


println("---------------------------------------------------------")

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
            println("---------------------------------------------------------")
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
    "NLOa"   => [[uwreal([-21.61536,0.1161895],"Spacelike"),"Spacelike - Mainz 2025"]],
    "NLOb"   => [[uwreal([11.118475,0.0791056],"Spacelike"),"Spacelike - Mainz 2025"]],
    "NLOc"   => [
        # [uwreal([0.39710,0.00966],"Spacelike"),"Spacelike - Mainz 2022"],
        # [uwreal([0.37685,0.00837],"Spacelike"),"Spacelike - Mainz 2025"],
        [uwreal([0.3735,0.0081],"Spacelike"),"Spacelike - Mainz 2025"],
        ],
    "NLO"    => [[uwreal([-10.13129,0.03680297],"Spacelike"),"Spacelike - Mainz 2025"]],
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
    for res in reverse(RES_sl[diag])
        i -= 1
        uwerr(res[1])
        ax.errorbar(10*(res[1].mean+AMUIB[diag].mean), xerr=10*sqrt(res[1].err^2+AMUIB[diag].err^2), i, 0.0, fmt="d", color="green", ms=10, capsize=2)
        push!(y_ticks,res[2])
    end
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

[uwerr(AMU_pQCD[key]) for key in keys(AMU_pQCD)]

fig, (ax1, ax2) = subplots(1, 2,
    gridspec_kw = Dict("width_ratios" => [3, 3], "wspace" => 0),
    figsize = (12, 6)
)

factor = 10*4/9

ax1.errorbar(factor*AMU["NLOa"]["CC"].mean, xerr=factor*AMU["NLOa"]["CC"].err, 0.0, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2)
ax1.errorbar(factor*AMU["NLOa"]["CC"].mean, xerr=factor*AMUERR["NLOa"]["CC"], 0.0, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2)
ax1.errorbar(AMU_pQCD["NLOa"].mean, xerr=AMU_pQCD["NLOa"].err, -1.0, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2)

ax2.errorbar(factor*AMU["NLOb"]["CC"].mean, xerr=factor*AMU["NLOb"]["CC"].err, 0.0, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2)
ax2.errorbar(factor*AMU["NLOb"]["CC"].mean, xerr=factor*AMUERR["NLOb"]["CC"], 0.0, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2)
ax2.errorbar(AMU_pQCD["NLOb"].mean, xerr=AMU_pQCD["NLOb"].err, -1.0, 0.0, fmt="o", mfc="none", color="black", ms=10, capsize=2)

ax1.set_ylim(-1.5,0.5)
ax2.set_ylim(-1.5,0.5)

ax1.set_yticks([0.0,-1.0], ["Lattice. Mainz","pQCD (Toni)"], rotation = 0, fontsize=15)
ax2.set_yticks([])  # Remove redundant y-ticks on the second plot

ax1.set_xlabel(latexstring("(4/9)\\ a_\\mu^{c,c}[\\mathrm{NLOa}]\\ \\times\\ 10^{-11}"))
ax2.set_xlabel(latexstring("(4/9)\\ a_\\mu^{c,c}[\\mathrm{NLOb}]\\ \\times\\ 10^{-11}"))

tight_layout()
display(gcf())

close()
