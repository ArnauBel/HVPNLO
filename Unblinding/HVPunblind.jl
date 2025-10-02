# Import packages

using Revise

include("../HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs

using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

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

charge_factor = Dict(
    "g33" => 1., "g88" => 1/3., "gCCconn" => 4/9., "∆ls_amu" => 1/3., "∆lc_b" => 4/9.,
    "g33s" => "", "g88s" => "(1/3)", "gCCconns" => "(4/9)", "∆ls_amus" => "(1/3)", "∆lc_bs" => "(4/9)"
)

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

BL_factor = 1.5

BLIND_LD  = Any[true,BL_factor]

Q33 = 5.0
QCC = 5.0

STD_DERIV  = false
tl_IMPR    = true
VREF       = true
RESC       = false

AMU      = Dict()
AMUSYST  = Dict()
AMUT0ERR = Dict()
AMUERR   = Dict()
if VREF
    AMUFVC = Dict()
end

for diag in ["NLOa","NLOb","NLOa&b"]

    # For the SD window

    amu      = Dict()
    amusyst  = Dict()
    amut0err = Dict()
    amuerr   = Dict()
    if VREF
        amufvc = Dict()
    end

    path_bdio = path_bdio_dict["local"]

    amu33sub, info33sub = BDIOread_MAtot(path_bdio,diag,"SDsub","g33",StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,Q=Q33)
    b33Pert = TXTread_bQ(path_bPert,diag)[Qlist .== Q33][1]; uwerr(b33Pert)
    if VREF
        FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,"SDsub",Q=Q33)
        amu33sub += FVC_ChPT
        amu33SDfvc = abs(0.1*FVC_ChPT)
    end
    amu33SD = amu33sub + Window("SD")(0) * value(b33Pert); uwerr(amu33SD)
    amu33SDsyst  = sqrt(info33sub["syst"]^2 + Window("SD")(0)^2 * err(b33Pert)^2)
    # amu33SDt0err = get_t0err([amu33SD],sqrtt0_ph_Regensburg)[1]
    # amu33SDerr   = VREF ? sqrt(err(amu33SD)^2 + amu33SDsyst^2 + amu33SDt0err^2) : sqrt(err(amu33SD)^2 + amu33SDsyst^2 + amu33SDt0err^2 + amu33SDfvc^2)

    ∆ls_amu, info∆ls_amu = BDIOread_MAtot(path_bdio,diag,"SDsub","∆ls_amu",StdDer=STD_DERIV)

    amu88SD = amu33SD + ∆ls_amu; uwerr(amu88SD)
    amu88SDsyst  = amu33SDsyst + info∆ls_amu["syst"]
    # amu88SDt0err = get_t0err([amu88SD],sqrtt0_ph_Regensburg)[1]
    # amu88SDerr   = sqrt(err(amu88SD)^2 + amu88SDsyst^2 + amu88SDt0err^2)

    amuCCsub, infoCCsub = BDIOread_MAtot(path_bdio,diag,"SDsub","gCCconn",StdDer=STD_DERIV,Q=QCC)
    ∆lc_b, info∆lc_b    = BDIOread_MAtot(path_bdio,diag,"SDsub","∆lc_b",StdDer=STD_DERIV,Q=QCC)
    amuCCSD = amuCCsub + Window("SD")(0) * (2*value(b33Pert) + ∆lc_b); uwerr(amuCCSD)
    der_mDs = mchist(amuCCSD, "MD_ph [GeV]")[1] / artificial_err
    MD_ph_prime, mDs_SU3 = BDIOread_mDs(path_bdio)
    amuCCSD += value(MD_ph - MD_ph_prime) * der_mDs; uwerr(amuCCSD)
    amuCCSDsyst  = sqrt(infoCCsub["syst"]^2 + Window("SD")(0)^2 * (2*err(b33Pert)^2 + info∆lc_b["syst"]^2))
    # amuCCSDt0err = get_t0err([amuCCSD],sqrtt0_ph_Regensburg)[1]
    # amuCCSDerr   = sqrt(err(amuCCSD)^2 + amuCCSDsyst^2 + amuCCSDt0err^2)

    amuSD = amu33SD + (1/3) * amu88SD + (4/9) * amuCCSD; uwerr(amuSD)
    amuSDsyst  = sqrt(amu33SDsyst^2 + 1/9 * amu88SDsyst^2 + 16/81 * amuCCSDsyst^2)
    amuSDt0err = get_t0err([amuSD],sqrtt0_ph_Regensburg)[1]
    amuSDerr   = VREF ? sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2) : sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2 + amu33SDfvc^2)


    amu["SD"]      = amuSD
    amusyst["SD"]  = amuSDsyst
    amut0err["SD"] = amuSDt0err
    amuerr["SD"]   = amuSDerr
    if VREF
        amufvc["SD"] = amu33SDfvc
    end

    # For the I&LD windows

    for wind in ["ID","LD"]
        BLIND = wind=="LD" ? BLIND_LD[1] : false

        amu33, info33 = BDIOread_MAtot(path_bdio,diag,wind,"g33",StdDer=STD_DERIV,BLIND=BLIND,Vref=VREF)
        if VREF
            FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,wind)
            amu33 += FVC_ChPT
            amu33fvc = abs(0.1*FVC_ChPT)
        end
        amu33syst  = info33["syst"]
        # amu33t0err = get_t0err([amu33],sqrtt0_ph_Regensburg)[1]
        # amu33err   = VREF ? sqrt(err(amu33)^2 + amu33syst^2 + amu33t0err^2) : sqrt(err(amu33)^2 + amu33syst^2 + amu33t0err^2 + amu33fvc^2)

        amu88, info88 = BDIOread_MAtot(path_bdio,diag,wind,"g88",StdDer=STD_DERIV,BLIND=BLIND)
        amu88syst  = info88["syst"]
        # amu88t0err = get_t0err([amu88],sqrtt0_ph_Regensburg)[1]
        # amu88err   = sqrt(err(amu88)^2 + amu88syst^2 + amu88t0err^2)

        amuCC, infoCC = BDIOread_MAtot(path_bdio,diag,wind,"gCCconn",StdDer=STD_DERIV,BLIND=false); uwerr(amuCC)
        der_mDs = mchist(amuCC, "MD_ph [GeV]")[1] / artificial_err
        MD_ph_prime, mDs_SU3 = BDIOread_mDs(path_bdio)
        amuCC = amuCC + value(MD_ph - MD_ph_prime) * der_mDs; uwerr(amuCC)
        amuCCsyst  = infoCC["syst"]
        # amuCCt0err = get_t0err([amuCC],sqrtt0_ph_Regensburg)[1]
        # amuCCerr   = sqrt(err(amuCC)^2 + amuCCsyst^2 + amuCCt0err^2)

        amuW = amu33 + (1/3) * amu88 + (4/9) * amuCC; uwerr(amuW)
        amuWsyst  = sqrt(amu33syst^2 + 1/9 * amu88syst^2 + 16/81 * amuCCsyst^2)
        amuWt0err = get_t0err([amuW],sqrtt0_ph_Regensburg)[1]
        amuWerr   = VREF ? sqrt(err(amuW)^2 + amuWsyst^2 + amuWt0err^2) : sqrt(err(amuW)^2 + amuWsyst^2 + amuWt0err^2 + amu33fvc^2)

        amu[wind]      = amuW
        amusyst[wind]  = amuWsyst
        amut0err[wind] = amuWt0err
        amuerr[wind]   = amuWerr
        if VREF
            amufvc[wind] = amu33fvc
        end
    end


    BLIND_LD_factor = BLIND_LD[1] ? BLIND_LD[2] : 1.0

    AMU[diag]      = amu["SD"] + amu["ID"] + amu["LD"]/BLIND_LD_factor 
    AMUSYST[diag]  = sqrt(amusyst["SD"]^2 + amusyst["ID"]^2 + (amusyst["LD"]/BLIND_LD_factor)^2)
    AMUT0ERR[diag] = get_t0err([AMU[diag]],sqrtt0_ph_Regensburg)[1]
    AMUERR[diag]   = sqrt(AMU[diag].err^2 + AMUSYST[diag]^2 + AMUT0ERR[diag]^2)
    if VREF
        AMUFVC[diag] = abs(0.1*JDL2read_FVC_ChPT(path_FVCcont,diag,"NW"))
    end
end

