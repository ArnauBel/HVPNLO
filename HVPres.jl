# Import packages

using Revise

include("HVPtool/HVPtool.jl")
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
    "local" => joinpath(julia_script_directory, "..", "ObsBDIO"),
    # "local" => joinpath(julia_script_directory, "..", "ObsCrosschecks", "Obs4LD(t0=Regensburg)", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "ObsExternal", "PortableSSD", "ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_bPert   = joinpath(julia_script_directory, "..", "PertSD")
path_FVCcont = joinpath(julia_script_directory, "..", "FVCcont")

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

##==========================> INDIVIDUAL RESULTS <==========================##


diag = "LO"  # LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
wind = "LD"  # NW  SD  ID  LD  ILD
comp = "g33"  # g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b

IMPR = "all"  # "all"  "1"  "1old"  "2"

Q = 5.0  # virtuality for SDsub

BLIND = false

STD_DERIV  = false
tl_IMPR    = false
VREF       = true
RESC       = true

if wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]
    @info(" RES [diag. $diag; wind. $wind; comp. $comp] \n - Vref  : $VREF  \n - Rescal: $RESC \n - Q: $Q")
else
    @info(" RES [diag. $diag; wind. $wind; comp. $comp] \n - Vref  : $VREF  \n - Rescal: $RESC")
end

path_bdio = path_bdio_dict["local"]

if IMPR == "all"
    AMU, INFO = BDIOread_MAtot(path_bdio,diag,wind,comp,resc=RESC,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q); uwerr(AMU)
    SYSTerr  = INFO["syst"]
else
    amu, info = BDIOread_MAtot(path_bdio,diag,wind,comp,resc=RESC,read="impr",impr_set=IMPR,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)
    
    weight = info["weight"]
    res_tot = amu["res_tot"]
    
    mykeys = DictComptoKey[comp]
    WEIGHT  = vcat([weight[key] for key in mykeys]...)./(length(mykeys))
    RES_TOT = vcat([res_tot[key] for key in mykeys]...)

    AMU, SYSTerr = model_average(RES_TOT, WEIGHT); AMU = AMU[1]
end

scale_ph = !RESC ? sqrtt0_ph_Regensburg :  fPi_ph_PDGFLAG
SCALEerr = get_t0err([AMU],scale_ph,resc=RESC)[1]

if VREF
    FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,wind)
    AMU += FVC_ChPT
end
uwerr(AMU)
dig = comp in ["gCCdisc","gC8disc"] ? 7 : 5

factor = diag == "LO" ? charge_factor[comp] : charge_factor[comp]*10
if !VREF
    println("     ⟹ $(charge_factor[comp*"s"]) amu[$diag($wind)|$comp] = $(round((factor)*value(AMU),digits=dig))($(round((factor)*err(AMU),digits=dig)))($(round((factor)*SYSTerr,digits=dig)))($(round((factor)*SCALEerr,digits=dig)))[$(round((factor)*sqrt(err(AMU)^2+SYSTerr^2+SCALEerr^2),digits=dig))]")
else
    println("     ⟹ $(charge_factor[comp*"s"]) amu[$diag($wind)|$comp] = $(round((factor)*value(AMU),digits=dig))($(round((factor)*err(AMU),digits=dig)))($(round((factor)*SYSTerr,digits=dig)))($(round((factor)*SCALEerr,digits=dig)))($(round((factor)*0.1*abs(FVC_ChPT),digits=dig)))[$(round((factor)*sqrt(err(AMU)^2+SYSTerr^2+SCALEerr^2+(0.1*FVC_ChPT)^2),digits=dig))]")
end

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> SD RESULTS <==========================##

diag = "NLOa&b"  #  LO  NLOa  NLOb  NLOa&b

Q33 = 5.0
QCC = 5.0

VREF      = false
tl_IMPR   = true
STD_DERIV = false

path_bdio = path_bdio_dict["local"]

amu33sub, info33sub = BDIOread_MAtot(path_bdio,diag,"SDsub","g33",StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,Q=Q33)
b33Pert = TXTread_bQ(path_bPert,diag)[Qlist .== Q33][1]; uwerr(b33Pert)

amu33SD = amu33sub + Window("SD")(0) * value(b33Pert); uwerr(amu33SD)
amu33SDsyst  = sqrt(info33sub["syst"]^2 + Window("SD")(0)^2 * err(b33Pert)^2)
amu33SDt0err = get_t0err([amu33SD],sqrtt0_ph_Regensburg)[1]
amu33SDerr   = sqrt(err(amu33SD)^2 + amu33SDsyst^2 + amu33SDt0err^2)

∆ls_amu, info∆ls_amu = BDIOread_MAtot(path_bdio,diag,"SDsub","∆ls_amu",StdDer=STD_DERIV,Vref=VREF)

amu88SD = amu33SD + ∆ls_amu; uwerr(amu88SD)
amu88SDsyst  = amu33SDsyst + info∆ls_amu["syst"]
amu88SDt0err = get_t0err([amu88SD],sqrtt0_ph_Regensburg)[1]
amu88SDerr   = sqrt(err(amu88SD)^2 + amu88SDsyst^2 + amu88SDt0err^2)

amuCCsub, infoCCsub = BDIOread_MAtot(path_bdio,diag,"SDsub","gCCconn",StdDer=STD_DERIV,Vref=VREF,Q=QCC)
∆lc_b, info∆lc_b    = BDIOread_MAtot(path_bdio,diag,"SDsub","∆lc_b",StdDer=STD_DERIV,Vref=VREF,Q=QCC)
amuCCSD = amuCCsub + Window("SD")(0) * (2*value(b33Pert) + ∆lc_b); uwerr(amuCCSD)
der_mDs = mchist(amuCCSD, "MD_ph [GeV]")[1] / artificial_err
MD_ph_prime = BDIOread_mDs(path_bdio)["mDs ph"]
amuCCSD += value(MD_ph - MD_ph_prime) * der_mDs; uwerr(amuCCSD)
amuCCSDsyst  = sqrt(infoCCsub["syst"]^2 + Window("SD")(0)^2 * (2*err(b33Pert)^2 + info∆lc_b["syst"]^2))
amuCCSDt0err = get_t0err([amuCCSD],sqrtt0_ph_Regensburg)[1]
amuCCSDerr   = sqrt(err(amuCCSD)^2 + amuCCSDsyst^2 + amuCCSDt0err^2)

amuSD = amu33SD + (1/3) * amu88SD + (4/9) * amuCCSD; uwerr(amuSD)
amuSDsyst  = sqrt(amu33SDsyst^2 + 1/9 * amu88SDsyst^2 + 16/81 * amuCCSDsyst^2)
amuSDt0err = get_t0err([amuSD],sqrtt0_ph_Regensburg)[1]
amuSDerr   = sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2)

digits = 4
factor = diag == "LO" ? 1 : 10

println(" Iso-spin decomposition:")
println("---------------------------------------------------------")
println(" => amu($diag;SD;33) = $(round(value(amu33SD)*factor,digits=digits))($(round(err(amu33SD)*factor,digits=digits)))($(round(amu33SDsyst*factor,digits=digits)))($(round(amu33SDt0err*factor,digits=digits)))[$(round(amu33SDerr*factor,digits=digits))]")
println(" => (1/3)amu($diag;SD;88) = $(round(value(amu88SD)*factor/3,digits=digits))($(round(err(amu88SD)*factor/3,digits=digits)))($(round(amu88SDsyst*factor/3,digits=digits)))($(round(amu88SDt0err*factor/3,digits=digits)))[$(round(amu88SDerr*factor/3,digits=digits))]")
println(" => (4/9)amu($diag;SD;CC) = $(round(value(amuCCSD)*factor*4/9,digits=digits))($(round(err(amuCCSD)*factor*4/9,digits=digits)))($(round(amuCCSDsyst*factor*4/9,digits=digits)))($(round(amuCCSDt0err*factor*4/9,digits=digits)))[$(round(amuCCSDerr*factor*4/9,digits=digits))] \n")

println(" Final estimation for diag $diag SD window")
println("---------------------------------------------------------")
println(" => amu($diag;SD) = $(round(value(amuSD)*factor,digits=digits))($(round(err(amuSD)*factor,digits=digits)))($(round(amuSDsyst*factor,digits=digits)))($(round(amuSDt0err*factor,digits=digits)))[$(round(amuSDerr*factor,digits=digits))]")

##==========================> ID & LD RESULTS <==========================##

diag = "NLOa&b"  #  LO  NLOa  NLOb  NLOa&b
wind = "LD"  #  ID  LD  ILD

BLIND = true

VREF      = false
STD_DERIV = false

path_bdio = path_bdio_dict["local"]

amu33, info33 = BDIOread_MAtot(path_bdio,diag,wind,"g33",StdDer=STD_DERIV,BLIND=BLIND,Vref=VREF)
amu33syst  = info33["syst"]
amu33t0err = get_t0err([amu33],sqrtt0_ph_Regensburg)[1]
amu33err   = sqrt(err(amu33)^2 + amu33syst^2 + amu33t0err^2)

amu88, info88 = BDIOread_MAtot(path_bdio,diag,wind,"g88",StdDer=STD_DERIV,BLIND=BLIND)
amu88syst  = info88["syst"]
amu88t0err = get_t0err([amu88],sqrtt0_ph_Regensburg)[1]
amu88err   = sqrt(err(amu88)^2 + amu88syst^2 + amu88t0err^2)

amuCC, infoCC = BDIOread_MAtot(path_bdio,diag,wind,"gCCconn",StdDer=STD_DERIV,BLIND=false); uwerr(amuCC)
der_mDs = mchist(amuCC, "MD_ph [GeV]")[1] / artificial_err
MD_ph_prime = BDIOread_mDs(path_bdio)["mDs ph"]
amuCC = amuCC + value(MD_ph - MD_ph_prime) * der_mDs; uwerr(amuCC)
amuCCsyst  = infoCC["syst"]
amuCCt0err = get_t0err([amuCC],sqrtt0_ph_Regensburg)[1]
amuCCerr   = sqrt(err(amuCC)^2 + amuCCsyst^2 + amuCCt0err^2)

amu = amu33 + (1/3) * amu88 + (4/9) * amuCC; uwerr(amu)
amusyst  = sqrt(amu33syst^2 + 1/9 * amu88syst^2 + 16/81 * amuCCsyst^2)
amut0err = get_t0err([amu],sqrtt0_ph_Regensburg)[1]
amuerr   = sqrt(err(amu)^2 + amusyst^2 + amut0err^2)

digits = 6
factor = diag == "LO" ? 1 : 10

println(" Iso-spin decomposition:")
println("---------------------------------------------------------")
println(" => amu($diag;$wind;33) = $(round(value(amu33)*factor,digits=digits))($(round(err(amu33)*factor,digits=digits)))($(round(amu33syst*factor,digits=digits)))($(round(amu33t0err*factor,digits=digits)))[$(round(amu33err*factor,digits=digits))]")
println(" => (1/3)amu($diag;$wind;88) = $(round(value(amu88)*factor/3,digits=digits))($(round(err(amu88)*factor/3,digits=digits)))($(round(amu88syst*factor/3,digits=digits)))($(round(amu88t0err*factor/3,digits=digits)))[$(round(amu88err*factor/3,digits=digits))]")
println(" => (4/9)amu($diag;$wind;CC) = $(round(value(amuCC)*factor*4/9,digits=digits))($(round(err(amuCC)*factor*4/9,digits=digits)))($(round(amuCCsyst*factor*4/9,digits=digits)))($(round(amuCCt0err*factor*4/9,digits=digits)))[$(round(amuCCerr*factor*4/9,digits=digits))] \n")

println(" Final estimation for diag $diag $wind window")
println("---------------------------------------------------------")
println(" => amu($diag;$wind) = $(round(value(amu)*factor,digits=digits))($(round(err(amu)*factor,digits=digits)))($(round(amusyst*factor,digits=digits)))($(round(amut0err*factor,digits=digits)))[$(round(amuerr*factor,digits=digits))]")


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> COMP RESULTS <==========================##

diag = "NLOb"  #  LO  NLOa  NLOb  NLOa&b
comp = "gCCconn"  #  g33  g88  gCCconn

Q33 = 5.0
QCC = 5.0

VREF      = false
tl_IMPR   = true
STD_DERIV = false


BLIND = Any[false,1.]  #  [false,1.0]  [true,1.5]

path_bdio = path_bdio_dict["local"]

# SD

b33Pert = TXTread_bQ(path_bPert,diag)[Qlist .== Q33][1]; uwerr(b33Pert)

if comp == "g33"
    amuSDsub, infoSDsub = BDIOread_MAtot(path_bdio,diag,"SDsub",comp,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,Q=Q33)
    amuSD = amuSDsub + Window("SD")(0) * value(b33Pert); uwerr(amuSD)
    amuSDsyst  = sqrt(infoSDsub["syst"]^2 + Window("SD")(0)^2 * err(b33Pert)^2)
    amuSDt0err = get_t0err([amuSD],sqrtt0_ph_Regensburg)[1]
    amuSDerr   = sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2)
elseif comp == "88"
    ∆ls_amu, info∆ls_amu = BDIOread_MAtot(path_bdio,diag,"SDsub","∆ls_amu",StdDer=STD_DERIV)
    amu33sub, info33sub = BDIOread_MAtot(path_bdio,diag,"SDsub","g33",StdDer=STD_DERIV,tlImpr=tl_IMPR,Q=Q33)
    amu33 = amu33sub + Window("SD")(0) * value(b33Pert)

    amuSD = amu33 + ∆ls_amu; uwerr(amuSD)
    amuSDsyst  = info33["syst"] + info∆ls_amu["syst"]
    amuSDt0err = get_t0err([amuSD],sqrtt0_ph_Regensburg)[1]
    amuSDerr   = sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2)
elseif comp == "gCCconn"
    amuSDsub, infoSDsub = BDIOread_MAtot(path_bdio,diag,"SDsub",comp,StdDer=STD_DERIV,tlImpr=tl_IMPR,Q=QCC)
    ∆lc_b, info∆lc_b    = BDIOread_MAtot(path_bdio,diag,"SDsub","∆lc_b",StdDer=STD_DERIV,Vref=VREF,Q=QCC)
    amuSD = amuSDsub + Window("SD")(0) * (2*value(b33Pert) + ∆lc_b); uwerr(amuSD)
    amuSDsyst  = sqrt(infoSDsub["syst"]^2 + Window("SD")(0)^2 * (2*err(b33Pert)^2 + info∆lc_b["syst"]^2))
    amuSDt0err = get_t0err([amuSD],sqrtt0_ph_Regensburg)[1]
    amuSDerr   = sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2)
end

# ID

amuID, infoID = BDIOread_MAtot(path_bdio,diag,"ID",comp,StdDer=STD_DERIV,Vref=VREF); uwerr(amuID)
amuIDsyst  = infoID["syst"]
amuIDt0err = get_t0err([amuID],sqrtt0_ph_Regensburg)[1]
amuIDerr   = sqrt(err(amuID)^2 + amuIDsyst^2 + amuIDt0err^2)

# LD

amuLD, infoLD = BDIOread_MAtot(path_bdio,diag,"LD",comp,StdDer=STD_DERIV,Vref=VREF,BLIND=BLIND[1])
amuLD /= BLIND[2]; uwerr(amuLD)
amuLDsyst  = infoLD["syst"] / BLIND[2]
amuLDt0err = get_t0err([amuLD],sqrtt0_ph_Regensburg)[1]
amuLDerr   = sqrt(err(amuLD)^2 + amuLDsyst^2 + amuLDt0err^2)


amu = amuSD + amuID + amuLD; uwerr(amu)
amusyst = sqrt(amuSDsyst^2+amuIDsyst^2+amuLDsyst^2)
amut0err = get_t0err([amu],sqrtt0_ph_Regensburg)[1]
amuerr   = sqrt(err(amu)^2 + amusyst^2 + amut0err^2)

if comp == "gCCconn"
    MD_ph_prime = BDIOread_mDs(path_bdio)["mDs ph"]

    der_mDs    = mchist(amu  , "MD_ph [GeV]")[1] / artificial_err
    der_mDs_SD = mchist(amuSD, "MD_ph [GeV]")[1] / artificial_err
    der_mDs_ID = mchist(amuID, "MD_ph [GeV]")[1] / artificial_err
    der_mDs_LD = mchist(amuLD, "MD_ph [GeV]")[1] / artificial_err

    amu   += value(MD_ph - MD_ph_prime) * der_mDs   ; uwerr(amu)
    amuSD += value(MD_ph - MD_ph_prime) * der_mDs_SD; uwerr(amuSD)
    amuID += value(MD_ph - MD_ph_prime) * der_mDs_ID; uwerr(amuID)
    amuLD += value(MD_ph - MD_ph_prime) * der_mDs_LD; uwerr(amuLD)
end

digits = 4
factor = (diag == "LO" ? 1 : 10) * (charge_factor[comp])


println(" Window decomposition:")
println("---------------------------------------------------------")
println(" => $(charge_factor[comp*"s"]) amu($diag;SD;$comp) = $(round(value(amuSD)*factor,digits=digits))($(round(err(amuSD)*factor,digits=digits)))($(round(amuSDsyst*factor,digits=digits)))($(round(amuSDt0err*factor,digits=digits)))[$(round(amuSDerr*factor,digits=digits))]")
println(" => $(charge_factor[comp*"s"]) amu($diag;ID;$comp) = $(round(value(amuID)*factor,digits=digits))($(round(err(amuID)*factor,digits=digits)))($(round(amuIDsyst*factor,digits=digits)))($(round(amuIDt0err*factor,digits=digits)))[$(round(amuIDerr*factor,digits=digits))]")
println(" => $(charge_factor[comp*"s"]) amu($diag;LD;$comp) = $(round(value(amuLD)*factor,digits=digits))($(round(err(amuLD)*factor,digits=digits)))($(round(amuLDsyst*factor,digits=digits)))($(round(amuLDt0err*factor,digits=digits)))[$(round(amuLDerr*factor,digits=digits))]")


println(" Final estimation for diag $diag comp $comp")
println("---------------------------------------------------------")
println(" => $(charge_factor[comp*"s"]) amu($diag;$comp) = $(round(value(amu)*factor,digits=digits))($(round(err(amu)*factor,digits=digits)))($(round(amusyst*factor,digits=digits)))($(round(amut0err*factor,digits=digits)))[$(round(amuerr*factor,digits=digits))]")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> FULL DIAG <==========================##

diag = "NLOa&b"  #  LO  NLOa  NLOb  NLOa&b

Q33 = 5.0
QCC = 5.0

VREF      = false
tl_IMPR   = true
STD_DERIV = false

BLIND_LD  = Any[true,1.5]

# For the SD window

amu      = Dict()
amusyst  = Dict()
amut0err = Dict()
amuerr   = Dict()

path_bdio = path_bdio_dict["local"]

amu33sub, info33sub = BDIOread_MAtot(path_bdio,diag,"SDsub","g33",StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,Q=Q33)
b33Pert = TXTread_bQ(path_bPert,diag)[Qlist .== Q33][1]; uwerr(b33Pert)

amu33SD = amu33sub + Window("SD")(0) * value(b33Pert); uwerr(amu33SD)
amu33SDsyst  = sqrt(info33sub["syst"]^2 + Window("SD")(0)^2 * err(b33Pert)^2)
amu33SDt0err = get_t0err([amu33SD],sqrtt0_ph_Regensburg)[1]
amu33SDerr   = sqrt(err(amu33SD)^2 + amu33SDsyst^2 + amu33SDt0err^2)

∆ls_amu, info∆ls_amu = BDIOread_MAtot(path_bdio,diag,"SDsub","∆ls_amu",StdDer=STD_DERIV,Vref=VREF)

amu88SD = amu33SD + ∆ls_amu; uwerr(amu88SD)
amu88SDsyst  = amu33SDsyst + info∆ls_amu["syst"]
amu88SDt0err = get_t0err([amu88SD],sqrtt0_ph_Regensburg)[1]
amu88SDerr   = sqrt(err(amu88SD)^2 + amu88SDsyst^2 + amu88SDt0err^2)

amuCCsub, infoCCsub = BDIOread_MAtot(path_bdio,diag,"SDsub","gCCconn",StdDer=STD_DERIV,Vref=VREF,Q=QCC)
∆lc_b, info∆lc_b    = BDIOread_MAtot(path_bdio,diag,"SDsub","∆lc_b",StdDer=STD_DERIV,Vref=VREF,Q=QCC)
amuCCSD = amuCCsub + Window("SD")(0) * (2*value(b33Pert) + ∆lc_b); uwerr(amuCCSD)
der_mDs = mchist(amuCCSD, "MD_ph [GeV]")[1] / artificial_err
MD_ph_prime = BDIOread_mDs(path_bdio)["mDs ph"]
amuCCSD += value(MD_ph - MD_ph_prime) * der_mDs; uwerr(amuCCSD)
amuCCSDsyst  = sqrt(infoCCsub["syst"]^2 + Window("SD")(0)^2 * (2*err(b33Pert)^2 + info∆lc_b["syst"]^2))
amuCCSDt0err = get_t0err([amuCCSD],sqrtt0_ph_Regensburg)[1]
amuCCSDerr   = sqrt(err(amuCCSD)^2 + amuCCSDsyst^2 + amuCCSDt0err^2)

amuSD = amu33SD + (1/3) * amu88SD + (4/9) * amuCCSD; uwerr(amuSD)
amuSDsyst  = sqrt(amu33SDsyst^2 + 1/9 * amu88SDsyst^2 + 16/81 * amuCCSDsyst^2)
amuSDt0err = get_t0err([amuSD],sqrtt0_ph_Regensburg)[1]
amuSDerr   = sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2)


amu["SD"]      = amuSD
amusyst["SD"]  = amuSDsyst
amut0err["SD"] = amuSDt0err
amuerr["SD"]   = amuSDerr

# For the I&LD windows

for wind in ["ID","LD"]
    BLIND = wind=="LD" ? BLIND_LD[1] : false

    amu33, info33 = BDIOread_MAtot(path_bdio,diag,wind,"g33",StdDer=STD_DERIV,BLIND=BLIND,Vref=VREF)
    amu33syst  = info33["syst"]
    amu33t0err = get_t0err([amu33],sqrtt0_ph_Regensburg)[1]
    amu33err   = sqrt(err(amu33)^2 + amu33syst^2 + amu33t0err^2)

    amu88, info88 = BDIOread_MAtot(path_bdio,diag,wind,"g88",StdDer=STD_DERIV,BLIND=BLIND)
    amu88syst  = info88["syst"]
    amu88t0err = get_t0err([amu88],sqrtt0_ph_Regensburg)[1]
    amu88err   = sqrt(err(amu88)^2 + amu88syst^2 + amu88t0err^2)

    amuCC, infoCC = BDIOread_MAtot(path_bdio,diag,wind,"gCCconn",StdDer=STD_DERIV,BLIND=false); uwerr(amuCC)
    der_mDs = mchist(amuCC, "MD_ph [GeV]")[1] / artificial_err
    MD_ph_prime = BDIOread_mDs(path_bdio)["mDs ph"]
    amuCC = amuCC + value(MD_ph - MD_ph_prime) * der_mDs; uwerr(amuCC)
    amuCCsyst  = infoCC["syst"]
    amuCCt0err = get_t0err([amuCC],sqrtt0_ph_Regensburg)[1]
    amuCCerr   = sqrt(err(amuCC)^2 + amuCCsyst^2 + amuCCt0err^2)

    amuW = amu33 + (1/3) * amu88 + (4/9) * amuCC; uwerr(amuW)
    amuWsyst  = sqrt(amu33syst^2 + 1/9 * amu88syst^2 + 16/81 * amuCCsyst^2)
    amuWt0err = get_t0err([amuW],sqrtt0_ph_Regensburg)[1]
    amuWerr   = sqrt(err(amuW)^2 + amuWsyst^2 + amuWt0err^2)

    amu[wind]      = amuW
    amusyst[wind]  = amuWsyst
    amut0err[wind] = amuWt0err
    amuerr[wind]   = amuWerr
end


BLIND_LD_factor = BLIND_LD[1] ? BLIND_LD[2] : 1.0

AMU      = amu["SD"] + amu["ID"] + amu["LD"]/BLIND_LD_factor; uwerr(AMU)
AMUSYST  = sqrt(amusyst["SD"]^2 + amusyst["ID"]^2 + (amusyst["LD"]/BLIND_LD_factor)^2)
AMUT0ERR = get_t0err([AMU],sqrtt0_ph_Regensburg)[1]
AMUERR   = sqrt(err(AMU)^2 + AMUSYST^2 + AMUT0ERR^2)

digits = 4
factor = (diag == "LO" ? 1 : 10)

println(" Window decomposition:")
println("---------------------------------------------------------")
println(" => amu($diag;SD) = $(round(value(amu["SD"])*factor,digits=digits))($(round(err(amu["SD"])*factor,digits=digits)))($(round(amusyst["SD"]*factor,digits=digits)))($(round(amut0err["SD"]*factor,digits=digits)))[$(round(amuerr["SD"]*factor,digits=digits))]")
println(" => amu($diag;ID) = $(round(value(amu["ID"])*factor,digits=digits))($(round(err(amu["ID"])*factor,digits=digits)))($(round(amusyst["ID"]*factor,digits=digits)))($(round(amut0err["ID"]*factor,digits=digits)))[$(round(amuerr["ID"]*factor,digits=digits))]")
println(" => amu($diag;LD) = $(round(value(amu["LD"])*factor,digits=digits))($(round(err(amu["LD"])*factor,digits=digits)))($(round(amusyst["LD"]*factor,digits=digits)))($(round(amut0err["LD"]*factor,digits=digits)))[$(round(amuerr["LD"]*factor,digits=digits))]")


println(" Final estimation for diag $diag")
println("---------------------------------------------------------")
println(" => amu($diag;$comp) = $(round(value(AMU)*factor,digits=digits))($(round(err(AMU)*factor,digits=digits)))($(round(AMUSYST*factor,digits=digits)))($(round(AMUT0ERR*factor,digits=digits)))[$(round(AMUERR*factor,digits=digits))]")