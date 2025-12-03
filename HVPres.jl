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
    "g33s" => "", "g88s" => "(1/3)", "gCCconns" => "(4/9)", "∆ls_amus" => "(1/3)", "∆lc_bs" => "(4/9)",
    "g3333" => 1.,  "g3388" => 2/3., "g33CC" => 8/9., "g8888" => 1/9., "g88CC" => 8/27., "gCCCC" => 16/81., 
    "g3333s" => "",  "g3388s" => "(2/3)", "g33CCs" => "(8/9)", "g8888s" => "(1/9)", "g88CCs" => "(8/27)", "gCCCCs" => "(16/81)",
)

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> INDIVIDUAL RESULTS <==========================##


diag = "NLOc"  # LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
wind = "NW"  # NW  SD  SDsub  ID  LD  ILD
comp = "g3333"  # g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b

IMPR = "all"  # "all"  "1"  "1old"  "2"

Q = 5.0  # virtuality for SDsub

BLIND = false

STD_DERIV  = false
tl_IMPR    = false
VREF       = false
RESC       = false

path_bdio = path_bdio_dict["clust"]


if wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]
    @info(" RES [diag. $diag; wind. $wind; comp. $comp] \n - Vref  : $VREF  \n - Rescal: $RESC \n - Q: $Q")
else
    @info(" RES [diag. $diag; wind. $wind; comp. $comp] \n - Vref  : $VREF  \n - Rescal: $RESC")
end

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

scale_ph = !RESC ? sqrtt0_ph_Madrid :  fPi_ph_PDGFLAG
SCALEerr = get_t0err([AMU],scale_ph,resc=RESC)[1]

# if VREF && !BLIND
#     FVC_ChPT = 1.5*JDL2read_FVC_ChPT(path_FVCcont,diag,wind,Q=Q)
#     AMU += FVC_ChPT
# end
if VREF
    FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,wind,Q=Q)
    AMU += FVC_ChPT
end
uwerr(AMU)
# dig = comp in ["gCCdisc","gC8disc"] ? 7 : 5

factor = diag == "LO" ? charge_factor[comp] : charge_factor[comp]*10
sys_vec = !VREF ? [SYSTerr,SCALEerr] : [SYSTerr,SCALEerr,0.1*abs(FVC_ChPT)]

println("     ⟹ $(charge_factor[comp*"s"]) amu[$diag($wind)|$comp] = $(print_uwreal(factor*AMU,factor*sys_vec,total=true))")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> SD RESULTS <==========================##

diag = "NLOa"  #  LO  NLOa  NLOb  NLOa&b

Q33 = 5.0
QCC = 5.0

STD_DERIV = false
tl_IMPR33 = true
VREF33    = true
RESC      = false

path_bdio = path_bdio_dict["clust"]

amu33sub, info33sub = BDIOread_MAtot(path_bdio,diag,"SDsub","g33",resc=RESC,StdDer=STD_DERIV,tlImpr=tl_IMPR33,Vref=VREF33,Q=Q33)
b33Pert = TXTread_bQ(path_bPert,diag)[Qlist .== Q33][1]; uwerr(b33Pert)

if VREF33
    FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,"SDsub",Q=Q33)[1]
    amu33sub += FVC_ChPT
    amuSDfvc = abs(0.1*FVC_ChPT)
end

amu33SD = amu33sub + Window("SD")(0) * b33Pert.mean; uwerr(amu33SD)
amu33SDsyst  = sqrt(info33sub["syst"]^2 + Window("SD")(0)^2 * err(b33Pert)^2)
amu33SDt0err = get_t0err([amu33SD],sqrtt0_ph_Madrid)[1]
amu33SDerr   = !VREF33 ? sqrt(err(amu33SD)^2 + amu33SDsyst^2 + amu33SDt0err^2) :  sqrt(err(amu33SD)^2 + amu33SDsyst^2 + amu33SDt0err^2 + amuSDfvc^2)

∆ls_amu, info∆ls_amu = BDIOread_MAtot(path_bdio,diag,"SDsub","∆ls_amu",resc=RESC,StdDer=STD_DERIV)

amu88SD = amu33SD + ∆ls_amu; uwerr(amu88SD)
amu88SDsyst  = amu33SDsyst + info∆ls_amu["syst"]
amu88SDt0err = get_t0err([amu88SD],sqrtt0_ph_Madrid)[1]
amu88SDerr   = sqrt(err(amu88SD)^2 + amu88SDsyst^2 + amu88SDt0err^2)

amuCCsub, infoCCsub = BDIOread_MAtot(path_bdio,diag,"SDsub","gCCconn",resc=RESC,StdDer=STD_DERIV,Q=QCC)
∆lc_b, info∆lc_b    = BDIOread_MAtot(path_bdio,diag,"SDsub","∆lc_b",resc=RESC,StdDer=STD_DERIV,Q=QCC)
amuCCSD = amuCCsub + Window("SD")(0) * (2*b33Pert.mean + ∆lc_b); uwerr(amuCCSD)
der_mDs = mchist(amuCCSD, "MD_ph [GeV]")[1] / artificial_err
MD_ph_prime, mDs_SU3 = BDIOread_mDs(path_bdio)
amuCCSD += value(MD_ph - MD_ph_prime) * der_mDs; uwerr(amuCCSD) # there's a mistake here!!!
amuCCSDsyst  = sqrt(infoCCsub["syst"]^2 + Window("SD")(0)^2 * (2*err(b33Pert)^2 + info∆lc_b["syst"]^2))
amuCCSDt0err = get_t0err([amuCCSD],sqrtt0_ph_Madrid)[1]
amuCCSDerr   = sqrt(err(amuCCSD)^2 + amuCCSDsyst^2 + amuCCSDt0err^2)

amuSD = amu33SD + (1/3) * amu88SD + (4/9) * amuCCSD; uwerr(amuSD)
amuSDsyst  = sqrt(amu33SDsyst^2 + 1/9 * amu88SDsyst^2 + 16/81 * amuCCSDsyst^2)
amuSDt0err = get_t0err([amuSD],sqrtt0_ph_Madrid)[1]
amuSDerr   = !VREF33 ? sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2) : sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2 + amuSDfvc^2)

# digits = 4
factor = diag == "LO" ? 1 : 10

println(" Iso-spin decomposition:")
println("---------------------------------------------------------")

sys33SD_vec = !VREF33 ? [amu33SDsyst,amu33SDt0err] : [amu33SDsyst,amu33SDt0err,amuSDfvc]

println(" => amu($diag;SD;33) = $(print_uwreal(factor*amu33SD,factor*sys33SD_vec,total=true))")
println(" => (1/3)amu($diag;SD;88) = $(print_uwreal((factor/3)*amu88SD,(factor/3)*[amu88SDsyst,amu88SDt0err],total=true))")
println(" => (4/9)amu($diag;SD;CC) = $(print_uwreal((factor*4/9)*amuCCSD,(factor*4/9)*[amuCCSDsyst,amuCCSDt0err],total=true))")

println(" Final estimation for diag $diag SD window")
println("---------------------------------------------------------")

sysSD_vec = !VREF33 ? [amuSDsyst,amuSDt0err] : [amuSDsyst,amuSDt0err,amuSDfvc]

println(" => amu($diag;SD) = $(print_uwreal(factor*amuSD,factor*sysSD_vec,total=true))")

##==========================> ID & LD RESULTS <==========================##

diag = "NLOa"  #  LO  NLOa  NLOb  NLOa&b
wind = "LD"  #  ID  LD  ILD

STD_DERIV  = false
VREF33     = true
RESC       = false

BLIND = true

path_bdio = path_bdio_dict["clust"]

amu33, info33 = BDIOread_MAtot(path_bdio,diag,wind,"g33",StdDer=STD_DERIV,BLIND=BLIND,Vref=VREF33)
amu33syst  = info33["syst"]
amu33t0err = get_t0err([amu33],sqrtt0_ph_Madrid)[1]

if VREF33
    FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,wind)
    amufvc = abs(0.1*FVC_ChPT)
    if !BLIND
        amu33 += FVC_ChPT; uwerr(amu33)
    end
    # amu33 += 1.5*FVC_ChPT; uwerr(amu33)
end
amu33err = (!VREF33 || BLIND) ? sqrt(err(amu33)^2 + amu33syst^2 + amu33t0err^2) : sqrt(err(amu33)^2 + amu33syst^2 + amu33t0err^2 + amufvc^2)

amu88, info88 = BDIOread_MAtot(path_bdio,diag,wind,"g88",StdDer=STD_DERIV,BLIND=BLIND)
amu88syst  = info88["syst"]
amu88t0err = get_t0err([amu88],sqrtt0_ph_Madrid)[1]
amu88err   = sqrt(err(amu88)^2 + amu88syst^2 + amu88t0err^2)

amuCC, infoCC = BDIOread_MAtot(path_bdio,diag,wind,"gCCconn",StdDer=STD_DERIV,BLIND=false); uwerr(amuCC)
der_mDs = mchist(amuCC, "MD_ph [GeV]")[1] / artificial_err
MD_ph_prime, mDs_SU3 = BDIOread_mDs(path_bdio)
amuCC = amuCC + value(MD_ph - MD_ph_prime) * der_mDs; uwerr(amuCC) # there's a mistake here!!!
amuCCsyst  = infoCC["syst"]
amuCCt0err = get_t0err([amuCC],sqrtt0_ph_Madrid)[1]
amuCCerr   = sqrt(err(amuCC)^2 + amuCCsyst^2 + amuCCt0err^2)

amu = amu33 + (1/3) * amu88 + (4/9) * amuCC; uwerr(amu)
amusyst  = sqrt(amu33syst^2 + 1/9 * amu88syst^2 + 16/81 * amuCCsyst^2)
amut0err = get_t0err([amu],sqrtt0_ph_Madrid)[1]
amuerr   = (!VREF33 || BLIND) ? sqrt(err(amu)^2 + amusyst^2 + amut0err^2) : sqrt(err(amu)^2 + amusyst^2 + amut0err^2 + amufvc^2)

# digits = 6
factor = diag == "LO" ? 1 : 10

println(" Iso-spin decomposition:")
println("---------------------------------------------------------")

sys33_vec = !VREF33 ? [amu33syst,amu33t0err] : [amu33syst,amu33t0err,amufvc]

println(" => amu($diag;$wind;33) = $(print_uwreal(factor*amu33,factor*sys33_vec,total=true))")
println(" => (1/3)amu($diag;$wind;88) = $(print_uwreal((factor/3)*amu88,(factor/3)*[amu88syst,amu88t0err],total=true))")
println(" => (4/9)amu($diag;$wind;CC) = $(print_uwreal((factor*4/9)*amuCC,(factor*4/9)*[amuCCsyst,amuCCt0err],total=true))")

println(" Final estimation for diag $diag $wind window")
println("---------------------------------------------------------")

sys_vec = !VREF33 ? [amusyst,amut0err] : [amusyst,amut0err,amufvc]

println(" => amu($diag;$wind) = $(print_uwreal(factor*amu,factor*sys_vec,total=true))")

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

STD_DERIV  = false
tl_IMPR    = true
VREF       = false
RESC       = false


BLIND = Any[false,1.]  #  [false,1.0]  [true,1.5]

path_bdio = path_bdio_dict["local"]

if VREF && comp != "g33"
    error("Vref cannot be set for comp. $comp")
end
if BLIND[1] && comp != "g33"
    error("Blind only proceeds for comp. g33")
end

# SD

b33Pert = TXTread_bQ(path_bPert,diag)[Qlist .== Q33][1]; uwerr(b33Pert)

if comp == "g33"
    amuSDsub, infoSDsub = BDIOread_MAtot(path_bdio,diag,"SDsub",comp,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,Q=Q33)
    if VREF
        FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,"SDsub",Q=Q33)
        amuSDsub += FVC_ChPT
        amuSDfvc = abs(0.1*FVC_ChPT)
    end
    amuSD = amuSDsub + Window("SD")(0) * value(b33Pert); uwerr(amuSD)
    amuSDsyst  = sqrt(infoSDsub["syst"]^2 + Window("SD")(0)^2 * err(b33Pert)^2)
    amuSDt0err = get_t0err([amuSD],sqrtt0_ph_Madrid)[1]

    amuSDerr   = VREF ? sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2) : sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2 + amuSDfvc^2)
elseif comp == "88"
    ∆ls_amu, info∆ls_amu = BDIOread_MAtot(path_bdio,diag,"SDsub","∆ls_amu",StdDer=STD_DERIV)
    amu33sub, info33sub = BDIOread_MAtot(path_bdio,diag,"SDsub","g33",StdDer=STD_DERIV,tlImpr=tl_IMPR,Q=Q33)
    amu33 = amu33sub + Window("SD")(0) * value(b33Pert)

    amuSD = amu33 + ∆ls_amu; uwerr(amuSD)
    amuSDsyst  = info33["syst"] + info∆ls_amu["syst"]
    amuSDt0err = get_t0err([amuSD],sqrtt0_ph_Madrid)[1]
    amuSDerr   = sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2)
elseif comp == "gCCconn"
    amuSDsub, infoSDsub = BDIOread_MAtot(path_bdio,diag,"SDsub",comp,StdDer=STD_DERIV,tlImpr=tl_IMPR,Q=QCC)
    ∆lc_b, info∆lc_b    = BDIOread_MAtot(path_bdio,diag,"SDsub","∆lc_b",StdDer=STD_DERIV,Q=QCC)
    amuSD = amuSDsub + Window("SD")(0) * (2*value(b33Pert) + ∆lc_b); uwerr(amuSD)
    amuSDsyst  = sqrt(infoSDsub["syst"]^2 + Window("SD")(0)^2 * (2*err(b33Pert)^2 + info∆lc_b["syst"]^2))
    amuSDt0err = get_t0err([amuSD],sqrtt0_ph_Madrid)[1]
    amuSDerr   = sqrt(err(amuSD)^2 + amuSDsyst^2 + amuSDt0err^2)
end

# ID

amuID, infoID = BDIOread_MAtot(path_bdio,diag,"ID",comp,StdDer=STD_DERIV,Vref=VREF); uwerr(amuID)
if VREF
    FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,"ID")
    amuID += FVC_ChPT; uwerr(amuID)
    amuIDfvc = abs(0.1*FVC_ChPT)
end
amuIDsyst  = infoID["syst"]
amuIDt0err = get_t0err([amuID],sqrtt0_ph_Madrid)[1]
amuIDerr   = !VREF ? sqrt(err(amuID)^2 + amuIDsyst^2 + amuIDt0err^2) : sqrt(err(amuID)^2 + amuIDsyst^2 + amuIDt0err^2 + amuIDfvc^2)

# LD

amuLD, infoLD = BDIOread_MAtot(path_bdio,diag,"LD",comp,StdDer=STD_DERIV,Vref=VREF,BLIND=BLIND[1])
amuLD /= BLIND[2]; uwerr(amuLD)
if VREF
    FVC_ChPT = JDL2read_FVC_ChPT(path_FVCcont,diag,"LD")
    amuLD += FVC_ChPT; uwerr(amuLD)
    amuLDfvc = abs(0.1*FVC_ChPT)
end
amuLDsyst  = infoLD["syst"] / BLIND[2]
amuLDt0err = get_t0err([amuLD],sqrtt0_ph_Madrid)[1]
amuLDerr   = !VREF ? sqrt(err(amuLD)^2 + amuLDsyst^2 + amuLDt0err^2) : sqrt(err(amuLD)^2 + amuLDsyst^2 + amuLDt0err^2 + amuLDfvc^2)


amu = amuSD + amuID + amuLD; uwerr(amu)
amusyst  = sqrt(amuSDsyst^2+amuIDsyst^2+amuLDsyst^2)
amut0err = get_t0err([amu],sqrtt0_ph_Madrid)[1]
if VREF
    amufvc = 0.1*JDL2read_FVC_ChPT(path_FVCcont,diag,"NW")
end
amuerr   = !VREF ? sqrt(err(amu)^2 + amusyst^2 + amut0err^2) : sqrt(err(amu)^2 + amusyst^2 + amut0err^2 + amufvc^2)

if comp == "gCCconn"
    MD_ph_prime, mDs_SU3 = BDIOread_mDs(path_bdio)

    der_mDs    = mchist(amu  , "MD_ph [GeV]")[1] / artificial_err
    der_mDs_SD = mchist(amuSD, "MD_ph [GeV]")[1] / artificial_err
    der_mDs_ID = mchist(amuID, "MD_ph [GeV]")[1] / artificial_err
    der_mDs_LD = mchist(amuLD, "MD_ph [GeV]")[1] / artificial_err

    amu   += value(MD_ph - MD_ph_prime) * der_mDs    # ; uwerr(amu)
    amuSD += value(MD_ph - MD_ph_prime) * der_mDs_SD # ; uwerr(amuSD)
    amuID += value(MD_ph - MD_ph_prime) * der_mDs_ID # ; uwerr(amuID)
    amuLD += value(MD_ph - MD_ph_prime) * der_mDs_LD # ; uwerr(amuLD)
end

# digits = 4
factor = (diag == "LO" ? 1 : 10) * (charge_factor[comp])

if !VREF
    systSD = [amuSDsyst,amuSDt0err]
    systID = [amuIDsyst,amuIDt0err]
    systLD = [amuLDsyst,amuLDt0err]
    syst   = [amusyst,amut0err]
else
    systSD = [amuSDsyst,amuSDt0err,amuSDfvc]
    systID = [amuIDsyst,amuIDt0err,amuIDfvc]
    systLD = [amuLDsyst,amuLDt0err,amuLDfvc]
    syst   = [amusyst,amut0err,amufvc]
end

println(" Window decomposition:")
println("---------------------------------------------------------")

println(" => $(charge_factor[comp*"s"]) amu($diag;SD;$comp) = $(print_uwreal(factor*amuSD,factor*systSD,total=true))")
println(" => $(charge_factor[comp*"s"]) amu($diag;ID;$comp) = $(print_uwreal(factor*amuID,factor*systID,total=true))")
println(" => $(charge_factor[comp*"s"]) amu($diag;LD;$comp) = $(print_uwreal(factor*amuLD,factor*systLD,total=true)) \n")


println(" Final estimation for diag $diag comp $comp")
println("---------------------------------------------------------")
println(" => $(charge_factor[comp*"s"]) amu($diag;$comp) = $(print_uwreal(factor*amu,factor*syst,total=true))")


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> NLOc COMP RES <==========================##

diag = "NLOc"

STD_DERIV = false

path_bdio = path_bdio_dict["local"]


amu = Dict(); info = Dict()
amusyt = Dict(); amut0 = Dict()

for comp in ["g3333","g3388","g33CC","g8888","g88CC","gCCCC"]
    amu[comp], info[comp] = BDIOread_MAtot(path_bdio,diag,"NW",comp,resc=false,StdDer=STD_DERIV,tlImpr=false,Vref=false,BLIND=false)
    amusyt[comp] = info[comp]["syst"]
    amut0[comp]  = get_t0err([amu[comp]],sqrtt0_ph_Madrid)[1]
end

MD_ph_prime, mDs_SU3 = BDIOread_mDs(path_bdio)

for comp in ["g33CC","g88CC","gCCCC"]
    der_mDs   = mchist(amu[comp]  , "MD_ph [GeV]")[1] / artificial_err
    amu[comp] += value(MD_ph - MD_ph_prime) * der_mDs
end

AMU = uwreal(0.0); AMUSYST2 = 0.0
println(" Iso-spin decomposition:")
println("---------------------------------------------------------")
for comp in ["g3333","g3388","g33CC","g8888","g88CC","gCCCC"]
    println("$(charge_factor[comp*"s"]) amu($diag;$comp) = $(print_uwreal((10*charge_factor[comp])*amu[comp],(10*charge_factor[comp])*[amusyt[comp],amut0[comp]],total=true))")
    AMU      += (10*charge_factor[comp])*amu[comp]
    AMUSYST2 += ((10*charge_factor[comp])*amusyt[comp])^2
end
println("\n")

AMUSYST = sqrt(AMUSYST2)
AMUt0   = get_t0err([AMU],sqrtt0_ph_Madrid)[1]

println(" Final estimation for diag $diag")
println("---------------------------------------------------------")
println("amu($diag) = $(print_uwreal(AMU,[AMUSYST,AMUt0],total=true))")


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> FULL DIAG <==========================##

diag = "NLOa"  #  LO  NLOa  NLOb  NLOa&b

Q33 = 5.0
QCC = 5.0

STD_DERIV  = false
tl_IMPR    = true
VREF       = true
RESC       = false

BLIND_LD  = Any[true,1.0]  #  Any[true,1.5]  Any[false,1.0]

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
# amu33SDt0err = get_t0err([amu33SD],sqrtt0_ph_Madrid)[1]
# amu33SDerr   = VREF ? sqrt(err(amu33SD)^2 + amu33SDsyst^2 + amu33SDt0err^2) : sqrt(err(amu33SD)^2 + amu33SDsyst^2 + amu33SDt0err^2 + amu33SDfvc^2)

∆ls_amu, info∆ls_amu = BDIOread_MAtot(path_bdio,diag,"SDsub","∆ls_amu",StdDer=STD_DERIV)

amu88SD = amu33SD + ∆ls_amu; uwerr(amu88SD)
amu88SDsyst  = amu33SDsyst + info∆ls_amu["syst"]
# amu88SDt0err = get_t0err([amu88SD],sqrtt0_ph_Madrid)[1]
# amu88SDerr   = sqrt(err(amu88SD)^2 + amu88SDsyst^2 + amu88SDt0err^2)

amuCCsub, infoCCsub = BDIOread_MAtot(path_bdio,diag,"SDsub","gCCconn",StdDer=STD_DERIV,Q=QCC)
∆lc_b, info∆lc_b    = BDIOread_MAtot(path_bdio,diag,"SDsub","∆lc_b",StdDer=STD_DERIV,Q=QCC)
amuCCSD = amuCCsub + Window("SD")(0) * (2*value(b33Pert) + ∆lc_b); uwerr(amuCCSD)
der_mDs = mchist(amuCCSD, "MD_ph [GeV]")[1] / artificial_err
MD_ph_prime, mDs_SU3 = BDIOread_mDs(path_bdio)
amuCCSD += value(MD_ph - MD_ph_prime) * der_mDs; uwerr(amuCCSD)
amuCCSDsyst  = sqrt(infoCCsub["syst"]^2 + Window("SD")(0)^2 * (2*err(b33Pert)^2 + info∆lc_b["syst"]^2))
# amuCCSDt0err = get_t0err([amuCCSD],sqrtt0_ph_Madrid)[1]
# amuCCSDerr   = sqrt(err(amuCCSD)^2 + amuCCSDsyst^2 + amuCCSDt0err^2)

amuSD = amu33SD + (1/3) * amu88SD + (4/9) * amuCCSD; uwerr(amuSD)
amuSDsyst  = sqrt(amu33SDsyst^2 + 1/9 * amu88SDsyst^2 + 16/81 * amuCCSDsyst^2)
amuSDt0err = get_t0err([amuSD],sqrtt0_ph_Madrid)[1]
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
    # amu33t0err = get_t0err([amu33],sqrtt0_ph_Madrid)[1]
    # amu33err   = VREF ? sqrt(err(amu33)^2 + amu33syst^2 + amu33t0err^2) : sqrt(err(amu33)^2 + amu33syst^2 + amu33t0err^2 + amu33fvc^2)

    amu88, info88 = BDIOread_MAtot(path_bdio,diag,wind,"g88",StdDer=STD_DERIV,BLIND=BLIND)
    amu88syst  = info88["syst"]
    # amu88t0err = get_t0err([amu88],sqrtt0_ph_Madrid)[1]
    # amu88err   = sqrt(err(amu88)^2 + amu88syst^2 + amu88t0err^2)

    amuCC, infoCC = BDIOread_MAtot(path_bdio,diag,wind,"gCCconn",StdDer=STD_DERIV,BLIND=false); uwerr(amuCC)
    der_mDs = mchist(amuCC, "MD_ph [GeV]")[1] / artificial_err
    MD_ph_prime, mDs_SU3 = BDIOread_mDs(path_bdio)
    amuCC = amuCC + value(MD_ph - MD_ph_prime) * der_mDs; uwerr(amuCC)
    amuCCsyst  = infoCC["syst"]
    # amuCCt0err = get_t0err([amuCC],sqrtt0_ph_Madrid)[1]
    # amuCCerr   = sqrt(err(amuCC)^2 + amuCCsyst^2 + amuCCt0err^2)

    amuW = amu33 + (1/3) * amu88 + (4/9) * amuCC; uwerr(amuW)
    amuWsyst  = sqrt(amu33syst^2 + 1/9 * amu88syst^2 + 16/81 * amuCCsyst^2)
    amuWt0err = get_t0err([amuW],sqrtt0_ph_Madrid)[1]
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

AMU      = amu["SD"] + amu["ID"] + amu["LD"]/BLIND_LD_factor; uwerr(AMU)
AMUSYST  = sqrt(amusyst["SD"]^2 + amusyst["ID"]^2 + (amusyst["LD"]/BLIND_LD_factor)^2)
AMUT0ERR = get_t0err([AMU],sqrtt0_ph_Madrid)[1]
AMUERR   = sqrt(err(AMU)^2 + AMUSYST^2 + AMUT0ERR^2)
if VREF
    AMUFVC = abs(0.1*JDL2read_FVC_ChPT(path_FVCcont,diag,"NW"))
end

# digits = 4
factor = (diag == "LO" ? 1 : 10)

if !VREF
    systSD = [amusyst["SD"],amut0err["SD"]]
    systID = [amusyst["ID"],amut0err["ID"]]
    systLD = [amusyst["LD"],amut0err["LD"]]/BLIND_LD_factor
    syst   = [AMUSYST,AMUT0ERR]
else
    systSD = [amusyst["SD"],amut0err["SD"],amufvc["SD"]]
    systID = [amusyst["ID"],amut0err["ID"],amufvc["ID"]]
    systLD = [amusyst["LD"],amut0err["LD"],amufvc["LD"]]/BLIND_LD_factor
    syst   = [AMUSYST,AMUT0ERR,AMUFVC]
end

println(" Window decomposition:")
println("---------------------------------------------------------")

println(" => amu($diag;SD) = $(print_uwreal(factor*amuSD,factor*systSD,total=true))")
println(" => amu($diag;ID) = $(print_uwreal(factor*amuID,factor*systID,total=true))")
println(" => amu($diag;LD) = $(print_uwreal(factor*amuLD,factor*systLD,total=true))")


println(" Final estimation for diag $diag")
println("---------------------------------------------------------")
println(" => amu($diag) = $(print_uwreal(factor*amu,factor*syst,total=true))")
