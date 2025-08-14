# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

# Include Isovector Model and Model Average tools

include("isovModel.jl")

include("HVPtools/Fit&MA.jl")

include("HVPtools/Writer.jl")

include("HVPtools/Reader.jl")

include("HVPtools/Utils.jl")

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")

# Dict to relate contributions to their Dict keys

DictComptoKey = Dict{String,Vector{String}}(
    "g33"      => ["g33_ll","g33_lc"],
    "g88"      => ["g88_ll","g88_lc"],
    # only interested in the lc (local-conserved) discr. for the cc conn
    # "CCconn"  => ["gCCconn_lc"], # ["gCCconn_ll","gCCconn_lc"]
    "gCCconn"  => ["gCCconn_SU3_lc"], # ["gCCconn_SU3_ll","gCCconn_SU3_lc"]

    "∆ls_amu" => ["∆ls_amu_ll","∆ls_amu_lc"],
    "∆lc_b"   => ["∆lc_b_ll","∆lc_b_lc"],

    # only interested in the cc (conserved-conserved) discr. for the cc disc  & c8 disc
    "gCCdisc" => ["gCCdisc_cc"],
    "gC8disc" => ["gC8disc_cc"],

    "g3333"    => ["g3333_ll","g3333_lc"],
    "g8888"    => ["g8888_ll","g8888_lc"],
    "gCCCC"    => ["gCCCC_ll","gCCCC_lc"],
    "g3388"    => ["g3388_ll","g3388_lc"],
    "g33CC"    => ["g33CC_ll","g33CC_lc"],
    "g88CC"    => ["g88CC_ll","g88CC_lc"]
)

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> SCALE SETTING ERROR (& SHIFT TO BRUNO) <==========================##

wind = "ID"  # NW  SD  ID  LD  ILD
diag = "LO"  # LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
comp = "g33"  # g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b
key  = "g33_ll"

impr_set = "1old"

Q = 5.0

model_var_list = [a3, a2phi2, a2phi4, a3phi2, phi2sqr, phi2log, phi2inv, phi4sqr, phi4log]
MultFunc = nothing  # nothing  deltaphi

FitCUT = "ß>3.34"  # "All"  "ß>3.34"  "mπ<400"  "ß>3.34&mπ<400"

STD_DERIV  = false
SimpleBase = false

t0_SHIFT = true

sqrtt0_ph_TAR = 0.1443  # 0.1443  sqrtt0_ph_CLS


factor = Dict(
    "33" => 1., "88" => 1/3., "cc conn" => 4/9.,
    "33s" => "", "88s" => "(1/3)", "cc conns" => "(4/9)"
)

res, info = BDIOread_MA(path_bdio,diag,wind,comp,model_var_list,FitCUT,impr_set,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,Q=Q)
amu = Dict(); amu[key] = factor[comp]*res["res"][key]; uwerr(amu[key])

der = mchist(amu[key], "sqrtt0 [fm]")[1] / artificial_err

sqrt_t0_final = t0_SHIFT ? sqrtt0_ph_TAR : sqrtt0_ph_Regensburg

amu_t0shift = amu[key] + value(sqrt_t0_final - sqrtt0_ph_Regensburg) * der

amu_err = amu_t0shift + 0.0; obs = [amu_err]
obs[1] = obs[1] + uwreal([0.0,factor[comp]*info["syst"]],"MA$comp systematics")

t0err = get_t0err([amu_t0shift],sqrt_t0_final)[1]
add_t0_err!(obs, sqrt_t0_final)
uwerr(obs[1])

println("  IMPR. SET.: $IMPR_SET")
println("  ⟹ $(factor[comp*"s"]) amu[$diag;$comp] = $(round(value(amu_t0shift),digits=3))($(round(err(amu_t0shift),digits=3)))($(round(factor[comp]*info["syst"],digits=3)))($(round(t0err,digits=3)))[$(round(err(obs[1]),digits=3))]")


## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> CHARM CONNECRED mDs SHIFT (& SHIFT TO BRUNO) <==========================##

wind = "SDsub"  # NW  SD  ID  LD  ILD
diag = "LO"  # LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
comp = "cc conn"  # cc conn  ∆lc_b

impr_set = "1"

Q = 5.0

model_var_list = [a3, a4, a2phi2, a2phi4, phi2sqr, phi2log]

FitCUT = "ß>3.34"  # "All"  "ß>3.34"  "mπ<400"  "ß>3.34&mπ<400"


STD_DERIV  = false
SimpleBase = false

t0_SHIFT = true

sqrtt0_ph_TAR = uwreal([0.1443,0.0],"faket0")  # 0.1443  sqrtt0_ph_CLS

res, info = BDIOread_MA(path_bdio,diag,wind,comp,model_var_list,FitCUT,impr_set,SimpleBase=SimpleBase,MultFunc=nothing,StdDer=STD_DERIV,Q=Q)
amu = 4/9*res["res"][DictComptoKey[comp][1]]; uwerr(amu)

dermDs  = mchist(amu, "MD_ph [GeV]")[1] / artificial_err

amuC = deepcopy(amu)

# (4/9) * amuC
# (4/9) * value(MD_ph) * der

MD_ph_prime = BDIOread_mDs(path_bdio)["mDs ph"]

amuC_mDscorr = amuC + value(MD_ph - MD_ph_prime) * dermDs; uwerr(amuC_mDscorr)

dert0 = mchist(amuC_mDscorr, "sqrtt0 [fm]")[1] / artificial_err

sqrt_t0_final = t0_SHIFT ? sqrtt0_ph_TAR : sqrtt0_ph_Regensburg

amuC_mDscorr_t0corr = amuC_mDscorr + value(sqrt_t0_final - sqrtt0_ph_Regensburg) * dert0

t0err = get_t0err([amuC_mDscorr_t0corr],sqrt_t0_final)[1]

amuC_err = amuC_mDscorr_t0corr + 0.0; obs = [amuC_err]
obs[1] = obs[1] + uwreal([0.0,4/9*info["syst"]["gcc_lc_conn_beta"]],"MAcc systematics")
add_t0_err!(obs, sqrt_t0_final); uwerr(amuC_mDscorr_t0corr[1])

println("  ⟹ (4/9) amu[$diag:comp] = $(round(value(amuC_mDscorr_t0corr[1]),digits=3))($(round(err(amuC_mDscorr_t0corr[1]),digits=3)))($(round((4/9)*info["syst"]["gcc_lc_conn_beta"],digits=3)))($(round(t0err,digits=3)))[$(round(err(obs[1]),digits=3))]")

(4/9)*info["syst"]["gcc_lc_conn_beta"]

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> RESULT <==========================##

wind  = "NW"  # NW  SD  ID  LD  ILD
order = "NLO"  # LO  NLO  LO&NLO  diagram
diag  = ""    # NLOa  MLOb  NLOc
a_b   = true

type_basemodel = "phi4"  # phi4  simple
type_DA        = "All-(a4,a2loga)"  # All  All-(a4)  All-(a4,a2phi4)  ||  Only(phi2inv,logphi2)  Only(a3cutoff,phi2inv,logphi2)

# Needed dict

DictChargeFac = Dict{String,Dict}(
    "33"      => Dict{String,Any}("n" => 1.0, "s" => ""),
    "88"      => Dict{String,Any}("n" => 1/3., "s" => "1/3"),
    "cc conn" => Dict{String,Any}("n" => 4/9., "s" => "4/9"),
    "cc disc" => Dict{String,Any}("n" => 4/9., "s" => "4/9"),
    "c8 disc" => Dict{String,Any}("n" => 2/(3*sqrt(3)), "s" => "2/3√3"),

    "3333"    => Dict{String,Any}("n" => 1.0, "s" => ""),
    "8888"    => Dict{String,Any}("n" => 1/9., "s" => "1/9"),
    "CCCC"    => Dict{String,Any}("n" => 16/81., "s" => "16/81"),
    "3388"    => Dict{String,Any}("n" => 2/3., "s" => "2/3"),
    "33CC"    => Dict{String,Any}("n" => 8/9., "s" => "8/9"),
    "88CC"    => Dict{String,Any}("n" => 8/27., "s" => "8/27")
)

# start code

if order == "LO"
    DIAG = ["LO"]
elseif order == "NLO"
    if a_b
        DIAG = ["NLOa&b","NLOc"]
    else
        DIAG = ["NLOa","NLOb","NLOc"]
    end
elseif order == "LO&NLO"
    if a_b
        DIAG = ["LO","NLOa&b","NLOc"]
    else
        DIAG = ["LO","NLOa","NLOb","NLOc"]
    end
elseif order == "diagram"
    DIAG = [diag]
end

pDA   = joinpath(path_bdio,"DA")
pDAt0 = joinpath(path_bdio,"DA_t0")

pDAtype   = joinpath(pDA,wind,"base[$type_basemodel]",type_DA)
pDAtypet0 = joinpath(pDAt0,wind,"base[$type_basemodel]",type_DA)

amu   = Dict{String, uwreal}(); syst2   = Dict{String,Float64}()
amut0 = Dict{String, uwreal}(); syst2t0 = Dict{String,Float64}()
for diag in DIAG
    if diag in ["LO", "NLOa", "NLOb", "NLOa&b"]
        if wind in ["NW","SD"]
            COMP = ["33", "88", "cc conn", "cc disc", "c8 disc"]
        elseif wind in ["ID","LD","ILD"]
            COMP = ["33", "88", "cc conn"]
        end
    elseif diag == "NLOc"
        COMP = ["3333", "8888", "CCCC", "3388", "33CC", "88CC"]
    end

    p_ = joinpath(pDAtype,"MA",diag)
    pt0_ = joinpath(pDAtypet0,"MA",diag)

    amu[diag] = uwreal(0.0); syst2[diag] = 0.0
    print(" ⟹ amu($wind)[$diag] =")
    for (ic,comp) in enumerate(COMP)

        # ERRt0 = false

        fb = BDIO_open(joinpath(p_,comp,"MA"),"r")
        res = Dict{String, Any}()
        i=0; mykeys = ["res_tot","amu"]
        while ALPHAdobs_next_p(fb)
            i += 1
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            res[mykeys[i]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)[mykeys[i]][i == 1 ? (1:end) : 1]
        end
        BDIO_close!(fb)

        resinfo = load(joinpath(p_,comp,"MA_info.jld2"), "MAinfo")

        # ERRt0 = true

        fbt0 = BDIO_open(joinpath(pt0_,comp,"MA"),"r")
        rest0 = Dict{String, Any}()
        i=0; mykeys = ["res_tot","amu"]
        while ALPHAdobs_next_p(fbt0)
            i += 1
            d = ALPHAdobs_read_parameters(fbt0)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            rest0[mykeys[i]] = ALPHAdobs_read_next(fbt0, size=sz, keys=ks)[mykeys[i]][i == 1 ? (1:end) : 1]
        end
        BDIO_close!(fbt0)

        # resinfot0 = load(joinpath(pt0_,comp,"MA_info.jld2"), "MAinfo")

        myres = rest0["amu"] * (value(res["amu"])/value(rest0["amu"]))

        amu[diag]   += DictChargeFac[comp]["n"] * myres
        syst2[diag] += (DictChargeFac[comp]["n"] * resinfo["syst"])^2
        if comp in ["cc disc","c8 disc"]
            print(" $(DictChargeFac[comp]["s"]) ["); format_obs(myres*1e5,syst=resinfo["syst"]*1e5); print("x10^(-5)] ")
        else
            print(" $(DictChargeFac[comp]["s"]) ["); format_obs(myres,syst=resinfo["syst"]); print("] ")
        end
        ic != length(COMP) ? print("+") : print("=")
    end
    print(" "); format_obs(amu[diag],syst=sqrt(syst2[diag])); print("\n")
end

if order in ["LO&NLO","NLO"]
    amu_O = uwreal(0.0); syst2_O = 0.0
    print("\n ⟹ amu($wind)[$order] =")
    for (id,diag) in enumerate(DIAG)

        amu_O   += amu[diag]
        syst2_O += syst2[diag]

        print(" "); format_obs(amu[diag],syst=sqrt(syst2[diag]))
        id != length(DIAG) ? print("+") : print("=")
    end
    print(" "); format_obs(amu_O,syst=sqrt(syst2_O)); print("\n")
else
    nothing
end

