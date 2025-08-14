# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

# Include Isovector Model and Model Average tools

include("isovModel.jl")

include("HVPtools/Fit&MA.jl")

include("HVPtools/Reader.jl")

include("HVPtools/Utils.jl")

# BDIO path definition

julia_script_directory = @__DIR__

STD_DERIV = true

if STD_DERIV
    @info("Standard derivative has been chosen !!")
    path_bdio = joinpath(julia_script_directory, "..", "ObsBDIOstd")
else
    path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")
end

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

# ensList = ["H101","B450","N202","N300","H102","N101","C101","S400","N203","N200","D200","N302","E250","J303","E300","J500","A654","N451","D452","J501"]

DictComptoKey = Dict{String,Vector{String}}(
    "33"      => ["g33_ll","g33_lc"],
    "88"      => ["g88_ll","g88_lc"],

    # only interested in the lc (local-conserved) discr for the cc conn
    # "cc conn" => ["gcc_ll_conn","gcc_lc_conn"],
    # "cc conn" => ["gcc_lc_conn"],
    "cc conn" => ["gcc_lc_conn_beta"],

    "cc disc" => ["gcc_cc_disc"],
    "c8 disc" => ["gc8_cc_disc"],

    "3333"    => ["g3333_ll","g3333_lc"],
    "8888"    => ["g8888_ll","g8888_lc"],
    "CCCC"    => ["gCCCC_ll","gCCCC_lc"],
    "3388"    => ["g3388_ll","g3388_lc"],
    "33CC"    => ["g33CC_ll","g33CC_lc"],
    "88CC"    => ["g88CC_ll","g88CC_lc"]
)


## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> Model Average <==========================##

wind = "ID"  # NW, SD, ID, LD, ILD
diag = "LO"  # LO  NLOa  NLOb  NLOa&b  NLOc

base_model = "phi4"  # phi4  simple
perm_model = "All-(a4,a2loga)"  # All  All-(a4)  All-(a4,a2loga)  All-(a4,a2phi4)  ||  Only(phi2inv,logphi2)  Only(a3cutoff,phi2inv,logphi2)

IMPR_SET = ["1old","2"] # ["1","2"] ["1old","2"]

# comp = [33  88  cc conn  cc disc  c8 disc  ||  3333  8888  CCCC  3388  33CC  88CC]
if diag in ["LO", "NLOa", "NLOb", "NLOa&b", "NLOa&b(+)"]
    if wind in ["NW","SD"]
        COMP = ["33", "88", "cc conn", "cc disc", "c8 disc"]
    elseif wind in ["ID","LD","ILD"]
        COMP = ["33", "88", "cc conn"]
    end
elseif diag == "NLOc"
    COMP = ["3333", "8888", "CCCC", "3388", "33CC", "88CC"]
end

@info(" STARTING MA [diag. $diag; wind. $wind]")

for comp in COMP
    IMPR_SET = comp in ["cc disc","c8 disc"] ? IMPR_SET[1] : IMPR_SET
    mykeys = DictComptoKey[comp]

    @info(" MA for comp. $comp")

    println("- Reading Fit data")

    println("   - Reading X & Y data...")

    XYdata = BDIOread_XYdata(path_bdio,diag,wind,comp,base_model,perm_model,IMPR_SET=IMPR_SET)

    xdata = XYdata["xdata"]
    ydata = XYdata["ydata"]

    println("   - Reading FitRes...")

    FitResVec = JDL2read_FitRes(path_bdio,diag,wind,comp,base_model,perm_model)

    println("   - Reading results and parameters...")

    res, par = BDIOread_res(path_bdio,diag,wind,comp,base_model,perm_model;param=true,IMPR_SET=IMPR_SET)

    println("   - Reading model information...")

    modelinfo = JDL2read_ModelInfo(path_bdio,diag,wind,comp,base_model,perm_model)

    println("- Computing Model Average")

    fitcat     = Dict{String, Dict}()
    res_vec    = Vector{Vector{uwreal}}()
    fitcat_vec = Vector{FitCat}()
    for impr_set in IMPR_SET
        fitcat[impr_set] = Dict{String,FitCat}()
        for key in mykeys

            str = "$(diag)_$(mykeys)_set$(impr_set)"

            fitcat[impr_set][key] = FitCat(xdata, ydata[impr_set][key],str)
            FitResVec_ = FitResVec[impr_set][key]
            for fit in FitResVec_
                push!(fitcat[impr_set][key].fit,fit)
            end
            push!(res_vec, res[impr_set][key])
            push!(fitcat_vec, fitcat[impr_set][key])
        end
    end

    println("- Computing weights & final result...")

    result_tot = vcat(res_vec...)
    weight_tot = get_w_from_fitcat(fitcat_vec, norm=false)
    renorm_w!(weight_tot,modelinfo["length"])

    result, syst = model_average(result_tot, weight_tot); uwerr(result)

    digits = comp in ["cc disc","c8 disc"] ? 7 : 3
    println("  ⟹ amu[$diag:$comp] = $(round(value(result),digits=digits))($(round(err(result),digits=digits)))($(round(syst,digits=digits)))")

    println("- Writing BDIO & JLD2...")

    pType  = joinpath(path_bdio,"Fit&MA",wind,"base[$base_model]",perm_model)
    
    pMA = joinpath(pType,"MA")
    !ispath(pMA) ? mkdir(pMA) : nothing

    pdiag = joinpath(pMA,diag)
    !ispath(pdiag) ? mkdir(pdiag) : nothing

    pcomp = joinpath(pdiag,comp)
    !ispath(pcomp) ? mkdir(pcomp) : nothing

    MADict = Dict{String, Any}(
        "res_tot" => result_tot,
        "weight_tot" => weight_tot,
        "amu" => result,
        "syst" => syst
    )

    pMA = joinpath(pcomp,"MAinfo.jld2")
    save(pMA,"MAinfo",MADict)

    res_tot = Dict{String, Array{uwreal}}(
        "res_tot" => result_tot,
    )
    amu = Dict{String, Array{uwreal}}(
        "amu" => [result],
    )

    io = IOBuffer()
    write(io, "$comp MA")
    fb = ALPHAdobs_create(joinpath(pcomp,"MA"), io)

    extra = Dict{String, Any}("comp" => comp)
    ALPHAdobs_write(fb, res_tot, extra=extra)
    ALPHAdobs_write(fb, amu, extra=extra)

    ALPHAdobs_close(fb)

    println("\n")
end

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> SCALE SETTING ERROR (& SHIFT) <==========================##

wind = "ID"  # NW  SD  ID  LD  ILD
diag = "LO"  # LO  NLOa  NLOb  NLOa&b  NLOc

base_model = "phi4"  # phi4  simple
perm_model = "All-(a4,a2loga)"  # All  All-(a4)  All-(a4,a2phi4)  ||  Only(phi2inv,logphi2)  Only(a3cutoff,phi2inv,logphi2)

comp = "88"  # 33  88  cc conn

factor = Dict(
    "33" => 1., "88" => 1/3., "cc conn" => 4/9.,
    "33s" => "", "88s" => "(1/3)", "cc conns" => "(4/9)"
)

amu, info = BDIOread_amu(path_bdio,diag,wind,comp,base_model,perm_model); amu["amu"] = factor[comp]*amu["amu"]; uwerr(amu["amu"])
der = mchist(amu["amu"], "sqrtt0 [fm]")[1] / artificial_err

amu_t0shift = amu["amu"] + value(sqrtt0_ph_CLS - sqrtt0_ph_Regensburg) * der

amu_err = amu_t0shift + 0.0; obs = [amu_err]
obs[1] = obs[1] + uwreal([0.0,factor[comp]*info["syst"]],"MA$comp systematics")

# t0err = get_t0err([amu_t0shift],sqrtt0_ph_Regensburg)[1]
# add_t0_err!(obs, sqrtt0_ph_Regensburg)
t0err = get_t0err([amu_t0shift],sqrtt0_ph_CLS)[1]
add_t0_err!(obs, sqrtt0_ph_CLS)
uwerr(obs[1])

amu_t0shift

println("  ⟹ $(factor[comp*"s"]) amu[$diag;$comp] = $(round(value(amu_t0shift),digits=3))($(round(err(amu_t0shift),digits=3)))($(round(factor[comp]*info["syst"],digits=3)))($(round(t0err,digits=3)))[$(round(err(obs[1]),digits=3))]")


## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> CHARM CONNECRED mDs SHIFT (mDs' -> mDs) <==========================##

wind = "ID"  # NW  SD  ID  LD  ILD
diag = "LO"  # LO  NLOa  NLOb  NLOa&b  NLOc

base_model = "phi4"  # phi4  simple
perm_model = "All-(a4,a2loga)"  # All  All-(a4)  All-(a4,a2loga)  All-(a4,a2phi4)  ||  Only(phi2inv,logphi2)  Only(a3cutoff,phi2inv,logphi2)

amu, info = BDIOread_amu(path_bdio,diag,wind,"cc conn",base_model,perm_model); uwerr(amu["amu"])
dermDs  = mchist(amu["amu"], "MD_ph [GeV]")[1] / artificial_err

amuC = amu["amu"] #+ uwreal([0.0,info["syst"]],"MA syst charm")

# (4/9) * amuC
# (4/9) * value(MD_ph) * der

MD_ph_prime = BDIOread_mDs(path_bdio)["mDs ph"]

amuC_mDscorr = amuC + (MD_ph - MD_ph_prime) * dermDs; uwerr(amuC_mDscorr)

dert0 = mchist(amuC_mDscorr, "sqrtt0 [fm]")[1] / artificial_err

amuC_mDscorr_t0corr = amuC_mDscorr + value(sqrtt0_ph_CLS - sqrtt0_ph_Regensburg) * dert0

factor_amuC_mDscorr_t0corr = 4/9 * amuC_mDscorr_t0corr; uwerr(amuC_mDscorr_t0corr)

t0err = get_t0err([factor_amuC_mDscorr_t0corr],sqrtt0_ph_CLS)[1]

amuC_err = factor_amuC_mDscorr_t0corr + 0.0; obs = [amuC_err]
obs[1] = obs[1] + uwreal([0.0,4/9*info["syst"]],"MAcc systematics")
add_t0_err!(obs, sqrtt0_ph_CLS); uwerr(factor_amuC_mDscorr_t0corr[1])

println("  ⟹ (4/9) amu[$diag:cc conn] = $(round(value(factor_amuC_mDscorr_t0corr),digits=3))($(round(err(factor_amuC_mDscorr_t0corr),digits=3)))($(round((4/9)*info["syst"],digits=3)))($(round((4/9)*t0err,digits=3)))[$(round((4/9)*info["syst"],digits=3))]")


## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> READING TEST <==========================##

wind = "NW"  # NW  SD  ID  LD  ILD
diag = "NLOa&b"  # LO  NLOa  NLOb  NLOa&b  NLOc
comp = "cc conn" # 33  88  cc conn  cc disc  c8 disc

type_basemodel = "phi4"  # phi4  simple
type_DA        = "All-(a4,a2loga)"  # All  All-(a4)  All-(a4,a2phi4)  ||  Only(phi2inv,logphi2)  Only(a3cutoff,phi2inv,logphi2)

data = "res"  # res  info



pMA = joinpath(path_bdio,"DA",wind,"base[$type_basemodel]",type_DA,"MA",diag,comp)
if data == "res"
    fb = BDIO_open(joinpath(pMA,"MA"),"r")
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
elseif data == "info"
    res = load(joinpath(pMA,"MA_info.jld2"), "MAinfo")
end

res


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

