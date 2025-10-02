# Import packages

using Revise

include("HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs

using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

using Setfield

# include uwreal constants

# include("HVPtool/uwConst.jl")

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_bPert = joinpath(julia_script_directory, "..", "PertSD")

# Dict to relate contributions to their Dict keys

DictComptoKey = Dict{String,Vector{String}}(
    "g33"      => ["g33_ll","g33_lc"],
    "g88"      => ["g88_ll","g88_lc"],
    # only interested in the lc (local-conserved) discr. for the cc conn
    "gCCconn"  => ["gCCconn_SU3_lc"],
    # "gCCconn"  => ["gCCconn_SU3_ll","gCCconn_SU3_lc"],

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

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> Model Average <==========================##

diag = "NLOa&b"  #  LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
wind = "ID"  #  NW  SD  SDsub  ID  LD  ILD
comp = "g88"  #  g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b

Q = 5.0  # virtuality for SDsub

model_var_list = Function[a3,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log,phi2inv,logphi2]
# model_var_list = Function[a3,a2y,ysqr,ylog,yinv,logy]
MultFunc = nothing  #  nothing  deltaphi

IMPR_SET = ["1","2"]  #  ["1"]  ["1old"]  ["2"]  ["1","2"]  ["1old","2"]  ["1","1old","2"]

FITCUT = ["None","mass","beta","beta&mass"]  #  ["None","mass","beta","beta&mass"]

BLIND = false

STD_DERIV  = false
tl_IMPR    = false
VREF       = false
RESC       = false

SimpleBase = false

WRITE      = true
OVERWRITE  = false

path_bdio_r = path_bdio_dict["local"]
path_bdio_w = path_bdio_dict["local"]

if comp != "g33" && VREF
    @error("Cannot project to Vref for chosen iso-spin")
end

for FitCut in FITCUT
    if wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]
        @info(" STARTING MA [diag. $diag; wind. $wind; comp. $comp] \n - Vref: $VREF  \n - Rescal: $RESC  \n - Fit type: $FitCut \n - Q: $Q")
    else
        @info(" STARTING MA [diag. $diag; wind. $wind; comp. $comp] \n - Vref: $VREF  \n - Rescal: $RESC  \n - Fit type: $FitCut")
    end

    mykeys = DictComptoKey[comp]

    xdata = Matrix[]; ydata = Dict()
    FitResVec = Dict()
    res =  Dict(); par = Dict()
    modelinfo = Dict()
    for impr_set in IMPR_SET
        println("- Reading Fit data for impr. set $impr_set")

        println("   - Reading X & Y data...")

        xdata, ydata[impr_set] = BDIOread_XYdata(path_bdio_r,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

        println("   - Reading FitRes...")

        FitResVec[impr_set] = JDL2read_FitRes(path_bdio_r,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

        println("   - Reading results and parameters...")

        res[impr_set] = BDIOread_res(path_bdio_r,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,param=false,Q=Q)

        println("   - Reading model information...")

        modelinfo = JDL2read_ModelInfo(path_bdio_r,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)
    end

    # @info(" ⟹ Data ready to compute MA")
    println("- Computing MA")

    # Data ready to compute the MA

    for impr_set in IMPR_SET
        fitcat = Dict()
        weight = Dict()
        amu = Dict{String,Array{uwreal}}()
        syst = Dict()
        for key in mykeys
            println("   - $key")
            str = "$(diag)_$(key)_set$(impr_set)"

            fitcat[key] = FitCat(xdata, ydata[impr_set][key],str)
            for fit in FitResVec[impr_set][key]
                push!(fitcat[key].fit,fit)
            end
            weight[key] = get_w_from_fitcat([fitcat[key]], norm=false)
            amu[key], syst[key] = model_average(res[impr_set][key], weight[key]); uwerr.(amu[key])

            digits = comp in ["cc disc","c8 disc"] ? 7 : 3
            println("      ⟹ amu[$diag($wind)|$comp:$impr_set:$(key[end-1:end])] = $(round(value(amu[key][1]),digits=digits))($(round(err(amu[key][1]),digits=digits)))($(round(syst[key],digits=digits)))")
        end

        if WRITE
            println("- Writing BDIO & JLD2...")

            model_str = func_str(model_var_list,Order=true)

            tl_str  = (tl_IMPR && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
            SIMstr  = SimpleBase ? "SIMPLE" : ""
            MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
            SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
            IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
            RESstr = RESC ? "_resc" : ""
            VREFstr = VREF ? "_Vref" : ""
            DERstr  = STD_DERIV ? "_std" : "" 
            BLINstr = BLIND ? "Blind" : ""
            JDL2info = create_path(path_bdio_w,["Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCut)]","MA","$(BLINstr)MAinfo$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jld2"],OVERWRITE=OVERWRITE)
            pBDIOMA  = create_path(path_bdio_w,["Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCut)]","MA","$(BLINstr)MAres$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr)"],OVERWRITE=OVERWRITE)

            MADict = Dict{String, Any}(
                # "res_vec_jdl2" => res[impr_set],
                # "amu_jld2" => amu,
                "weight" => weight,
                "syst" => syst
            )

            save(JDL2info,"MAinfo",MADict)

            io = IOBuffer()
            write(io, "MA")
            fb = ALPHAdobs_create(pBDIOMA, io)

            extra = Dict{String, Any}("diag" => diag, "wind" => wind, "comp" => comp, "impr set" => impr_set)
            ALPHAdobs_write(fb, res[impr_set], extra=extra)
            ALPHAdobs_write(fb, amu, extra=extra)

            ALPHAdobs_close(fb)
        end
    end # end impr_set loop
    println("\n ")
end # end FitCut loop


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> ENSEMBLE CUT COMBINATION <==========================##

diag = "NLOa&b"  #  LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
wind = "ID"  #  NW  SD  SDsub  ID  LD  ILD
comp = "g88"  #  g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b

Q = 5.0  # virtuality for SDsub

BLIND = false

discr = "all"  #  all  ll  lc

MultFunc = nothing  #  nothing  deltaphi

IMPR_SET = ["1","2"]  #  ["1"]  ["1old"]  ["2"]  ["1","2"]  ["1old","2"]  ["1","1old","2"]

STD_DERIV  = false
tl_IMPR    = false
VREF       = false
RESC       = false

WRITE      = true
OVERWRITE  = false

path_bdio_r = path_bdio_dict["local"]
path_bdio_w = path_bdio_dict["local"]


# SD:
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [false,Function[a3,a4,a2phi2,phi2sqr,phi2log]],
#     "beta"      => [false,Function[a3,a4,a2phi2,phi2sqr,phi2log]],
#     "mass"      => [false,Function[a3,a4,a2phi2,phi2sqr,phi2log]],
#     "beta&mass" => [false,Function[a3,a4,a2phi2,phi2sqr,phi2log]],
# )
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [false,Function[a3,a2phi2,phi2sqr,phi2log]],
#     "beta"      => [false,Function[a3,a2phi2,phi2sqr,phi2log]],
# )
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [true,Function[a3,a4,a2phi2,phi2sqr,phi2log,phi4]],
#     "beta"      => [true,Function[a3,a2phi2,phi2sqr,phi2log,phi4]],
#     "mass"      => [true,Function[a3,a4,a2phi2,phi2sqr,phi2log,phi4]],
#     "beta&mass" => [true,Function[a3,a2phi2,phi2sqr,phi2log,phi4]],
# )
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [false,Function[a3,a4,a2phi2,phi2sqr,phi2log]],
#     "beta"      => [false,Function[a3,a4,a2phi2,phi2sqr,phi2log]],
# )

# ID:
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [false,Function[a3,a2phi2,phi2sqr,phi2log,phi2inv,logphi2]],
#     "beta"      => [false,Function[a3,a2phi2,phi2sqr,phi2log,phi2inv,logphi2]],
#     "mass"      => [false,Function[a3,a2phi2,phi2sqr,phi2log,phi2inv,logphi2]],
#     "beta&mass" => [false,Function[a3,a2phi2,phi2sqr,phi2log,phi2inv,logphi2]],
# )
FITCUTtoMODEL = Dict{String,Any}(
    "None"      => [false,Function[a3,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log,phi2inv,logphi2]],
    "beta"      => [false,Function[a3,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log,phi2inv,logphi2]],
    "mass"      => [false,Function[a3,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log,phi2inv,logphi2]],
    "beta&mass" => [false,Function[a3,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log,phi2inv,logphi2]],
)
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [false,Function[a3,a2y,ysqr,ylog,yinv,logy]],
#     "beta"      => [false,Function[a3,a2y,ysqr,ylog,yinv,logy]],
#     "mass"      => [false,Function[a3,a2y,ysqr,ylog,yinv,logy]],
#     "beta&mass" => [false,Function[a3,a2y,ysqr,ylog,yinv,logy]],
# )
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [false,Function[a3,a2phi2,phi2sqr,phi2log]],
#     "beta"      => [false,Function[a3,a2phi2,phi2sqr,phi2log]],
#     "mass"      => [false,Function[a3,a2phi2,phi2sqr,phi2log]],
#     "beta&mass" => [false,Function[a3,a2phi2,phi2sqr,phi2log]],
# )
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [false,Function[a3,a2y,ysqr,ylog]],
#     "beta"      => [false,Function[a3,a2y,ysqr,ylog]],
#     "mass"      => [false,Function[a3,a2y,ysqr,ylog]],
#     "beta&mass" => [false,Function[a3,a2y,ysqr,ylog]],
# )
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [true,Function[a3,a2phi2,phi2sqr,phi2log,phi4]],
#     "beta"      => [true,Function[a3,a2phi2,phi2sqr,phi2log,phi4]],
#     "mass"      => [true,Function[a3,a2phi2,phi2sqr,phi2log,phi4]],
#     "beta&mass" => [true,Function[a3,a2phi2,phi2sqr,phi2log,phi4]],
# )

# LD
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [false,Function[a3,a2phi2,phi2sqr,phi2log,phi2inv,logphi2]],
#     "beta"      => [false,Function[a3,a2phi2,phi2sqr,phi2log,phi2inv,logphi2]],
#     "mass"      => [false,Function[a3,a2phi2,phi2sqr,phi2log,phi2inv,logphi2]],
#     "beta&mass" => [false,Function[a3,a2phi2,phi2sqr,phi2log,phi2inv,logphi2]],
# )
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [false,Function[a2y,ysqr,ylog,yinv,logy]],
#     "beta"      => [false,Function[a2y,ysqr,ylog,yinv,logy]],
#     "mass"      => [false,Function[a2y,ysqr,ylog,yinv,logy]],
#     "beta&mass" => [false,Function[a2y,ysqr,ylog,yinv,logy]],
# )
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [false,Function[a3,a2phi2,phi2sqr,phi2log]],
#     "beta"      => [false,Function[a3,a2phi2,phi2sqr,phi2log]],
#     "mass"      => [false,Function[a3,a2phi2,phi2sqr,phi2log]],
#     "beta&mass" => [false,Function[a3,a2phi2,phi2sqr,phi2log]],
# )
# FITCUTtoMODEL = Dict{String,Any}(
#     "None"      => [true,Function[a3,a2phi2,phi2sqr,phi2log,phi4]],
#     "beta"      => [true,Function[a3,a2phi2,phi2sqr,phi2log,phi4]],
#     "mass"      => [true,Function[a3,a2phi2,phi2sqr,phi2log,phi4]],
#     "beta&mass" => [true,Function[a3,a2phi2,phi2sqr,phi2log,phi4]],
# )


if wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]
    @info(" STARTING MA [diag. $diag; wind. $wind; comp. $comp] \n - Vref  : $VREF  \n - Rescal: $RESC \n - Fit cut combination: $(keys(FITCUTtoMODEL)) \n - Q: $Q")
else
    @info(" STARTING MA [diag. $diag; wind. $wind; comp. $comp] \n - Vref  : $VREF  \n - Rescal: $RESC \n - Fit cut combination: $(keys(FITCUTtoMODEL))")
end

mykeys = DictComptoKey[comp]
# mykeys = discr == "all" ? DictComptoKey[comp] : DictComptoKey[comp][discr .== ["ll","lc"]]


factor = Dict(
    "g33" => 1., "g88" => 1/3., "gCCconn" => 4/9., "∆ls_amu" => 1/3., "∆lc_b" => 4/9.,
    "g33s" => "", "g88s" => "(1/3)", "gCCconns" => "(4/9)", "∆ls_amus" => "(1/3)", "∆lc_bs" => "(4/9)"
)

if comp != "g33" && VREF
    @error("Cannot project to Vref for chosen iso-spin")
end

if comp == "gCCconn"
    MD_ph_prime = BDIOread_mDs(path_bdio_r)
end

# FITCUT = keys(FITCUTtoMODEL)
FITCUT = sort(collect(keys(FITCUTtoMODEL)), by = x -> Dict(s => i for (i, s) in enumerate(["None","beta","mass","beta&mass"]))[x])

xdata     = Dict()
ydata     = Dict()
FitResVec = Dict()
res       = Dict()
modelinfo = Dict()

res_tot = Dict()
weight  = Dict()
amu     = Dict()
syst    = Dict()

for impr_set in IMPR_SET
    println("- Impr. set $impr_set")

    xdata[impr_set]     = Matrix[]
    ydata[impr_set]     = Dict()
    FitResVec[impr_set] = Dict()
    res[impr_set]       =  Dict()
    modelinfo[impr_set] = Dict()
    for FitCut in FITCUT
        println("   - Fit cut $FitCut")

        println("      - Reading X & Y data...")

        xdata[FitCut], ydata[impr_set][FitCut] = BDIOread_XYdata(path_bdio_r,diag,wind,comp,FITCUTtoMODEL[FitCut][2],FitCut,impr_set;resc=RESC,SimpleBase=FITCUTtoMODEL[FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

        println("      - Reading FitRes...")

        FitResVec[impr_set][FitCut] = JDL2read_FitRes(path_bdio_r,diag,wind,comp,FITCUTtoMODEL[FitCut][2],FitCut,impr_set;resc=RESC,SimpleBase=FITCUTtoMODEL[FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

        println("      - Reading results...")
        try
            res[impr_set][FitCut], info = BDIOread_MA(path_bdio_r,diag,wind,comp,FITCUTtoMODEL[FitCut][2],FitCut,impr_set,resc=RESC,SimpleBase=FITCUTtoMODEL[FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)
        catch
            @warn("     Reading results from fits because MA file could not be used, this will slow down the reading")
            res[impr_set][FitCut] = BDIOread_res(path_bdio_r,diag,wind,comp,FITCUTtoMODEL[FitCut][2],FitCut,impr_set;resc=RESC,SimpleBase=FITCUTtoMODEL[FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,param=false,Q=Q)
        end

        println("      - Reading model information...")

        modelinfo[FitCut] = JDL2read_ModelInfo(path_bdio_r,diag,wind,comp,FITCUTtoMODEL[FitCut][2],FitCut,impr_set;resc=RESC,SimpleBase=FITCUTtoMODEL[FitCut][1],MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)
    end

    println("- Computing MA")

    fitcat = Dict()
    res_tot[impr_set] = Dict{String,Array{uwreal}}()
    weight[impr_set]  = Dict()
    amu[impr_set]     = Dict{String,Array{uwreal}}()
    syst[impr_set]    = Dict()
    for key in mykeys
        println("   - $key")
        str = "$(diag)_$(mykeys)_set$(impr_set)"

        fitcat[key] = Dict()
        for FitCut in FITCUT
            fitcat[key][FitCut] = FitCat(xdata[FitCut], ydata[impr_set][FitCut][key],str)
            for fit in FitResVec[impr_set][FitCut][key]
                push!(fitcat[key][FitCut].fit,fit)
            end
        end
        weight[impr_set][key]  = get_w_from_fitcat([fitcat[key][FitCut] for FitCut in FITCUT], norm=false)
        res_tot[impr_set][key] = vcat([res[impr_set][FitCut]["res_tot"][key] for FitCut in FITCUT]...)
        amu[impr_set][key], syst[impr_set][key] = model_average(res_tot[impr_set][key], weight[impr_set][key]); uwerr.(amu[impr_set][key])

        # if comp == "gCCconn"
        #     der_mDs = mchist(amu[impr_set][key][1], "MD_ph [GeV]")[1] / artificial_err
        #     amu[impr_set][key][1] = amu[impr_set][key][1] + value(MD_ph - MD_ph_prime) * der_mDs; uwerr(amu[impr_set][key][1])
        # end

        digits = comp in ["gCCdisc","gC8disc"] ? 7 : 5
        println("     ⟹ $(factor[comp*"s"]) amu[$diag($wind)|$comp:$impr_set:$(key[end-1:end])] = $(round(value((factor[comp])*amu[impr_set][key][1]),digits=digits))($(round((factor[comp])*err(amu[impr_set][key][1]),digits=digits)))($((factor[comp])*round(syst[impr_set][key],digits=digits)))")
    end


end
println("----------------------------------------------------------------------------")
println("- COMPUTING FINAL RESULT:")

WEIGHT  = vcat([vcat([weight[impr_set][key] for key in mykeys]...)  for impr_set in IMPR_SET]...)./(length(mykeys)*length(IMPR_SET))
RES_TOT = vcat([vcat([res_tot[impr_set][key] for key in mykeys]...) for impr_set in IMPR_SET]...)

AMU, SYST = model_average(RES_TOT, WEIGHT); uwerr.(AMU)

# if comp == "gCCconn"
#     der_mDs = mchist(AMU[1], "MD_ph [GeV]")[1] / artificial_err
#     AMU[1] = AMU[1] + value(MD_ph - MD_ph_prime) * der_mDs; uwerr(AMU[1])
# end

dig = comp in ["gCCdisc","gC8disc"] ? 7 : 5
println("     ⟹ $(factor[comp*"s"]) amu[$diag($wind)|$comp] = $(round((factor[comp])*value(AMU[1]),digits=dig))($(round((factor[comp])*err(AMU[1]),digits=dig)))($(round((factor[comp])*SYST,digits=dig)))[$(round((factor[comp])*sqrt(err(AMU[1])^2+SYST^2),digits=dig))]")

if WRITE
    println("- Writing BDIO & JLD2...")
    println("   - Writing partial results...")
    for impr_set in IMPR_SET

        tl_str  = (tl_IMPR && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
        SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
        IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
        RESstr = RESC ? "_resc" : ""
        VREFstr = VREF ? "_Vref" : ""
        DERstr  = STD_DERIV ? "_std" : "" 
        BLINstr = BLIND ? "Blind" : ""
        JDL2info = create_path(path_bdio_w,["Fit&MA",diag,wind,comp*tl_str,"MAtot","$(BLINstr)MAinfo$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jld2"],OVERWRITE=OVERWRITE)
        pBDIOMA  = create_path(path_bdio_w,["Fit&MA",diag,wind,comp*tl_str,"MAtot","$(BLINstr)MAres$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr)"],OVERWRITE=OVERWRITE)

        MADict = Dict{String, Any}(
            # "res_vec_jdl2" => res[impr_set],
            "amu_jld2" => amu[impr_set],
            "weight" => weight[impr_set],
            "syst" => syst[impr_set],
            "FITCUTtoMODEL" => FITCUTtoMODEL,
            "MultFunc" => MultFunc,
            "Keys" => mykeys
        )

        save(JDL2info,"MAinfo",MADict)

        io = IOBuffer()
        write(io, "MA")
        fb = ALPHAdobs_create(pBDIOMA, io)

        extra = Dict{String, Any}("diag" => diag, "wind" => wind, "comp" => comp, "impr set" => impr_set)
        ALPHAdobs_write(fb, res_tot[impr_set], extra=extra)
        ALPHAdobs_write(fb, amu[impr_set], extra=extra)

        ALPHAdobs_close(fb)
    end
    println("   - Writing FINAL results...")

    tl_str  = (tl_IMPR && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    RESstr = RESC ? "_resc" : ""
    VREFstr = VREF ? "_Vref" : ""
    DERstr  = STD_DERIV ? "_std" : "" 
    BLINstr = BLIND ? "Blind" : ""
    JDL2info = create_path(path_bdio_w,["Fit&MA",diag,wind,comp*tl_str,"MAtot","$(BLINstr)MAINFO$(SUBQstr)$(RESstr)$(VREFstr)$(DERstr).jld2"],OVERWRITE=OVERWRITE)
    pBDIOMA  = create_path(path_bdio_w,["Fit&MA",diag,wind,comp*tl_str,"MAtot","$(BLINstr)MARES$(SUBQstr)$(RESstr)$(VREFstr)$(DERstr)"],OVERWRITE=OVERWRITE)

    MADict = Dict{String, Any}(
        # "res_vec_jdl2" => RES_TOT,
        "amu_jld2" => AMU,
        "weight" => WEIGHT,
        "syst" => SYST,
        "FITCUTtoMODEL" => FITCUTtoMODEL,
        "MultFunc" => MultFunc,
        "Keys" => mykeys,
        "IMPR_SET" => IMPR_SET
    )

    save(JDL2info,"MAinfo",MADict)

    io = IOBuffer()
    write(io, "MA")
    fb = ALPHAdobs_create(pBDIOMA, io)

    extra = Dict{String, Any}("diag" => diag, "wind" => wind, "comp" => comp)
    # ALPHAdobs_write(fb, RES_TOT, extra=extra)
    ALPHAdobs_write(fb, AMU, extra=extra)

    ALPHAdobs_close(fb)
end

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> READING TEST <==========================##

diag = "NLOa"  # LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
wind = "ID"  # NW  SD  ID  LD  ILD
comp = "g33"  #  g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b

Q = 5.0  # virtuality for SDsub

BLIND = false

model_var_list = Function[a3,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log,phi4]
MultFunc = nothing # nothing  deltaphi

impr_set = "2"

FitCut = "beta&mass"  # "None"  "beta"  "mass"  "beta&mass"

STD_DERIV  = false
tl_IMPR    = false
VREF       = true
RESC       = false

SimpleBase = true

amu, info = BDIOread_MA(path_bdio_w,diag,wind,comp,model_var_list,FitCut,impr_set,resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)
res = amu["res"]["gCCconn_SU3_lc"]
syst = info["syst"]["gCCconn_SU3_lc"]
uwerr(res); 10 .* [value(res),err(res),syst]

##--- Full MA

diag = "NLOa&b"  # LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
wind = "ID"  # NW  SD  SDsub  ID  LD  ILD
comp = "g33"  #  g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b

Q = 5.0  # virtuality for SDsub

BLIND = false

STD_DERIV  = false
tl_IMPR    = false
VREF       = true
RESC       = false

path_bdio_r = path_bdio_dict["local"]

AMU, INFO = BDIOread_MAtot(path_bdio_r,diag,wind,comp,resc=RESC,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

uwerr(AMU); [AMU.mean,AMU.err,INFO["syst"]].*10

# get_t0err([AMU],sqrtt0_ph_Regensburg)[1]

# sqrt(INFO["syst"]^2 +  get_t0err([AMU],sqrtt0_ph_Regensburg)[1]^2 + err(AMU)^2)

##

impr_set = "2"

AMU_set, INFO_set = BDIOread_MAtot(path_bdio_r,diag,wind,comp,read="impr",impr_set=impr_set,StdDer=STD_DERIV,tlImpr=tl_IMPR,BLIND=BLIND,Q=Q)

uwerr(AMU_set); [AMU_set,INFO_set["syst"]]
AMU_set["res"]
INFO_set["syst"]


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> Fit + MA check <==========================##

# include("HVPtools/Reader.jl")

diag = "NLOa&b"  # LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
wind = "ID"  # NW  SD  SDsub  ID  LD  ILD
comp = "g33"  #  g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b

impr_set = "1"

Q = 5.0  # virtuality for SDsub

BLIND = false

model_var_list = [a3,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log,phi2inv,logphi2]
# model_var_list = [a3,a4,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log,phi4]
# model_var_list = [a3,a2y,ysqr,ylog,yinv,logy]
MultFunc = nothing  #  nothing  deltaphi

FitCut = "None"  # "None"  "beta"  "mass"  "beta&mass"

STD_DERIV  = false
tl_IMPR    = false
VREF       = true
RESC       = false

SimpleBase = false

path_bdio = path_bdio_dict["local"]

fitRes = JDL2read_FitRes(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

modelInfo = JDL2read_ModelInfo(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

res, par = BDIOread_res(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q,param=true)

# Fit check

amu, info = BDIOread_MA(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set,resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

@info("All data ready")

##

discr = "ll"  #  ll  lc  cc

key = length(DictComptoKey[comp])==1 ? DictComptoKey[comp][1] : DictComptoKey[comp][discr .== ["ll","lc"]][1]

w     = info["weight"][key]
chir  = getfield.(fitRes[key],:chi2)./getfield.(fitRes[key],:chi2exp)
amu   = res[key]
param = par[key]

argw = sortperm(w, rev=true)

@info("diag: $diag\t wind: $wind\t comp: $comp\t FitCut: $FitCut\t set: $impr_set\t discr: $discr")

println("aµ;\t\t chi2/chi2exp => w;\t Model label")
println("----------------------------------------------------------------------------------------------------------------------------------------------------------------------")
for i in argw[w[argw] .> w[argw[1]]/40]
    uwerr(amu[i]); uwerr.(param[i])
    par_str = "$(round(abs(value(param[i][1])/err(param[i][1])),digits=2))"
    for j=2:length(param[i])
        par_str *= ",$(round(abs(value(param[i][j])/err(param[i][j])),digits=2))"
    end
    println("$(round(value(amu[i]),digits=3))±$(round(err(amu[i]),digits=3))  \t $(round(chir[i],digits=2)) => $(round(w[i],digits=3))   \t $(modelInfo["label_tot_isov"][i]) => [$par_str]")
end

# for i in argw[w[argw] .> w[argw[1]]/4]
#     uwerr(amu[i])
#     println("$(round(value(amu[i]),digits=3))±$(round(err(amu[i]),digits=3))  \t $(round(chir[i],digits=2)) => $(round(w[i],digits=3))   \t $(modelInfo["label_tot_isov"][i])")
# end
