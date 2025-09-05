# Import packages

using Revise

include("HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs

using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

using QuadGK

using ProgressBars
using Suppressor

using StatsBase
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

path_coef  = joinpath(julia_script_directory, "..", "KernelCoeff")

# Blind analysis (Simon K.) safe ensembles: 
# SU(3) sym.  H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ensList = ["A653","A654","B450","C101","C102","D150","D200","D201","D251","D450","D451","D452","E250","E300","F300","H101","H200","J303","J304","J306","J307","J500","J501","N101","N200","N202","N203","N302","N451","N452","S400"] # "H102","N300"

ensInfo = EnsInfo.(ensList)

# We do not have charm or disconnected data for some of the ensembles
ensNOcharm = ["C102","D150","D201","D251","D451","F300","H200","J304","J306","J307","J501","N451","N452"]
ensNOdisc  = ["F300","J306"]

# Dict to relate contributions to their Dict keys

DictComptoKey = Dict{String,Vector{String}}(
    "g33"      => ["g33_ll","g33_lc"],
    "g88"      => ["g88_ll","g88_lc"],
    # only interested in the lc (local-conserved) discr. for the cc conn
    "gCCconn"  => ["gCCconn_SU3_lc"], # ["gCCconn_ll","gCCconn_lc"]
    # "gCCconn"  => ["gCCconn_SU3_ll","gCCconn_SU3_lc"], # ["gCCconn_SU3_ll","gCCconn_SU3_lc"]

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

##==========================> Fits [LO, NLOa&b, NLOa, NLOb, NLOc] <==========================##

diag = "NLOa&b"  #  LO  NLOa  NLOb  NLOc  NLOa&b
wind = "LD"  #  NW  SD  SDsub  ID  LD  ILD

readIMPR_SET = ["1","2"] #  ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

BLIND = true

STD_DERIV  = false
tl_IMPR    = false
VREF       = false
RESC       = false

path_bdio_r = path_bdio_dict["local"]

@info(" Reading and preparing data [diag. $diag; wind. $wind; Resc. $RESC; Vref. $VREF]")
STD_DERIV ? @info("STANDARD DERIVATIVE is being employed in the IMPROVEMENT") : nothing

if tl_IMPR
    wind in ["SD","SDsub"] ? (amu3l = compute_HVPtl0(diag,wind,Qlist,path_coef)) : error("3l improvement cannot be applied for wind = $wind")
end

t0 = Dict(); HVP = Dict(); FVC = Dict()
for ens in ensInfo
    println("- Reading data ensemble: $(ens.id)")

    println("   - Reading t0...")

    t0[ens.id] = BDIOread_t0(path_bdio_r, ens)

    println("   - Reading HVP...")

    HVP[ens.id] = Dict()
    for impr_set in readIMPR_SET
        HVP[ens.id][impr_set], info = BDIOread_HVPens(path_bdio_r,diag,wind,ens,impr_set,info=true,resc=RESC,STD=STD_DERIV,BLIND=BLIND)
    
        apply_syst_HVP!(HVP[ens.id][impr_set],info["HVPsyst"],diag,wind,ens.id)
    end

    println("   - Reading FVC...")

    FVC[ens.id] = BDIOread_FVCens(path_bdio_r,diag,wind,ens,resc=RESC,STD=STD_DERIV,Vref=VREF,BLIND=BLIND)
    
    if !VREF; apply_syst_FVC!(FVC[ens.id],diag,wind,ens.id,IMPR_SET=readIMPR_SET,factor=0.1); end

    println("   - aµ = HVP + FVC")

    HVP_VolCorrect!(HVP[ens.id],FVC[ens.id],diag,IMPR_SET=readIMPR_SET)

    if wind in ["SD","SDsub"] && tl_IMPR
        println("   - Applying tree-level improvement...")

        HVP_3limpr!(HVP[ens.id],amu3l,IMPR_SET=readIMPR_SET,meth="prod")
    end
end

@info(" ⟹ Data ready to be fitted\n\n")


##==========================> Data ready to fit

comp = "gCCconn"  #  g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b

Q = 5.0  # virtuality for SDsub

model_var_list = Function[a3,a2phi2,phi2sqr,phi2log,phi2inv,logphi2]  #  [a3,a4,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log]  [a3,a2phi2,phi2sqr,phi2log]
# model_var_list = Function[a3,a4,a2y,ysqr,ylog]  #  [a3,a4,a2y,a2z,a3y,ysqr,ylog]  [a3,a2y,ysqr,ylog]

MultFunc = nothing  #  nothing  deltaphi

IMPR_SET = readIMPR_SET  #  readIMPR_SET  ["1"]  ["2"]  ["1","2"]  ["1old","2"]  ["1","1old","2"]

FITCUT = ["None","mass","beta","beta&mass"]  #  ["None","mass","beta","beta&mass"]
# FitCut = "beta"  #  None  beta  mass  beta&mass

SimpleBase = true
a2RESC     = true

FitINFO    = false

PVAL       = false

WRITE      = true
OVERWRITE  = false

mdof = 4  # minimum number of d.o.f. allowed

mykeys = DictComptoKey[comp]  #  [DictComptoKey[comp][1]]

path_bdio_w = path_bdio_dict["local"]

# Following the LD paper, when it comes to the iso-vector analysis, the 'untrusted' ensembles are: H105, H200, N300,  N302, S400
# ensExcl = ["H105","H200","N300","N302","S400","B450"] # Good for 33 & 88; ?B450 seems quite bad, not sure why
# ensExcl = ["H105","H200","N300","N302","S400"] # Good for ∆lc(b); the same as 33 but B450 seems to work now?
ensExcl = [] # handy for charmed contributions

# ensExcl = ["N300","N302","S400","H101"] # handy for CCconn LD ¿?
# ensExcl = ["N300","N302","S400","J303"] # handy for CCconn SD ¿?
# ensExcl = ["S400","A654"] # handy for CCconn SD ¿?

if comp != "g33" && VREF
    @error("Cannot project to Vref for chosen iso-spin")
end

for FitCut in FITCUT

    aDP      = FitCut in ["beta","beta&mass"] ? 5 : 6
    aDOF     = comp in ["gCCconn"] ? 1 : 2
    na_max   = aDP - aDOF
    nmPi_max = ((FitCut in ["mass","beta&mass"] || comp == "gCCconn") ? 1 : 2)
    nmK_max  = 1
    # add a limit to "pure" phi2 terms
    if wind in ["LD","ILD"]
        nmPi_max = [nmPi_max,1]
    end

    if wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]
        @info(" STARTING FITS [diag. $diag; wind. $wind; comp. $comp] \n - Vref   : $VREF  \n - Rescal : $RESC \n - Fit cut: $FitCut \n - Q: $Q")
    else
        @info(" STARTING FITS [diag. $diag; wind. $wind; comp. $comp] \n - Vref   : $VREF  \n - Rescal : $RESC \n - Fit cut: $FitCut")
    end

    if FitCut == "None"
        argExcl = []
    elseif FitCut == "beta"
        argExcl = getfield.(ensInfo,:beta) .== 3.34
    elseif FitCut == "mass"
        argExcl = [meson_ens[ensid]["mPi"] * (hbarc * sqrt(t0[ensid])/sqrtt0_ph) * 1e3 for ensid in getfield.(ensInfo,:id)] .> 400
    elseif FitCut == "beta&mass"
        argExcl1 = getfield.(ensInfo,:beta) .== 3.34
        argExcl2 = [meson_ens[ensid]["mPi"] * (hbarc * sqrt(t0[ensid])/sqrtt0_ph) * 1e3 for ensid in getfield.(ensInfo,:id)] .> 400
        argExcl = argExcl1 .| argExcl2
    else
        error("Fit Cut $FitCut was not recognised")
    end

    if (phi4 in model_var_list && !SimpleBase) || any(values(countmap(model_var_list)).>1)
        error(" Same term cannot be fitted twice!")
    end

    ensExcl_type = getfield.(ensInfo[argExcl],:id)

    xdata = []
    ydata = Dict(); [ydata[impr_set] = Dict() for impr_set in IMPR_SET]  # initialize dict

    ensListFit = Vector{EnsInfo}()
    println("- Creating 'x' & 'y' data points...")
    i = 0
    for ens in ensInfo
        if ((comp ∈ ["gCCconn","∆lc_b"]) && (ens.id ∈ ensNOcharm)) || (comp ∈ ["g88","∆ls_amu"] && ens.id ∈ ensNOdisc) || (comp ∈ ["gCCdisc","gC8disc"] && (ens.id ∈ ensNOdisc || ens.kappa_l == ens.kappa_s)) || (comp == "∆ls_amu" && ens.kappa_l == ens.kappa_s && MultFunc == deltaphi) || (ens.id ∈ union(ensExcl,ensExcl_type))
            println("   - [$(ens.id): Either excluded or 0 for contribution $comp]")
            continue
        end
        # println("   [$(ens.id):]")
        i += 1

        mP = BDIOread_mPP(path_bdio_w,ens)
        mπ = mP["mPi"]
        mK = ens.kappa_l != ens.kappa_s ? mP["mK"] : mπ
        if !RESC
            push!(xdata, [1 / (8*t0[ens.id]), 8*t0[ens.id]*mπ^2, 8*t0[ens.id]*(mK^2+0.5*mπ^2)])
        else
            fP = BDIOread_fPS(path_bdio_w,ens)
            fπ = fP["fPi"]
            fK = ens.kappa_l != ens.kappa_s ? fP["fK"] : fπ
            fKπ = 2/3*(fK+0.5*fπ)
            push!(xdata, [1 / (8*t0[ens.id]), mπ^2/(8π*fπ^2), (mK^2+0.5*mπ^2)/(8π*fKπ^2)])

        end

        # push!(xdata, [1 / (8*t0[ens.id]), 8*t0[ens.id]*m_ens[ens.id]["mPi"]^2, 8*t0[ens.id]*(m_ens[ens.id]["mK"]^2+0.5*m_ens[ens.id]["mPi"]^2)])

        for impr_set in IMPR_SET
            for key in mykeys
                i == 1 ? ydata[impr_set][key] = Vector{uwreal}() : nothing
                y = typeof(HVP[ens.id][impr_set][key]) == Vector{uwreal} ? HVP[ens.id][impr_set][key][findfirst(x -> x == Q, Qlist)] : HVP[ens.id][impr_set][key]
                push!(ydata[impr_set][key], y)
            end
        end
        push!(ensListFit,ens)
    end

    # nens = length(ensListFit) # number of data ensembles

    f_tot_isov, n_par_tot_isov, label_tot_isov = call_models(
        model_var_list,
        ensListFit,
        mdof,
        SimpleBase=SimpleBase,
        fPiresc=RESC,
        a2resc=a2RESC,
        MultFunc=MultFunc,
        na_max=na_max,
        nmPi_max=nmPi_max,
        nmK_max=nmK_max
    )

    println("- Fitting points")

    xdata = hcat(xdata...)
    xdata = [xdata[1,:] xdata[2,:] xdata[3,:]]

    res = Dict{String, Dict}(); [res[impr_set] = Dict{String, Array{uwreal}}() for impr_set in IMPR_SET]   # initialize dict
    par = Dict{String, Dict}(); [par[impr_set] = Dict{String, Any}() for impr_set in IMPR_SET]             # initialize dict
    fit = Dict{String, Dict}(); [fit[impr_set] = Dict{String, Any}() for impr_set in IMPR_SET]             # initialize dict

    for impr_set in IMPR_SET
        println("   - Starting set "*impr_set)
        for key in mykeys
            println("      - Fitting for comp. $key...")
            res[impr_set][key] = Vector{uwreal}()
            par[impr_set][key] = Vector{Vector{Vector{Float64}}}()
            fit[impr_set][key] = Vector{FitRes}()
            myFor = FitINFO ? collect(1:length(f_tot_isov)) : ProgressBar(collect(1:length(f_tot_isov)))
            for (n,i) in enumerate(myFor)
                myfit, fitresid = fit_routine(f_tot_isov[i], value.(xdata), ydata[impr_set][key], n_par_tot_isov[i], pval=PVAL, info=false, lineprint=FitINFO, fitRes=true)
                if FitINFO
                    println("n: $n")
                    println("---------------------------------------------------------------------")
                    println("Model permutative structure")
                    println("  $(label_tot_isov[i])")
                    println("---------------------------------------------------------------------")
                    println("Single ensemble contribution to chi2")
                    for i=collect(1:length(fitresid))
                        println("  $(ensListFit[i].id) => $((fitresid[i])^2)")
                    end
                    println("---------------------------------------------------------------------")
                end
                uwerr.(myfit.param)

                x_ph = !RESC ? [0.0 phi2_ph phi4_ph] : [0.0 y_ph z_ph]
                push!(res[impr_set][key], f_tot_isov[i](x_ph, myfit.param)[1])
                push!(par[impr_set][key], [[p.mean,p.err] for p in myfit.param])
                push!(fit[impr_set][key], @set myfit.param = uwreal[])  # param is already safed in a BDIO (or JDL2), to avoid writing the same info twice we save a copy of 'myfit' with a null param field
            end
        end
        if WRITE
            println("      - Printing BDIO & JDL2...")

            model_str = func_str(model_var_list,Order=true)

            tl_str  = (tl_IMPR && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
            SIMstr  = SimpleBase ? "SIMPLE" : ""
            MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
            SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
            IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
            VREFstr = VREF ? "_Vref" : ""
            RESstr = RESC ? "_resc" : ""
            DERstr  = STD_DERIV ? "_std" : ""
            BLINstr = (BLIND && comp != "gCCconn") ? "Blind_" : ""
            pBDIOres    = create_path(path_bdio_w,["Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCut)]","$(BLINstr)Res$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr)"],OVERWRITE=OVERWRITE)
            pBDIOxydata = create_path(path_bdio_w,["Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCut)]","$(BLINstr)XYdata$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr)"],OVERWRITE=OVERWRITE)
            pjld2param  = create_path(path_bdio_w,["Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCut)]","$(BLINstr)Param$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jdl2"],OVERWRITE=OVERWRITE)
            pjld2FitRes = create_path(path_bdio_w,["Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCut)]","$(BLINstr)FitRes$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jld2"],OVERWRITE=OVERWRITE)
            pjld2Model  = create_path(path_bdio_w,["Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCut)]","$(BLINstr)ModelInfo$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jld2"],OVERWRITE=OVERWRITE)

            io = IOBuffer()
            
            extra = Dict{String, Any}("diag" => diag, "wind" => wind, "comp" => comp, "impr_set" => impr_set)

            # res:
            fb = ALPHAdobs_create(pBDIOres, io)
            
            ALPHAdobs_write(fb, res[impr_set], extra=extra)

            ALPHAdobs_close(fb)

            # xydata:
            fb = ALPHAdobs_create(pBDIOxydata, io)

            xDict = Dict{String,Array{uwreal}}("xdata" => xdata)
            ALPHAdobs_write(fb, xDict, extra=extra)

            yDict = Dict{String,Array{uwreal}}()
            for key in mykeys
                yDict[key] = ydata[impr_set][key]
            end
            ALPHAdobs_write(fb, yDict, extra=extra)

            ALPHAdobs_close(fb)

            # param:
            @save pjld2param Param=par[impr_set]

            # fit
            @save pjld2FitRes FitRes=fit[impr_set]

            # model info:
            infoDict = Dict{String,Any}(
                "length" => length(f_tot_isov),
                "nens" => length(ensListFit),
                "ensList" => getfield.(ensListFit,:id),
                "n_par_tot_isov" => n_par_tot_isov,
                "label_tot_isov" => label_tot_isov,
                "fit_cut" => FitCut,
                "3l_impr" => tl_IMPR,
                "PVALs" => PVAL,
                "fPiRescaling" => RESC,
                "simpleBase" => SimpleBase,
                "a2Rescaling" => a2RESC,
                "multTerm" => !isnothing(MultFunc),
                "na_max" => na_max,
                "nmPi_max" => nmPi_max,
                "nmK_max" => nmK_max
            )

            save(pjld2Model,"Info",infoDict)
        end
    end # end impr_set loop
end # end FitCut loop

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> READING TEST <==========================##

diag = "NLOb"  # LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
wind = "SDsub"  # NW  SD  SDsub  ID  LD  ILD
comp = "gCCconn"  # g33  g88  ∆ls_amu  ∆lc_b  gCCconn  gCCdisc  gC8disc

BLIND = false

impr_set = "2"

Q = 5.0  # virtuality for SDsub

model_var_list = [a3,a4,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log,phi4]
MultFunc = nothing  # nothing  deltaphi

FitCut = "None"  # "None"  "beta"  "mass"  "beta&mass"

STD_DERIV  = false
tl_IMPR    = false
VREF       = false
RESC       = false

SimpleBase = true

path_bdio = path_bdio_dict["clust"]

xData, yData = BDIOread_XYdata(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

fitRes = JDL2read_FitRes(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

modelInfo = JDL2read_ModelInfo(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

res, par = BDIOread_res(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q,param=true)


##



#---- local

## xData (local)
#  "A653" -> 0.0576667 +/- 0.000242656  0.777284 +/- 0.00625858    1.16593 +/- 0.00938787
#  "A654" -> 0.0570975 +/- 0.000226954  0.489888 +/- 0.00567642    1.14885 +/- 0.00789529
#  "C101" -> 0.0429432 +/- 5.89993e-5   0.213307 +/- 0.00288089    1.09194 +/- 0.00401301
#  "C102" -> 0.0436319 +/- 7.49727e-5   0.212957 +/- 0.00341234    1.19366 +/- 0.00481404
#  "D150" -> 0.0424275 +/- 5.75878e-5   0.0753275 +/- 0.00254158   1.06004 +/- 0.00383333
#  "D200" -> 0.0241456 +/- 2.78394e-5   0.175146 +/- 0.00131217    1.09904 +/- 0.00188295
#  "D201" -> 0.0243349 +/- 3.10494e-5   0.173452 +/- 0.00229556    1.17966 +/- 0.00333037
#  "D251" -> 0.024196 +/- 2.57496e-5    0.350046 +/- 0.00105397    1.10953 +/- 0.00161672
#  "D450" -> 0.0337908 +/- 2.43447e-5   0.205891 +/- 0.00306971    1.10376 +/- 0.00450609
#  "D451" -> 0.0341158 +/- 3.04913e-5   0.202449 +/- 0.00235469    1.20053 +/- 0.00205573
#  "D452" -> 0.033564 +/- 5.31664e-5    0.105769 +/- 0.00202771    1.09027 +/- 0.00226589
#  "E250" -> 0.0240327 +/- 1.13641e-5   0.0747481 +/- 0.000777181  1.09375 +/- 0.00141521
#  "E300" -> 0.0145048 +/- 8.96085e-6   0.133477 +/- 0.000663372   1.12647 +/- 0.00184096
#  "F300" -> 0.0144398 +/- 6.45839e-6   0.0797501 +/- 0.000824834  1.09613 +/- 0.00237204
#  "H101" -> 0.0439493 +/- 7.01528e-5   0.759044 +/- 0.00494632    1.13857 +/- 0.00741948
#  "H102" -> 0.0434346 +/- 7.90816e-5   0.546402 +/- 0.00423306    1.11722 +/- 0.00660999
#  "J303" -> 0.0145108 +/- 2.14579e-5   0.288569 +/- 0.001821      1.12952 +/- 0.00347268
#  "J304" -> 0.0147172 +/- 1.88298e-5   0.290567 +/- 0.00146282    1.32588 +/- 0.0030472
#  "J306" -> 0.0145499 +/- 2.28764e-5   0.519102 +/- 0.00272993    1.14172 +/- 0.00596861
#  "J307" -> 0.0145263 +/- 4.07208e-5   0.761211 +/- 0.00355714    1.14182 +/- 0.0053357
#  "J500" -> 0.00895109 +/- 1.41351e-5  0.742845 +/- 0.00373566    1.11427 +/- 0.00560348
#  "J501" -> 0.00897281 +/- 2.31645e-5  0.482809 +/- 0.00210724    1.1034 +/- 0.0037426
#  "N101" -> 0.0432093 +/- 7.10286e-5   0.340089 +/- 0.00326093    1.1089 +/- 0.00403096
#  "N200" -> 0.0242147 +/- 2.48423e-5   0.352457 +/- 0.0022516     1.11413 +/- 0.0037365
#  "N202" -> 0.024196 +/- 6.57362e-5    0.744882 +/- 0.00311755    1.11732 +/- 0.00467633
#  "N203" -> 0.0243158 +/- 2.49161e-5   0.520954 +/- 0.00296591    1.11324 +/- 0.00473157
#  "N451" -> 0.0338844 +/- 6.83456e-5   0.361444 +/- 0.00316985    1.11799 +/- 0.00393712
#  "N452" -> 0.0340155 +/- 8.39537e-5   0.539543 +/- 0.00197107    1.12242 +/- 0.00294464

## xData (clust)
#  "A653" -> 0.0576667 +/- 0.000242656  0.778861 +/- 0.00744844   1.16829 +/- 0.00894931
#  "A654" -> 0.0570975 +/- 0.000226954  0.48535 +/- 0.0073146     1.1461 +/- 0.00912858
#  "C101" -> 0.0429432 +/- 5.89993e-5   0.213225 +/- 0.00326647   1.09279 +/- 0.00393942
#  "C102" -> 0.0436319 +/- 7.49727e-5   0.214357 +/- 0.0034773    1.19249 +/- 0.0054006
#  "D150" -> 0.0424275 +/- 5.75878e-5   0.0753467 +/- 0.00250742  1.06082 +/- 0.00393178
#  "D200" -> 0.0241456 +/- 2.78394e-5   0.175357 +/- 0.00152263   1.09944 +/- 0.00243875
#  "D201" -> 0.0243349 +/- 3.10494e-5   0.173565 +/- 0.00230741   1.17979 +/- 0.00373272
#  "D251" -> 0.024196 +/- 2.57496e-5    0.356918 +/- 0.00128639   1.11346 +/- 0.00200189
#  "D450" -> 0.0337536 +/- 4.30779e-5   0.205526 +/- 0.00213828   1.10405 +/- 0.00263739
#  "D451" -> 0.0341158 +/- 3.04913e-5   0.204811 +/- 0.00148146   1.20582 +/- 0.00205851
#  "D452" -> 0.033564 +/- 5.31664e-5    0.105159 +/- 0.00195417   1.08899 +/- 0.0025889
#  "E250" -> 0.0240327 +/- 1.13641e-5   0.072355 +/- 0.00142322   1.0913 +/- 0.0014815
#  "E300" -> 0.0145048 +/- 8.96085e-6   0.133049 +/- 0.000972648  1.12181 +/- 0.00190402
#  "F300" -> 0.0144398 +/- 6.45839e-6   0.0791642 +/- 0.00107765  0.118746 +/- 0.00120537
#  "H101" -> 0.0439493 +/- 7.01528e-5   0.755095 +/- 0.00527924   1.13264 +/- 0.00602417
#  "H102" -> 0.0434346 +/- 7.90816e-5   0.545661 +/- 0.00513017   1.11661 +/- 0.00597601
#  "J303" -> 0.0145108 +/- 2.14579e-5   0.288927 +/- 0.00166224   1.13072 +/- 0.00322501
#  "J304" -> 0.0147172 +/- 1.88298e-5   0.291513 +/- 0.00181889   1.32735 +/- 0.00359966
#  "J306" -> 0.0145499 +/- 2.28764e-5   0.519015 +/- 0.00241182   1.14255 +/- 0.00364398
#  "J307" -> 0.0145263 +/- 4.07208e-5   0.765776 +/- 0.00646566   1.14866 +/- 0.00754083
#  "J500" -> 0.00895109 +/- 1.41351e-5  0.742607 +/- 0.00365445   1.11391 +/- 0.00425076
#  "J501" -> 0.00897281 +/- 2.31645e-5  0.482822 +/- 0.00359719   1.10329 +/- 0.00543648
#  "N101" -> 0.0432093 +/- 7.10286e-5   0.341645 +/- 0.00314365   1.11123 +/- 0.00375428
#  "N200" -> 0.0242147 +/- 2.48423e-5   0.352128 +/- 0.00239174   1.11407 +/- 0.0033016
#  "N202" -> 0.024196 +/- 6.57362e-5    0.744767 +/- 0.00399076   1.11715 +/- 0.00489917
#  "N203" -> 0.0243158 +/- 2.49161e-5   0.520866 +/- 0.00228478   1.11345 +/- 0.00285456
#  "N451" -> 0.0338844 +/- 6.83456e-5   0.361787 +/- 0.00203083   1.11848 +/- 0.00309411
#  "N452" -> 0.0340155 +/- 8.39537e-5   0.539442 +/- 0.00273528   1.12243 +/- 0.003985

## yData "2" g33_ll (local)
#  "A653" -> 144.66411773391522 +/- 1.140411432445603
#  "A654" -> 183.57257825156546 +/- 1.5680671964901316
#  "C101" -> 277.21557678208507 +/- 5.537709757572617
#  "C102" -> 275.0811211830782 +/- 6.424067440986098
#  "D150" -> 380.1868106773942 +/- 8.205934459784624
#  "D200" -> 285.8962077491209 +/- 3.517505801944449
#  "D201" -> 286.0715141613569 +/- 4.212666941004851
#  "D251" -> 222.7063261129124 +/- 1.2762884764750435
#  "D450" -> 271.1384717290882 +/- 1.9756985251922714
#  "D451" -> 273.8488718740623 +/- 2.6988369319299323
#  "D452" -> 336.602240794872 +/- 3.6946938613817193
#  "E250" -> 383.4041870352148 +/- 2.269968358848745
#  "E300" -> 309.7412474054965 +/- 4.846430036631955
#  "F300" -> 377.4980115385536 +/- 9.07708547760968
#  "H101" -> 154.34543041342354 +/- 1.4570931682007944
#  "H102" -> 181.6169455179249 +/- 1.9585790911316436
#  "J303" -> 233.92619031502352 +/- 3.470526076376498
#  "J304" -> 234.8002407701209 +/- 3.033352737407814
#  "J306" -> 187.63433849885934 +/- 2.756545145213178
#  "J307" -> 152.1937839366353 +/- 2.0301274865867134
#  "J500" -> 159.20163392340928 +/- 3.784136513682057
#  "J501" -> 188.9478926866317 +/- 2.703596775500886
#  "N101" -> 227.17012173838407 +/- 3.1860372893850757
#  "N200" -> 223.1745732717938 +/- 3.735169414254189
#  "N202" -> 155.21277294239215 +/- 2.5178503148785003
#  "N203" -> 187.7458167273421 +/- 2.4635771147848984
#  "N451" -> 222.12060465618495 +/- 2.4984739003383867
#  "N452" -> 180.98757661428564 +/- 4.205689342571214

## yData "2" g33_ll (clust)
#  "A653" -> 144.66411773391525 +/- 1.0892011048627726
#  "A654" -> 185.9817875925889 +/- 1.4998283541788435
#  "C101" -> 277.3165603651992 +/- 5.513293594144572
#  "C102" -> 275.14517206696405 +/- 6.395252895149104
#  "D150" -> 380.18702050664245 +/- 8.20453779490911
#  "D200" -> 285.9327137732262 +/- 3.5175053186143344
#  "D201" -> 286.09478979925007 +/- 4.19353856425146
#  "D251" -> 223.54155468519306 +/- 1.26522599193757
#  "D450" -> 270.77178579205514 +/- 1.908550332052629
#  "D451" -> 273.8980384000386 +/- 2.6904945762004937
#  "D452" -> 336.6041002948077 +/- 3.686324045628494
#  "E250" -> 383.4043363705271 +/- 2.2699683722255264
#  "E300" -> 309.74803715066867 +/- 4.795821468180395
#  "F300" -> 377.4982481753033 +/- 8.893502721186255
#  "H101" -> 154.34543041342354 +/- 1.4508335832922723
#  "H102" -> 186.02176451081186 +/- 1.951600860627537
#  "J303" -> 234.27526754407648 +/- 3.470607051319889
#  "J304" -> 234.97819054717087 +/- 2.996277781505724
#  "J306" -> 191.01122849997603 +/- 2.7100367938143597
#  "J307" -> 152.1937839366353 +/- 2.0225428935214307
#  "J500" -> 159.20163392340928 +/- 3.783728668163358
#  "J501" -> 191.95888242727824 +/- 2.7100315399551222
#  "N101" -> 227.91743896916572 +/- 3.1529877643289344
#  "N200" -> 224.02305435620428 +/- 3.7313771935484996
#  "N202" -> 155.21277294239215 +/- 2.51034329240025
#  "N203" -> 191.50430694893734 +/- 2.448805583218909
#  "N451" -> 223.0464777865849 +/- 2.4720931336814735
#  "N452" -> 185.13092926880782 +/- 3.996543803919973


