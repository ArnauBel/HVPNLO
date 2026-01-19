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

ensList = ["A653","A654","B450","C101","C102","D150","D200","D201","D251","D450","D451","D452","E250","E300","F300","H101","H102","H200","H650","J303","J304","J306","J307","J500","J501","N101","N200","N202","N203","N300","N302","N451","N452","S400"]

ensInfo = EnsInfo.(ensList)

# We do not have charm or disconnected data for some of the ensembles
ensNOcharm = ["C102","D150","D201","D251","D451","F300","H200","H650","J304","J306","J307","J501","N451","N452"]
ensNOdisc  = ["F300","J306"]

# Dict to relate contributions to their Dict keys

DictComptoKey = Dict{String,Vector{String}}(
    "g33"      => ["g33_ll","g33_lc"],
    "g88"      => ["g88_ll","g88_lc"],
    "g88conn"  => ["g88conn_ll","g88conn_lc"],
    "gSS"      => ["gSS_ll","gSS_lc"],
    # only interested in the lc (local-conserved) discr. for the cc conn
    "gCCconn"  => ["gCCconn_SU3_lc"], # ["gCCconn_ll","gCCconn_lc"]
    # "gCCconn"  => ["gCCconn_SU3_ll","gCCconn_SU3_lc"], # ["gCCconn_SU3_ll","gCCconn_SU3_lc"]

    "∆ls_amu" => ["∆ls_amu_ll","∆ls_amu_lc"],
    "∆ls_amuconn" => ["∆ls_amuconn_ll","∆ls_amuconn_lc"],
    "∆lc_b"   => ["∆lc_b_lc"],

    # only interested in the cc (conserved-conserved) discr. for the cc disc  & c8 disc
    "gCCdisc" => ["gCCdisc_cc"],
    "gC8disc" => ["gC8disc_cc"],

    "g3333"    => ["g3333_llll","g3333_lclc"],
    "g8888"    => ["g8888_llll","g8888_lclc"],
    "gCCCC"    => ["gCCCC_lclc"],
    "g3388"    => ["g3388_llll","g3388_lclc"],
    "g33CC"    => ["g33CC_llll","g33CC_lclc"],
    "g88CC"    => ["g88CC_llll","g88CC_lclc"]
)

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> Fits [LO, NLOa&b, NLOa, NLOb, NLOc] <==========================##

diag = "NLOb"  #  LO  NLOa  NLOb  NLOc  NLOa&b
wind = "SDsub"  #  NW  SD  SDsub  ID  ILD  LD  LD1  LD2

readIMPR_SET = ["1","2"] #  ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

BLIND = false

STD_DERIV  = false
tl_IMPR    = true
VREF       = true
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
    
    # if !VREF; apply_syst_FVC!(FVC[ens.id],diag,wind,ens.id,IMPR_SET=readIMPR_SET,factor=0.1); end
    apply_syst_FVC!(FVC[ens.id],diag,wind,ens.id,IMPR_SET=readIMPR_SET,factor=0.1)

    println("   - aµ = HVP + FVC")

    HVP_VolCorrect!(HVP[ens.id],FVC[ens.id],diag,IMPR_SET=readIMPR_SET)

    if wind in ["SD","SDsub"] && tl_IMPR
        println("   - Applying tree-level improvement...")

        HVP_3limpr!(HVP[ens.id],amu3l,IMPR_SET=readIMPR_SET,meth="prod")
    end
end

@info(" ⟹ Data ready to be fitted\n\n")


##==========================> Data ready to fit

comp = "g33"  #  g33  g88  gSS  gCCconn  strange gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆ls_amuconn  ∆lc_b

Q = 5.0  # virtuality for SDsub

model_var_list = Function[a3,a4,a2phi2,a2phi4,phi2sqr,phi2log]  #  [a3,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log,phi2inv,logphi2]  [a3,a2phi2,a2phi4,a3phi2,phi2sqr,phi2log]  [a3,a2phi2,phi2sqr,phi2log,phi4]
# model_var_list = Function[a3,a2y,ysqr,ylog]  #  [a3,a4,a2y,a2z,a3y,ysqr,ylog]  [a3,a2y,ysqr,ylog]

MultFunc = nothing  #  nothing  deltaphi

IMPR_SET = [readIMPR_SET[2]]  #  readIMPR_SET  ["1"]  ["2"]  ["1","2"]  ["1old","2"]  ["1","1old","2"]

FITCUT = ["beta&mass"]  #  ["None","beta","mass","beta&mass","beta_ext"]  ["None","beta","mass","beta&mass"]  ["None","beta"]

SimpleBase = false
a2RESC     = false

FitINFO    = true

PVAL       = false

WRITE      = false
OVERWRITE  = false

mdof = 4  # minimum number of d.o.f. allowed

mykeys = [DictComptoKey[comp][1]]  #  [DictComptoKey[comp][1]]

path_bdio_w = path_bdio_dict["local"]

# Following the LD paper, when it comes to the iso-vector analysis, the 'untrusted' ensembles are: H105, H200, N300,  N302, S400
ensExcl = ["H105","H200","A654","N300","N302","S400"] # 33, 88, SS, ∆ls(aµ) (D201 problems for SD & ID)
# ensExcl = ["H105","H200","A654","N300","N302","S400","B450","A653"] # g3333, g8888, g3388
# ensExcl = ["H105","H200","S400"] # ∆lc(b)
# ensExcl = ["H105","H200"] # CC


if comp != "g33" && VREF
    error("Cannot project to Vref for chosen iso-spin")
end


for FitCut in FITCUT

    # aDP      = FitCut in ["beta","beta&mass"] ? 5 : 6
    # aDOF     = comp in ["gCCconn"] ? 1 : 2
    # na_max   = aDP - aDOF
    # nmPi_max = ((FitCut in ["mass","beta&mass"] || comp == "gCCconn") ? 1 : 2)
    # nmK_max  = 1
    # # add a limit to "pure" phi2 terms
    # if wind in ["LD","ILD"]
    #     nmPi_max = [nmPi_max,1]
    # end

    # following SD, ID, LD papers :
    na_max   = 2
    nmPi_max = wind in ["NW","ILD","LD","LD1","LD2"] ? 2 : 1
    nmK_max  = 1


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
    elseif FitCut == "beta_ext"
        argExcl1 = getfield.(ensInfo,:beta) .== 3.34
        argExcl2 = getfield.(ensInfo,:beta) .== 3.4
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
        if ((comp ∈ ["gCCconn","∆lc_b","g33CC","g88CC","gCCCC"]) && (ens.id ∈ ensNOcharm)) || (comp ∈ ["g88","∆ls_amu","g3388","g8888","g88CC"] && ens.id ∈ ensNOdisc) || (comp ∈ ["gCCdisc","gC8disc"] && (ens.id ∈ ensNOdisc || ens.kappa_l == ens.kappa_s)) || (comp in ["∆ls_amu","∆ls_amuconn"] && ens.kappa_l == ens.kappa_s && MultFunc == deltaphi) || (ens.id ∈ union(ensExcl,ensExcl_type))
            println("   - [$(ens.id): Either excluded or 0 for contribution $comp]")
            continue
        end
        # println("   - [$(ens.id): Accepted]")
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

diag = "NLOa&b"  # LO  NLOa  NLOb  NLOc  NLOa&b
wind = "SDsub"  # NW  SD  SDsub  ID  LD  ILD
comp = "∆ls_amu"  # g33  g88  ∆ls_amu  ∆lc_b  gCCconn  gCCdisc  gC8disc

BLIND = false

impr_set = "1"

Q = 5.0  # virtuality for SDsub

model_var_list = [a3,a2phi2,phi2sqr,phi2log]
MultFunc = deltaphi  # nothing  deltaphi

FitCut = "None"  # "None"  "beta"  "mass"  "beta&mass"

STD_DERIV  = false
tl_IMPR    = false
VREF       = false
RESC       = false

SimpleBase = false

path_bdio = path_bdio_dict["local"]

xData, yData = BDIOread_XYdata(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

fitRes = JDL2read_FitRes(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

modelInfo = JDL2read_ModelInfo(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

res, par = BDIOread_res(path_bdio,diag,wind,comp,model_var_list,FitCut,impr_set;resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q,param=true)

##

# impr_set = ""

# println("Ens\t HVP\t\t FVC")
# println("-------------------------------------")
# for key in sort!(collect(keys(FVC)))
#     println("$key \t $(print_uwreal(HVP[key][impr_set]["g33_ll"])) \t $(print_uwreal(FVC[key]["FVCg33"]))")
# end
