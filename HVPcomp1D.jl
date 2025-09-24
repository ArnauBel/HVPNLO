# Import packages

using Revise

include("HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

using ProgressBars
using Suppressor
using TimerOutputs

# include uwreal constants

# include("HVPtool/uwConst.jl")

# BDIO path definition (set 'STD_DERIV = true' to use the standard sym. derivative in the impr.)

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_tl   = joinpath(julia_script_directory, "..", "HVPData", "tree_level_improv")
path_spec = joinpath(julia_script_directory, "..", "HVPData", "spectroscopy")

pFVC_MLL  = joinpath(julia_script_directory, "..", "HVPData", "FSE_MLL")

# Blind analysis (Simon K.) safe ensembles: 
# SU(3) sym.  H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

# All considered ensembles are:
ensList = [
    "A653","A654","B450","C101","C102",
    "D150","D200","D201","D251","D450",
    "D451","D452","E250","E300","F300",
    "H101","H200","J303","J304","J306",
    "J307","J500","J501","N101","N200",
    "N202","N203","N302","N451","N452",
    "S400"] # "H102","N300"


ensInfo = EnsInfo.(ensList)

# We do not have charm or disconnected data for some of the ensembles
ensNOcharm  = ["C102","D150","D201","D251","D451","F300","H200","J304","J306","J307","J501","N451","N452"]
ensNOdisc   = ["F300","J306"]

ensSPECdata = ["D200","E250"]  # J303

# ensTAILfit  = ["A653","A654","B450"]

wpm = Dict{String, Vector{Float64}}()

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> 1D HVP computation (+ BM) [LO, NLOa, NLOb] <==========================##

diag = ""  # LO  NLOa  NLOb  NLOa&b
wind = ""  # NW  SD  ID  LD  ILD

IMPR_SET = ["1","2"] # ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

STD_DERIV = false
RESC      = false

OVERWRITE = false  # Allows to erase data and overwrite it with new data, use carefully !!

BLIND = false

path_bdio = path_bdio_dict["local"]


BM = wind in ["NW","LD","ILD"]


@info("STARTING HVP COMPUTATION [diag. $diag; wind. $wind; resc. $RESC]")
STD_DERIV ? @info("STANDARD DERIVATIVE is being employed in the IMPROVEMENT") : nothing

# LO => (alpha/pi)^2 || NLO => (alpha/pi)^3
exp_diag = diag == "LO" ? 2 : 3

if wind == "SD"
    corr33tl_ll, corr33tl_lc = read_tree_level_v33(path_tl, cons=true)
    corr33tl_v3s03_ll, corr33tl_v3s03_lc = read_tree_level_v3sig03(path_tl, cons=true, massless=true)

    corr33tl_ll = uwreal.(corr33tl_ll); corr33tl_lc = uwreal.(corr33tl_lc)
end


@time begin
    for ens in ensInfo

        @info("Computing HVP for ensemble $(ens.id)")
        ens.id ∈ ensNOcharm ? @info("  > NO CHARM DATA FOR $(ens.id)") : nothing
        ens.id ∈ ensNOdisc  ? @info("  > NO DISC. DATA FOR $(ens.id)") : nothing

        println("   - Reading TMR...")
        TMR = BDIOread_TMR(path_bdio,ens,diag,resc=RESC,beta=false,BLIND=BLIND)
        ((ens.id ∉ ensNOcharm) || (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)) ? TMRbeta = BDIOread_TMR(path_bdio,ens,diag,resc=RESC,beta=true) : nothing     
        

        if ens.id in ensSPECdata
            println("   - Reading spectral data...")
            E, Z, Z_impr = get_spectr_data(path_spec,ens)
        end

        T = HVPobs.Data.get_T(ens.id)
        sym_points = Int64(T/2+1)
        t = collect(1:sym_points); tBM = collect(1:length(TMR))
        if !RESC
            println("   - Reading t0...")
            t0 = BDIOread_t0(path_bdio,ens)
            aens = sqrtt0_ph / sqrt(t0); tfm = aens.*(collect(1:length(TMR)).-1)
            if (ens.id ∉ ensNOcharm) || (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)
                a_ß = sqrtt0_ph / sqrt(t0sym(ens.beta)); tfm_ß = a_ß.*(collect(1:length(TMRbeta)).-1)
            end
        else
            println("   - Reading fPi...")
            fPi = BDIOread_fPS(path_bdio,ens)["fPi"]
            aens = hbarc * fPi / fPi_ph; tfm = aens.*(collect(1:length(TMR)).-1)
            if (ens.id ∉ ensNOcharm) || (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)
                a_ß = hbarc * fPiph(ens.beta) / fPi_ph; tfm_ß = a_ß.*(collect(1:length(TMRbeta)).-1)
            end
        end

        # define windowed kernel
        if wind == "NW"
            TMRw  = TMR
            ((ens.id ∉ ensNOcharm) || (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)) ? TMRwbeta = TMRbeta : nothing
        else
            TMRw  = TMR .* Window(wind)(tfm)
            ((ens.id ∉ ensNOcharm) || (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)) ? TMRwbeta = TMRbeta .* Window(wind)(tfm_ß) : nothing
        end

        for impr_set in IMPR_SET

            println("   - Starting set "*impr_set)

            println("      - Reading corr...")
            corr = BDIOread_corr(path_bdio,ens,impr_set,STD=STD_DERIV)

            if ens.kappa_l == ens.kappa_s || ens.id ∈ ensNOdisc
                println("      - SU(3) flavour sym point or no disc. data available")

                light_keys = ["g33_ll","g33_lc"]
            else
                corr["g88_ll"] = corr["g88conn_ll"] .+ corr["g88disc_ll"] .+ (2).*corr["g08conn_ll"] .+ corr["g08disc_ll"] .+ corr["g80disc_ll"]
                corr["g88_lc"] = corr["g88conn_lc"] .+ corr["g88disc_lc"] .+ corr["g08conn_lc"] .+ corr["g08disc_lc"]

                light_keys = ["g33_ll","g33_lc","g88conn_ll","g88conn_lc","g88_ll","g88_lc"]
            end

            if wind == "SD"
                corr["g33tl_ll"] = corr33tl_ll[:]; corr["g33tl_lc"] = corr33tl_lc[:]

                beta = ens.beta
                if impr_set == "1"
                    cv_l = cv_loc(beta)
                    cv_c = cv_cons(beta)
                elseif impr_set =="1old"
                    cv_l = cv_loc_old(beta)
                    cv_c = cv_cons_old(beta)
                elseif impr_set =="2"
                    cv_l = cv_loc_set2(beta)
                    cv_c = cv_cons_set2(beta)
                end

                improve_corr_vkvk!(corr["g33tl_ll"], corr33tl_v3s03_ll, 2*cv_l, std=false, treelevel=true)
                improve_corr_vkvk_cons!(corr["g33tl_lc"], corr33tl_v3s03_ll, corr33tl_v3s03_lc, cv_l, cv_c, std=STD_DERIV, treelevel=true)

                corr["g33tl_ll"] ./= 2; corr["g33tl_lc"] ./= 2

                light_keys = union(light_keys,["g33tl_ll","g33tl_lc"])
            end

            println("      - Computing HVPs...")

            HVP = Dict{String, uwreal}()
            HVPsyst = Dict{String, Float64}()

            if BM
                PlatRec = Dict{String, Vector{Float64}}()

                tcut0 = 10
                tstep = 1
                tEeffFix = 1.2  # in  fm; point where the eff energy of the lower bound is fixed 

                uwerr(E0_ens[ens.id]["E0"]); E0 = E0_ens[ens.id]["E0"].mean + E0_ens[ens.id]["E0"].err

                for key in light_keys
                    
                    ub = Vector{uwreal}()
                    lb = Vector{uwreal}()

                    obs = corr[key]
                    int = obs[t] .* TMRw[t]

                    Eeff_ = uwreal(0.0)

                    if ens.id in ensSPECdata && key[1:3] == "g33"
                        println("         - $key (Rec w. spectroscopy data...)")

                        corrVec = reconstr_corr(ens,E,Z,Z_impr,impr_set=impr_set,nmax=Dict("E250"=>4,"D200"=>2)[ens.id],IMPR=true,RENORM=true,total=false)
                        # corrVec = [corr[t] for corr in corrVec]
                        uwerr.(corrVec[end]); uwerr.(obs[t])
                        trec = findfirst_uninterrupted(err.(corrVec[end]) .< err.(obs[t]))

                        if key[end-1:end] == "lc"
                            obsRec = corrVec[end] .* (corr["g33_lc"][t]./corr["g33_ll"][t])
                        elseif key[end-1:end] == "ll"
                            obsRec = corrVec[end]
                        end

                        intRec = obsRec[trec:end] .* TMRw[t][trec:end]

                        HVP[key]     = (alpha/pi)^exp_diag * (sum(int[1:trec-1])+sum(intRec)) * 1e10
                        HVPsyst[key] = 0.0
                        PlatRec[key] = [value(aens) * trec]
                    # elseif ens.id in ensTAILfit && key[1:3] == "g33"
                    #     println("         - $key (Rec w. cosh fit...)")
                        
                    #     @. exp_model(x0,p)   = p[2] * exp(-p[1] * x0)
                    #     @. cosh_model(x0,p)  = p[2] * (exp(-p[1] * x0) + exp(-p[1] * (T-x0)))

                    #     p0_tuple = p0_smallpbc_dict["$(ens.id)_set$(impr_set)_$(key)"]
                    #     fit_vec = []
                    #     for p0 in p0_tuple
                    #         data = obs[t][p0:end] 
                    #         fit  = fit_routine(cosh_model,collect(p0:sym_points).-1, data, 2, pval=true, info=false, lineprint=false)
                    #         push!(fit_vec,fit)
                    #     end
                        
                    #     w = get_w_from_fitres(vcat(fit_vec...), AIC=true)

                    #     param = getfield.(vcat(fit_vec...),:param)
                    #     p1_vec = [par[1] for par in param]
                    #     p2_vec = [par[2] for par in param]

                    #     p1_res, p1_sys = model_average(p1_vec, w)
                    #     p2_res, p2_sys = model_average(p2_vec, w)

                    #     p1 = p1_res[1] + uwreal([0.0,p1_sys],"p1 fit $(ens.id)"); uwerr(p1)
                    #     p2 = p2_res[1] + uwreal([0.0,p2_sys],"p2 fit $(ens.id)"); uwerr(p2)

                    #     trec = p0_tuple[1] + (argmax(w)-1)

                    #     obsRec = exp_model(t[trec:end].-1,[p1,p2])

                    #     intRec = obsRec .* TMRw[trec:end]

                    #     HVP[key]     = (alpha/pi)^exp_diag * (sum(int[1:trec-1])+sum(intRec)) * 1e10
                    #     HVPsyst[key] = 0.0
                    #     PlatRec[key] = [value(aens) * trec]
                    else
                        println("         - $key (BM...)")

                        for tcut in tcut0:tstep:t[end-1]  # compute the upper and lower corr bounds 

                            if aens*tcut < tEeffFix  # we fix the eff energy at some point
                                Eeff_=Eeff(tcut, obs)
                            end

                            UB = corr_bound(tBM, tcut, obs, E0)
                            LB = corr_bound(tBM, tcut, obs, Eeff_)
                    
                            UBint = UB .* TMRw[tcut+1:end]
                            LBint = LB .* TMRw[tcut+1:end]

                            upper_int = vcat(int[1:tcut],UBint...)
                            lower_int = vcat(int[1:tcut],LBint...)
                            
                            if diag == "NLOa&b" && aens*tBM[end] > NLOab_zero
                                # since there's a sign flip in ; the roles of the upper und lower bound is reversed for larger time slices
                                
                                bit_vec = value.(aens .* (tBM.-1)) .< NLOab_zero

                                upper_int, lower_int = [vcat(upper_int[bit_vec],lower_int[.!bit_vec]...),vcat(lower_int[bit_vec],upper_int[.!bit_vec]...)]
                            end

                            ub_ = (alpha/pi)^exp_diag * sum(upper_int) * 1e10
                            lb_ = (alpha/pi)^exp_diag * sum(lower_int) * 1e10

                            push!(ub, ub_)
                            push!(lb, lb_)
                        end

                        HVP[key], HVPsyst[key], PlatRec[key] = bounding_method(ub,lb,aens,PLAT=true,AVER=false,tcut0=tcut0)
                    end
                end

                if ens.kappa_l == ens.kappa_s
                    HVP["g88conn_ll"] = HVP["g88_ll"] = HVP["g33_ll"]
                    HVP["g88conn_lc"] = HVP["g88_lc"] = HVP["g33_lc"]

                    HVPsyst["g88conn_ll"] = HVPsyst["g88_ll"] = HVPsyst["g33_ll"]
                    HVPsyst["g88conn_lc"] = HVPsyst["g88_lc"] = HVPsyst["g33_lc"]

                    PlatRec["g88conn_ll"] = PlatRec["g88_ll"] = PlatRec["g33_ll"]
                    PlatRec["g88conn_lc"] = PlatRec["g88_lc"] = PlatRec["g33_lc"]
                end
            else
                for key in light_keys
                    println("         - $key")

                    obs = corr[key][t]

                    int = obs .* TMRw[t]
                    amu = (alpha/pi)^exp_diag * sum(int) * 1e10

                    HVP[key] = amu
                end

                if ens.kappa_l == ens.kappa_s
                    HVP["g88conn_ll"] = HVP["g88_ll"] = HVP["g33_ll"]
                    HVP["g88conn_lc"] = HVP["g88_lc"] = HVP["g33_lc"]
                end
            end

            if ens.id ∉ ensNOcharm
                for key in ["gCCconn_ll_sim","gCCconn_lc_sim","gCCconn_ll_sim+","gCCconn_lc_sim+"]
                    println("         - $key  &  $(key[1:8]*"SU3"*key[8:end])")

                    obs = corr[key][t]

                    int      = obs .* TMRw[t]
                    int_SU3  = obs .* TMRwbeta[t]
                    amu      = (alpha/pi)^exp_diag * sum(int) * 1e10
                    amu_SU3  = (alpha/pi)^exp_diag * sum(int_SU3) * 1e10

                    HVP[key] = amu
                    HVP[key[1:8]*"SU3"*key[8:end]] = amu_SU3
                end

                println("         - Reading kappaC target...")

                # kappaC_tar = BDIOread_KappaC_tar(path_bdio,ens)
                mDs_ph_prime, mDs_beta, Ds_dict = BDIOread_mDs_kappaC(path_bdio,ens)
                kappaC_tar = Ds_dict["kappaC"]

                println("         - Interpolating charm contributions...")

                # kappaC_tar = uwreal([kcd_in[ens.id]["kappaC"],kcd_in[ens.id]["kappaC_err"]], "kappaC target")   # extract KappaC from tables
                kappaC = [kcd_in[ens.id]["kappaC_sim"],kcd_in[ens.id]["kappaC_sim_plus"]]
                
                DkappaC = abs(kappaC[1]-kappaC[2])
                DkappaC_tar = minimum([abs(value(kappaC_tar)-kappaC[1]),abs(value(kappaC_tar)-kappaC[2])])

                for key in ["gCCconn_ll","gCCconn_lc","gCCconn_SU3_ll","gCCconn_SU3_lc"]
                    obs  = [HVP["$(key)_sim"],HVP["$(key)_sim+"]]
                    
                    par, _ = lin_fit(kappaC,obs,wpm=wpm,lineprint=false)
                    obs_tar = y_lin_fit(par,kappaC_tar)

                    if kappaC[1] <= value(kappaC_tar) <= kappaC[2]
                        HVP[key] = obs_tar
                        HVPsyst[key] = 0.0
                    else
                        HVP[key] = obs_tar # + 1.5 * (abs(DkappaC_tar) / DkappaC)^2 * uwreal([0.0,err(obs_tar)],"kappaC syst")
                        uwerr(obs_tar); HVPsyst[key] = 1.5 * (abs(DkappaC_tar) / DkappaC)^2 * err(obs_tar)
                    end
                end
            end

            if wind == "SD" && ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc
                for (k,key) in enumerate(["gCCdisc_cc","gC8disc_cc"])
                    println("         - $key")

                    obs = corr[key][t]

                    int = k==1 ? (obs .* TMRwbeta[t]) : (obs .* TMRw[t])
                    amu = (alpha/pi)^exp_diag * sum(int) * 1e10

                    HVP[key] = amu
                end
            end

            println("      - Writing BDIO & JDL2...")

            DERstr = STD_DERIV ? "_std" : "" 
            BLstr  = BLIND ? "Blind" : ""
            RESstr = RESC ? "_resc" : ""
            pBDIO  = create_path(path_bdio,["HVP&FVC",wind,ens.id,"$(ens.id)_$(BLstr)HVP$(diag)_set$(impr_set)$(RESstr)$(DERstr)"],OVERWRITE=OVERWRITE)
            pjdl2  = create_path(path_bdio,["HVP&FVC",wind,ens.id,"$(ens.id)_$(BLstr)HVP$(diag)_info_set$(impr_set)$(RESstr)$(DERstr).jld2"],OVERWRITE=OVERWRITE)

            io = IOBuffer()
            write(io, "$(ens.id) HVP")
            fb = ALPHAdobs_create(pBDIO, io)

            extra = Dict{String, Any}("ens" => ens.id, "impr_set" => impr_set, "diag" => diag, "wind" => wind)
            ALPHAdobs_write(fb, HVP, extra=extra)

            ALPHAdobs_close(fb)

            if BM
                HVPinfo = Dict{String,Dict}(
                    "HVPsyst" => HVPsyst,
                    "plateau/Reconstr" => PlatRec,
                )
            else
                HVPinfo = Dict{String,Dict}(
                    "HVPsyst" => HVPsyst,
                )            
            end
            save(pjdl2,"HVPinfo",HVPinfo)
        end # end impr_set loop
    end # end ens loop
end # end timer

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> 1D FVC computation [LO NLOa NLOb] <==========================##

diag = ""  # LO  NLOa  NLOb  NLOa&b
wind = ""  # NW  SD  ID  LD  ILD

RESC  = false
VREF  = false

BLIND = false

OVERWRITE = false  # Allows to erase data and overwrite it with new data, use carefully !!

path_bdio = path_bdio_dict["local"]

if !VREF
    @info("STARTING FVC COMPUTATION [diag. $diag; wind. $wind; resc. $RESC]")
else
    @info("STARTING FVC COMPUTATION (to Vref) [diag. $diag; wind. $wind; resc. $RESC]")
end

@time begin
    for ens in ensInfo

        @info("Computing for ensemble $(ens.id)")

        println("   - Reading TMR...")
        TMR = BDIOread_TMR(path_bdio,ens,diag,resc=RESC,beta=false,BLIND=BLIND)

        # define windowed kernel
        if wind in ["NW","LD","ILD"]
            t = collect(1:length(TMR))
        else
            T = HVPobs.Data.get_T(ens.id)
            t = collect(1:Int64(T/2+1))
            TMR = TMR[t]
        end

        if !RESC
            println("   - Reading t0...")
            t0 = BDIOread_t0(path_bdio,ens)
            aens = sqrtt0_ph / sqrt(t0); tfm = aens.*(t.-1)
        else
            println("   - Reading fPi...")
            fPi = BDIOread_fPS(path_bdio,ens)["fPi"]
            aens = hbarc * fPi / fPi_ph; tfm = aens.*(t.-1)
        end

        if wind == "NW" 
            TMRw = TMR
        else
            TMRw = TMR .* Window(wind)(tfm)
        end

        println("   - Reading fvc (HP)...")
        fvc_hp = BDIOread_FVCcorr(path_bdio,ens,Vref=VREF)

        if ens.id != "A653"
            println("   - Reading fvc (MLL)...")
            t_mll, fvc_mll = TXTread_FVCcorr_GS(pFVC_MLL,ens,Vref=VREF)
        end
        
        # LO => (alpha/pi)^2 || NLO => (alpha/pi)^3
        exp_diag = diag == "LO" ? 2 : 3

        println("   - Computing ∆aµ...")

        # tFVC = 1:Int64(length(fvc_hp["FVCPi1"])/2)
        Tover2 = Int64(HVPobs.Data.get_T(ens.id)/2+1)

        fvcPi_hp = vcat(TMRw[1],-fvc_hp["FVCPi6"])[t]
        fvcPi = ens.id != "A653" ? vcat(fvcPi_hp[1:Int64(t_mll[1])],-fvc_mll...)[t] : fvcPi_hp
        
        fvcK  = vcat(TMRw[1],-fvc_hp["FVCK6" ])[t]

        # fvcPi_hp = vcat(-fvc_hp["FVCPi6"][1:Tover2])
        # fvcPi = ens.id != "A653" ? vcat(fvcPi_hp[1:Int64(t_mll[1])],-fvc_mll...)[1:Tover2] : fvcPi_hp
        
        # fvcK  = vcat(-fvc_hp["FVCK6" ][1:Tover2])

        FVCPiint_hp = fvcPi_hp .* TMRw
        FVCPiint = fvcPi .* TMRw
        FVCKint  = fvcK  .* TMRw

        FVCPi_hp = (alpha/pi)^exp_diag * sum(FVCPiint_hp) * 1e10
        FVCPi = (alpha/pi)^exp_diag * sum(FVCPiint) * 1e10
        FVCK  = (alpha/pi)^exp_diag * sum(FVCKint ) * 1e10

        if ens.kappa_l == ens.kappa_s
            FVC33_hp = (3/2) * FVCPi_hp
            FVC33 = (3/2) * FVCPi
            FVC88 = (3/2) * FVCPi
        else
            FVC33_hp = FVCPi_hp + FVCK/2
            FVC33 = FVCPi + FVCK/2
            FVC88 = 2/3 * FVCK
        end


        FVC =  Dict{String, uwreal}(
            "FVCPi_hp" => FVCPi_hp, 
            "FVCPi" => FVCPi, 
            "FVCK"  => FVCK/2,
            "FVCg33_hp" => FVC33_hp, 
            "FVCg33" => FVC33, 
            "FVCg88" => FVC88
        )

        println("   - Writing BDIO...")

        REFstr = VREF ? "_Vref" : ""
        RESstr = RESC ? "_resc" : ""
        BLstr  = BLIND ? "Blind" : ""
        pBDIO = create_path(path_bdio,["HVP&FVC",wind,ens.id,"$(ens.id)_$(BLstr)FVC$(diag)$(RESstr)$(REFstr)"],OVERWRITE=OVERWRITE)

        io = IOBuffer()
        write(io, "$(ens.id) FVC")
        fb = ALPHAdobs_create(pBDIO, io)

        extra = Dict{String, Any}("ens" => ens.id, "diag" => diag, "wind" => wind, "Vref" => VREF)
        ALPHAdobs_write(fb, FVC, extra=extra)

        ALPHAdobs_close(fb)
    end
end # end timer


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##


##==========================> READING TEST <==========================##

diag     = "NLOa&b"  #  "LO"  "NLOa"  "NLOb"  "NLOc"
wind     = "ID"  #  "NW"  "SD"  "ID"  "LD"  "ILD"
ensid    = "H101"
impr_set = "2"

STD   = false
VREF  = true
RESC  = false

BLIND = false

path_bdio = path_bdio_dict["local"]

HVP, info = BDIOread_HVPens(path_bdio,diag,wind,ensid,impr_set,info=true,resc=RESC,STD=STD,BLIND=BLIND)

FVC = BDIOread_FVCens(path_bdio,diag,wind,ensid,Vref=VREF,resc=RESC,BLIND=BLIND)

##

for ens in ensInfo
    HVP, info = BDIOread_HVPens(path_bdio,diag,wind,ens.id,impr_set,info=true,resc=RESC,STD=STD,BLIND=BLIND)
    println("$(ens.id): $(print_uwreal(HVP["g33_ll"])) - $(print_uwreal(HVP["g33_lc"]))")
end

##==> Result computation in the 'flavour basis'

ensid    = "A653"
impr_set = "1"
discr    = "ll"

HVP = BDIOread_HVPens(path_bdio,diag,wind,ensid,impr_set,info=false)

light   = 5/9*(2*HVP["g33_$(discr)"]); uwerr(light)
strange = 1/9*(3*HVP["g88conn_$(discr)"]-HVP["g33_$(discr)"]); uwerr(strange)
charm   = 4/9*HVP["gCCconn_$(discr)"]; uwerr(charm)

println("Ensemble: $ensid [impr. set $impr_set & discr. $discr]")
println("- amu(l) = $light")
println("- amu(s) = $strange")
println("- amu(c) = $charm")

##

diag = ""  # LO  NLOa  NLOb  NLOa&b
wind = ""  # NW  SD  ID  LD  ILD

IMPR_SET = ["1","2"] # ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

STD_DERIV = false

OVERWRITE = false  # Allows to erase data and overwrite it with new data, use carefully !!

BLIND = false

path_bdio = path_bdio_dict["local"]

for ens in ensInfo

    HVP, info = BDIOread_HVPens(path_bdio,diag,wind,ensid,impr_set,info=true,BLIND=BLIND)

    DERstr = STD_DERIV ? "_std" : "" 
    BLstr  = BLIND ? "Blind" : ""
    pBDIO  = create_path(path_bdio,["HVP&FVC",wind,ens.id,"$(ens.id)_$(BLstr)HVP$(diag)_set$(impr_set)$(DERstr)"],OVERWRITE=true)
    pjdl2  = create_path(path_bdio,["HVP&FVC",wind,ens.id,"$(ens.id)_$(BLstr)HVP$(diag)_info_set$(impr_set)$(DERstr).jld2"],OVERWRITE=true)

    if BM
        HVPinfo = Dict{String,Dict}(
            "HVPsyst" => HVPsyst,
            "plateau/Reconstr" => PlatRec,
        )
    else
        HVPinfo = Dict{String,Dict}(
            "HVPsyst" => HVPsyst,
        )            
    end

    io = IOBuffer()
    write(io, "$(ens.id) HVP")
    fb = ALPHAdobs_create(pBDIO, io)

    extra = Dict{String, Any}("ens" => ens.id, "impr_set" => impr_set, "diag" => diag, "wind" => wind, "HVPinfo" => HVPinfo)
    ALPHAdobs_write(fb, HVP, extra=extra)

    ALPHAdobs_close(fb)
end
