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
    "S400"
    ] # "H102","N300"
ensList = ["E250","E300","F300","J303","J304","J306","J307","J500","J501"]
# E250 & D200 have spectroscopy data !!

ensInfo = EnsInfo.(ensList)

# We do not have charm or disconnected data for some of the ensembles
ensNOcharm  = ["C102","D150","D201","D251","D451","F300","H200","J304","J306","J307","J501","N451","N452"]
ensNOdisc   = ["F300","J306"]

ensSPECdata = ["D200","E250"]  # J303

# ensTAILfit  = ["A653","A654","B450"]

wpm = Dict{String, Vector{Float64}}()

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> 2D HVP computation (+ BM) [NLOc] <==========================##

# diag = "NLOc"
# wind = "NW"  #  to add: (SD  ID  LD  ILD)

IMPR_SET = ["1","2"] # ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

STD_DERIV = false
# RESC      = false

OVERWRITE = false  # Allows to erase data and overwrite it with new data, use carefully !!

path_bdio = path_bdio_dict["local"]

@info("STARTING HVP COMPUTATION [diag. NLOc; wind. NW]")
STD_DERIV ? @info("SRANDART DERIVATIVE is being employed in the IMPROVEMENT") : nothing

@time begin
    for ens in [ensInfo[1]]

        @info("Computing HVP for ensemble $(ens.id)")
        ens.id ∈ ensNOcharm ? @info("  > NO CHARM DATA FOR $(ens.id)") : nothing
        ens.id ∈ ensNOdisc  ? @info("  > NO DISC. DATA FOR $(ens.id)") : nothing

        println("   - Reading t0...")

        t0 = BDIOread_t0(path_bdio,ens)

        println("   - Reading TMR...")

        TMR = BDIOread_TMR(path_bdio,ens,"NLOc",resc=false,beta=false,BLIND=false)
    
        (ens.id ∉ ensNOcharm) ? TMRbeta = BDIOread_TMR(path_bdio,ens,"NLOc",resc=false,beta=true,BLIND=false) : nothing

        aens = sqrtt0_ph / sqrt(t0)
        sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1)
        t   = collect(1:sym_points)
        # tBM = collect(1:length(TMR[1,:]))
        # tfm = aens.*(t.-1)

        # define windowed kernel
        # if wind == "NW" 
        #     TMRw = TMR
        #     (ens.id ∉ ensNOcharm) ? TMRwbeta = TMRbeta : nothing
        # else
        #     TMRw = TMR .* Window2D(wind)(tfm,tfm)
        #     (ens.id ∉ ensNOcharm) ? TMRwbeta = TMRbeta .* Window2D(wind)(tfm,tfm) : nothing
        # end
        TMRw = TMR[t,t]
        (ens.id ∉ ensNOcharm) ? TMRwbeta = TMRbeta[t,t] : nothing

        for impr_set in IMPR_SET

            println("   - Starting set "*impr_set)
            println("      - Reading corr...")

            corr = BDIOread_corr(path_bdio,ens,impr_set,STD=STD_DERIV)

            if ens.kappa_l == ens.kappa_s || ens.id in ensNOdisc
                println("      - SU(3) flavour sym point or no disc. data available")
            else
                corr["g88_ll"] = corr["g88conn_ll"] .+ corr["g88disc_ll"] .+ (2).*corr["g08conn_ll"] .+ corr["g08disc_ll"] .+ corr["g80disc_ll"]
                corr["g88_lc"] = corr["g88conn_lc"] .+ corr["g88disc_lc"] .+ corr["g08conn_lc"] .+ corr["g08disc_lc"]
            end

            println("      - Computing HVPs...")

            tcut0 = 10
            tstep = 1
            tcut_fix = 1.2  # fm
 
            uwerr(E0_ens[ens.id]["E0"]); E0 = E0_ens[ens.id]["E0"].mean + E0_ens[ens.id]["E0"].err

            HVP     =  Dict{String, uwreal}()
            HVPsyst =  Dict{String, Float64}()
            PlatRec =  Dict{String, Vector{Float64}}()

            for discr in ["ll","lc"] # ["ll-ll","ll-lc","lc-lc"] or ["ll","lc","mixed"]
                println("         - discr. = $discr")

                ub3333 = Vector{uwreal}(); lb3333 = Vector{uwreal}()
                ub8888 = Vector{uwreal}(); lb8888 = Vector{uwreal}()
                ub3388 = Vector{uwreal}(); lb3388 = Vector{uwreal}()
                
                println("            - Light contributions...")
                println("               - Bounding correlators...")

                Eeff33 = uwreal(0.0); Eeff88 = uwreal(0.0)
                if ens.kappa_l == ens.kappa_s || ens.id ∈ ensNOdisc
                    obs33 = corr["g33_$discr"]

                    int3333 = (obs33[t].*hcat(obs33[t]...)) .* TMRw[t,t]

                    for tcut in ProgressBar(tcut0:tstep:t[end-1])

                        if aens*tcut < tcut_fix  # we fix the eff energy at some point
                            Eeff33=Eeff(tcut, obs33)
                        end
                        UB33 = corr_bound(t, tcut, obs33, E0)
                        LB33 = corr_bound(t, tcut, obs33, Eeff33)

    
                        UBInt3333_1 = (UB33 .* hcat(UB33...)) .* TMRw[tcut+1:end,tcut+1:end]
                        UBInt3333_2 = (obs33[1:tcut] .* hcat(UB33...)) .* TMRw[1:tcut,tcut+1:end]
    
                        LBInt3333_1 = (LB33 .* hcat(LB33...)) .* TMRw[tcut+1:end,tcut+1:end]
                        LBInt3333_2 = (obs33[1:tcut] .* hcat(LB33...)) .* TMRw[1:tcut,tcut+1:end]
    
                        ub_amuNLO3333 = (alpha/pi)^3 * (sum(int3333[1:tcut,1:tcut]) + sum(UBInt3333_1) + 2*sum(UBInt3333_2)) * 1e10
            
                        lb_amuNLO3333 = (alpha/pi)^3 * (sum(int3333[1:tcut,1:tcut]) + sum(LBInt3333_1) + 2*sum(LBInt3333_2)) * 1e10
            
                        push!(ub3333, ub_amuNLO3333); push!(lb3333, lb_amuNLO3333)
                        push!(ub8888, ub_amuNLO3333); push!(lb8888, lb_amuNLO3333)
                        push!(ub3388, ub_amuNLO3333); push!(lb3388, lb_amuNLO3333)
                    end
                else
                    obs33 = corr["g33_$discr"]
                    obs88 = corr["g88_$discr"]

                    int3333 = (obs33[t].*hcat(obs33[t]...)) .* TMRw[t,t]
                    int8888 = (obs88[t].*hcat(obs88[t]...)) .* TMRw[t,t]
                    int3388 = (obs33[t].*hcat(obs88[t]...)) .* TMRw[t,t]

                    for tcut in ProgressBar(tcut0:tstep:t[end-1])

                        if aens*tcut < tcut_fix  # we fix the eff energy at some point
                            Eeff33=Eeff(tcut, obs33)
                            Eeff88=Eeff(tcut, obs88)
                        end

                        UB33 = corr_bound(t, tcut, obs33, E0)
                        LB33 = corr_bound(t, tcut, obs33, Eeff33)
                        UB88 = corr_bound(t, tcut, obs88, E0)
                        LB88 = corr_bound(t, tcut, obs88, Eeff88)

                        UBInt3333_1 = (UB33 .* hcat(UB33...)) .* TMRw[tcut+1:end,tcut+1:end]
                        UBInt8888_1 = (UB88 .* hcat(UB88...)) .* TMRw[tcut+1:end,tcut+1:end]
                        UBInt3388_1 = (UB33 .* hcat(UB88...)) .* TMRw[tcut+1:end,tcut+1:end]
                        UBInt3333_2 = (obs33[1:tcut] .* hcat(UB33...)) .* TMRw[1:tcut,tcut+1:end]
                        UBInt8888_2 = (obs88[1:tcut] .* hcat(UB88...)) .* TMRw[1:tcut,tcut+1:end]
                        UBInt3388_2 = (obs33[1:tcut] .* hcat(UB88...)) .* TMRw[1:tcut,tcut+1:end]
                        UBInt8833_2 = (obs88[1:tcut] .* hcat(UB33...)) .* TMRw[1:tcut,tcut+1:end]

                        LBInt3333_1 = (LB33 .* hcat(LB33...)) .* TMRw[tcut+1:end,tcut+1:end]
                        LBInt8888_1 = (LB88 .* hcat(LB88...)) .* TMRw[tcut+1:end,tcut+1:end]
                        LBInt3388_1 = (LB33 .* hcat(LB88...)) .* TMRw[tcut+1:end,tcut+1:end]
                        LBInt3333_2 = (obs33[1:tcut] .* hcat(LB33...)) .* TMRw[1:tcut,tcut+1:end]
                        LBInt8888_2 = (obs88[1:tcut] .* hcat(LB88...)) .* TMRw[1:tcut,tcut+1:end]
                        LBInt3388_2 = (obs33[1:tcut] .* hcat(LB88...)) .* TMRw[1:tcut,tcut+1:end]
                        LBInt8833_2 = (obs88[1:tcut] .* hcat(LB33...)) .* TMRw[1:tcut,tcut+1:end]

                        INT3333 = sum(int3333[1:tcut,1:tcut])
                        INT8888 = sum(int8888[1:tcut,1:tcut])
                        INT3388 = sum(int3388[1:tcut,1:tcut])

                        ub_amuNLO3333 = (alpha/pi)^3 * (INT3333 + sum(UBInt3333_1) + 2*sum(UBInt3333_2)) * 1e10
                        ub_amuNLO8888 = (alpha/pi)^3 * (INT8888 + sum(UBInt8888_1) + 2*sum(UBInt8888_2)) * 1e10
                        ub_amuNLO3388 = (alpha/pi)^3 * (INT3388 + sum(UBInt3388_1) + sum(UBInt3388_2) + sum(UBInt8833_2)) * 1e10
            
                        lb_amuNLO3333 = (alpha/pi)^3 * (INT3333 + sum(LBInt3333_1) + 2*sum(LBInt3333_2)) * 1e10
                        lb_amuNLO8888 = (alpha/pi)^3 * (INT8888 + sum(LBInt8888_1) + 2*sum(LBInt8888_2)) * 1e10
                        lb_amuNLO3388 = (alpha/pi)^3 * (INT3388 + sum(LBInt3388_1) + sum(LBInt3388_2) + sum(LBInt8833_2)) * 1e10
            
                        push!(ub3333, ub_amuNLO3333); push!(lb3333, lb_amuNLO3333)
                        push!(ub8888, ub_amuNLO8888); push!(lb8888, lb_amuNLO8888)
                        push!(ub3388, ub_amuNLO3388); push!(lb3388, lb_amuNLO3388)
                    end
                end

                ub = [ub3333,ub8888,ub3388]
                lb = [lb3333,lb8888,lb3388]

                println("               - Applying BM...")

                for (i,comp) in enumerate(["3333","8888","3388"])
                    HVP["g$(comp)_$(discr)$(discr)"], HVPsyst["g$(comp)_$(discr)$(discr)"], PlatRec["g$(comp)_$(discr)$(discr)"] = bounding_method(ub[i],lb[i],aens,PLAT=true,AVER=false,tcut0=tcut0)
                end

                if ens.id ∉ ensNOcharm
                    println("            - Charmed contributions.")

                    obsCC  = corr["gCCconn_$(discr)_sim"]
                    obsCCp = corr["gCCconn_$(discr)_sim+"]

                    HVP["gCCCC_$(discr)$(discr)_sim"]  = (alpha/pi)^3 * sum((obsCC[t] .* hcat(obsCC[t]...)) .* TMRwbeta) * 1e10
                    HVP["gCCCC_$(discr)$(discr)_sim+"] = (alpha/pi)^3 * sum((obsCCp[t] .* hcat(obsCCp[t]...)) .* TMRw) * 1e10
                    HVP["g33CC_$(discr)$(discr)_sim"]  = (alpha/pi)^3 * sum((obs33[t] .* hcat(obsCC[t]...)) .* TMRw) * 1e10
                    HVP["g33CC_$(discr)$(discr)_sim+"] = (alpha/pi)^3 * sum((obs33[t] .* hcat(obsCCp[t]...)) .* TMRw) * 1e10
                    if ens.kappa_l != ens.kappa_s
                        HVP["g88CC_$(discr)$(discr)_sim"]  = (alpha/pi)^3 * sum((obs88[t] .* hcat(obsCC[t]...)) .* TMRw) * 1e10
                        HVP["g88CC_$(discr)$(discr)_sim+"] = (alpha/pi)^3 * sum((obs88[t] .* hcat(obsCCp[t]...)) .* TMRw) * 1e10
                    else
                        HVP["g88CC_$(discr)$(discr)_sim"]  = HVP["g33CC_$(discr)$(discr)_sim"] 
                        HVP["g88CC_$(discr)$(discr)_sim+"] = HVP["g33CC_$(discr)$(discr)_sim+"]
                    end

                    println("               - Reading kappaC target...")

                    # kappaC_tar = BDIOread_KappaC_tar(path_bdio,ens)
                    mDs_ph_prime, mDs_beta, Ds_dict = BDIOread_mDs_kappaC(path_bdio,ens)
                    kappaC_tar = Ds_dict["kappaC"]

                    println("               - Interpolating charm contributions...")

                    # kappaC_tar = uwreal([kcd_in[ens.id]["kappaC"],kcd_in[ens.id]["kappaC_err"]], "kappaC target")   # extract KappaC from tables
                    kappaC = [kcd_in[ens.id]["kappaC_sim"],kcd_in[ens.id]["kappaC_sim_plus"]]
                    
                    DkappaC = abs(kappaC[1]-kappaC[2])
                    DkappaC_tar = minimum([abs(value(kappaC_tar)-kappaC[1]),abs(value(kappaC_tar)-kappaC[2])])

                    for key in (ens.kappa_l != ens.kappa_s ? ["gCCCC_$(discr)$(discr)","g33CC_$(discr)$(discr)","g88CC_$(discr)$(discr)"] : ["gCCCC_$(discr)$(discr)","g33CC_$(discr)$(discr)"])
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
                    if ens.kappa_l == ens.kappa_s
                        HVP["g88CC_$(discr)$(discr)"]  = HVP["g33CC_$(discr)$(discr)"] 
                    end
                end
            end

            println("      - Writing BDIO & JDL2...")

            DERstr = STD_DERIV ? "_std" : "" 
            # RESstr = RESC ? "_resc" : ""
            pBDIO = create_path(path_bdio,["HVP&FVC","NW",ens.id,"$(ens.id)_HVPNLOc_set$(impr_set)$(DERstr)"],OVERWRITE=OVERWRITE)
            pjdl2 = create_path(path_bdio,["HVP&FVC","NW",ens.id,"$(ens.id)_HVPNLOc_info_set$(impr_set)$(DERstr).jld2"],OVERWRITE=OVERWRITE)

            io = IOBuffer()
            write(io, "$(ens.id) HVP")
            fb = ALPHAdobs_create(pBDIO, io)

            extra = Dict{String, Any}("ens" => ens.id, "impr_set" => impr_set, "diag" => "NLOc", "wind" => "NW")
            ALPHAdobs_write(fb, HVP, extra=extra)

            ALPHAdobs_close(fb)

            HVPinfo = Dict{String,Dict}(
                "HVPsyst" => HVPsyst,
                "plateau/Reconstr" => PlatRec,
            )
            save(pjdl2,"HVPinfo",HVPinfo)

        end # end impr_set loop
    end # end ens loop
end # end timer


## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> 2D FVC computation [NLOc] <==========================##

# diag = "NLOc"
# wind = "NW"  #  to add: (SD  ID  LD  ILD)

IMPR_SET = ["1","2"] # ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

STD_DERIV = false
# RESC      = false
# VREF      = false

OVERWRITE = false  # Allows to erase data and overwrite it with new data, use carefully !!

path_bdio = path_bdio_dict["local"]

@info("STARTING FVC COMPUTATION [diag. NLOc; wind. NW]")
STD_DERIV ? @info("SRANDART DERIVATIVE is being employed in the IMPROVEMENT") : nothing

@time begin
    for ens in [ensInfo[5]]

        @info("Computing for ensemble $(ens.id)")

        println("   - Reading TMR...")

        TMR = BDIOread_TMR(path_bdio,ens,"NLOc",resc=false,beta=false,BLIND=false)

        # aens = sqrtt0_ph / sqrt(t0)
        sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1)
        t   = collect(1:sym_points)
        # tBM = collect(1:length(TMR[1,:]))
        # tfm = aens.*(t.-1)

        # define windowed kernel
        # if wind == "NW" 
        #     TMRw = TMR
        # else
        #     TMRw = TMR .* Window2D(wind)(tfm,tfm)
        # end
        TMRw = TMR[t,t]

        println("   - Reading fvc (HP)...")
        fvc_hp = BDIOread_FVCcorr(path_bdio,ens,Vref=false)

        if ens.id != "A653"
            println("   - Reading fvc (MLL)...")
            t_mll, fvc_mll = TXTread_FVCcorr_GS(pFVC_MLL,ens,Vref=false)
        end

        fvcPi_hp = vcat(TMRw[1],-fvc_hp["FVCPi6"])[t]
        fvcPi = ens.id != "A653" ? vcat(fvcPi_hp[1:Int64(t_mll[1])],-fvc_mll...)[t] : fvcPi_hp
        fvcK  = vcat(TMRw[1],-fvc_hp["FVCK6" ])[t]

        println("   - Computing FVC...")

        for impr_set in IMPR_SET

            println("      - Starting set "*impr_set)
            println("         - Reading corr...")

            println("      - Reading corr...")
            corr = BDIOread_corr(path_bdio,ens,impr_set,STD=STD_DERIV)

            if ens.kappa_l == ens.kappa_s
                println("         - SU(3) flavour sym point!")

                corr["g88_ll"] = corr["g88conn_ll"]
                corr["g88_lc"] = corr["g88conn_lc"]

                fvc33_hp = (3/2) * fvcPi_hp
                fvc33 = (3/2) * fvcPi
            else
                corr["g88_ll"] = corr["g88conn_ll"] .+ corr["g88disc_ll"] .+ (2).*corr["g08conn_ll"] .+ corr["g08disc_ll"] .+ corr["g80disc_ll"]
                corr["g88_lc"] = corr["g88conn_lc"] .+ corr["g88disc_lc"] .+ corr["g08conn_lc"] .+ corr["g08disc_lc"]

                fvc33_hp = fvcPi_hp + fvcK/2
                fvc33 = fvcPi + fvcK/2
                fvc88 = 2/3 * fvcK
            end

            FVC     = Dict{String, uwreal}()
            FVCsyst = Dict{String, uwreal}()

            println("         - Computing FVC...")

            for discr in ["ll","lc"]
                println("            - discr. $discr")

                println("            -  Light contributions")

                int_g3333 = TMRw .* (corr["g33_$discr"][t].*hcat(fvc33...) .+ fvc33.*hcat(corr["g33_$discr"][t]...) .+ fvc33.*hcat(fvc33...))

                FVCg3333 = (alpha/pi)^3 * sum(int_g3333) * 1e10

                if ens.kappa_l != ens.kappa_s
                    int_g8888 = TMRw .* (corr["g88_$discr"][t].*hcat(fvc88...) .+ fvc88.*hcat(corr["g88_$discr"][t]...) .+ fvc88.*hcat(fvc88...))
                    int_g3388 = TMRw .* (corr["g33_$discr"][t].*hcat(fvc88...) .+ fvc88.*hcat(corr["g33_$discr"][t]...) .+ fvc33.*hcat(fvc88...))

                    FVCg8888 = (alpha/pi)^3 * sum(int_g8888) * 1e10
                    FVCg3388 = (alpha/pi)^3 * sum(int_g3388) * 1e10
                else
                    FVCg8888 = FVCg3388 = FVCg3333
                end

                FVC["FVCg3333_$(discr)$(discr)"] = FVCg3333
                FVC["FVCg8888_$(discr)$(discr)"] = FVCg8888
                FVC["FVCg3388_$(discr)$(discr)"] = FVCg3388

                # only HP; comment if too slow

                int_g3333_hp = TMRw .* (corr["g33_$discr"][t].*hcat(fvc33_hp...) .+ fvc33_hp.*hcat(corr["g33_$discr"][t]...) .+ fvc33_hp.*hcat(fvc33_hp...))

                FVCg3333_hp = (alpha/pi)^3 * sum(int_g3333_hp) * 1e10

                if ens.kappa_l != ens.kappa_s
                    int_g3388_hp = TMRw .* (corr["g33_$discr"][t].*hcat(fvc88...) .+ fvc88.*hcat(corr["g33_$discr"][t]...) .+ fvc33_hp.*hcat(fvc88...))
                    
                    FVCg3388_hp = (alpha/pi)^3 * sum(int_g3388_hp) * 1e10
                else
                    int_g3333_hp = TMRw .* (corr["g33_$discr"][t].*hcat(fvc33_hp...) .+ fvc33_hp.*hcat(corr["g33_$discr"][t]...) .+ fvc33_hp.*hcat(fvc33_hp...))

                    FVCg3388_hp = FVCg3333_hp
                end

                FVC["FVCg3333_hp_$(discr)$(discr)"] = FVCg3333_hp
                FVC["FVCg3388_hp_$(discr)$(discr)"] = FVCg3388_hp

                if ens.id ∉ ensNOcharm
                    println("            -  Charmed contributions")

                    int_g33CC_sim  = TMRw .* (fvc33.*hcat(corr["gCCconn_$(discr)_sim"][t]...))
                    int_g33CC_simp = TMRw .* (fvc33.*hcat(corr["gCCconn_$(discr)_sim+"][t]...))

                    FVCg33CC_sim  = (alpha/pi)^3 * sum(int_g33CC_sim ) * 1e10
                    FVCg33CC_simp = (alpha/pi)^3 * sum(int_g33CC_simp) * 1e10

                    if ens.kappa_l != ens.kappa_s
                        int_g88CC_sim  = TMRw .* (fvc88.*hcat(corr["gCCconn_$(discr)_sim"][t]...))
                        int_g88CC_simp = TMRw .* (fvc88.*hcat(corr["gCCconn_$(discr)_sim+"][t]...))

                        FVCg88CC_sim  = (alpha/pi)^3 * sum(int_g88CC_sim ) * 1e10
                        FVCg88CC_simp = (alpha/pi)^3 * sum(int_g88CC_simp) * 1e10
                    else
                        FVCg88CC_sim  = FVCg33CC_sim
                        FVCg88CC_simp = FVCg33CC_simp
                    end

                    FVC["FVCg33CC_$(discr)$(discr)_sim"]  = FVCg33CC_sim
                    FVC["FVCg33CC_$(discr)$(discr)_sim+"] = FVCg33CC_simp
                    FVC["FVCg88CC_$(discr)$(discr)_sim"]  = FVCg88CC_sim
                    FVC["FVCg88CC_$(discr)$(discr)_sim+"] = FVCg88CC_simp

                    # only HP; comment if too slow

                    int_g33CC_hp_sim  = TMRw .* (fvc33_hp.*hcat(corr["gCCconn_$(discr)_sim"][t]...))
                    int_g33CC_hp_simp = TMRw .* (fvc33_hp.*hcat(corr["gCCconn_$(discr)_sim+"][t]...))

                    FVCg33CC_hp_sim  = (alpha/pi)^3 * sum(int_g33CC_hp_sim ) * 1e10
                    FVCg33CC_hp_simp = (alpha/pi)^3 * sum(int_g33CC_hp_simp) * 1e10

                    FVC["FVCg33CC_hp_$(discr)$(discr)_sim"]  = FVCg33CC_hp_sim
                    FVC["FVCg33CC_hp_$(discr)$(discr)_sim+"] = FVCg33CC_hp_simp

                    println("               - Reading kappaC target...")

                    # kappaC_tar = BDIOread_KappaC_tar(path_bdio,ens)
                    mDs_ph_prime, mDs_beta, Ds_dict = BDIOread_mDs_kappaC(path_bdio,ens)
                    kappaC_tar = Ds_dict["kappaC"]

                    println("               - Interpolating charm contributions...")

                    # kappaC_tar = uwreal([kcd_in[ens.id]["kappaC"],kcd_in[ens.id]["kappaC_err"]], "kappaC target")   # extract KappaC from tables
                    kappaC = [kcd_in[ens.id]["kappaC_sim"],kcd_in[ens.id]["kappaC_sim_plus"]]
                    
                    DkappaC = abs(kappaC[1]-kappaC[2])
                    DkappaC_tar = minimum([abs(value(kappaC_tar)-kappaC[1]),abs(value(kappaC_tar)-kappaC[2])])

                    for key in (ens.kappa_l != ens.kappa_s ? ["FVCg33CC_$(discr)$(discr)","FVCg88CC_$(discr)$(discr)","FVCg33CC_hp_$(discr)$(discr)"] : ["FVCg33CC_$(discr)$(discr)","FVCg33CC_hp_$(discr)$(discr)"])
                        obs  = [FVC["$(key)_sim"],FVC["$(key)_sim+"]]
                        
                        par, _ = lin_fit(kappaC,obs,wpm=wpm,lineprint=false)
                        obs_tar = y_lin_fit(par,kappaC_tar)

                        if kappaC[1] <= value(kappaC_tar) <= kappaC[2]
                            FVC[key] = obs_tar
                            FVCsyst[key] = 0.0
                        else
                            FVC[key] = obs_tar # + 1.5 * (abs(DkappaC_tar) / DkappaC)^2 * uwreal([0.0,err(obs_tar)],"kappaC syst")
                            uwerr(obs_tar); FVCsyst[key] = 1.5 * (abs(DkappaC_tar) / DkappaC)^2 * err(obs_tar)
                        end
                    end
                    if ens.kappa_l == ens.kappa_s
                        FVC["FVCg88CC_$(discr)$(discr)"]  = FVC["FVCg33CC_$(discr)$(discr)"] 
                    end
                end
            end

            println("   - Writing BDIO...")

            DERstr = STD_DERIV ? "_std" : "" 
            # RESstr = RESC ? "_resc" : ""
            pBDIO = create_path(path_bdio,["HVP&FVC","NW",ens.id,"$(ens.id)_FVCNLOc_set$(impr_set)$(DERstr)"],OVERWRITE=OVERWRITE)
            pjdl2 = create_path(path_bdio,["HVP&FVC","NW",ens.id,"$(ens.id)_FVCNLOc_info_set$(impr_set)$(DERstr).jld2"],OVERWRITE=OVERWRITE)

            io = IOBuffer()
            write(io, "$(ens.id) FVC")
            fb = ALPHAdobs_create(pBDIO, io)

            extra = Dict{String, Any}("ens" => ens.id, "diag" => "NLOc", "impr_set" => impr_set, "wind" => "NW", "Vref" => false)
            ALPHAdobs_write(fb, FVC, extra=extra)

            ALPHAdobs_close(fb)

            FVCinfo = Dict{String,Dict}(
                "FVCsyst" => FVCsyst,
            )
            save(pjdl2,"FVCinfo",FVCinfo)
        end
    end
end # end timer


## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##


##==========================> READING TEST <==========================##

# diag     = "NLOc"
# wind     = "NW"  #  "NW"  "SD"  "ID"  "LD"  "ILD"
ensid    = "C101"
impr_set = "1"

STD   = false
# VREF  = false
# RESC  = false

path_bdio = path_bdio_dict["local"]

HVP, info = BDIOread_HVPens(path_bdio,"NLOc","NW",ensid,impr_set,info=true,resc=false,STD=STD,BLIND=false)

FVC = BDIOread_FVCens(path_bdio,"NLOc","NW",ensid,Vref=false,resc=false,BLIND=false)


HVP
FVC["1"]