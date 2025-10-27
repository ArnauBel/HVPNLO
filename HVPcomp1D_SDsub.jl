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
using TimerOutputs

# include uwreal constants

# include("HVPtool/uwConst.jl")

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_tl   = joinpath(julia_script_directory, "..", "HVPData", "tree_level_improv")
path_coef = joinpath(julia_script_directory, "..", "KernelCoeff")

# Blind analysis (Simon K.) safe ensembles: 
# SU(3) sym.  H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

# All considered ensembles are:
ensList = [
    "A653","A654","B450","C101","C102",
    "D150","D200","D201","D251","D450",
    "D451","D452","E250","E300","F300",
    "H101","H102","H200","J303","J304",
    "J306","J307","J500","J501","N101",
    "N200","N202","N203","N300","N302",
    "N451","N452","S400"]

ensInfo = EnsInfo.(ensList)

# We do not have charm or disconnected data for some of the ensembles
ensNOcharm  = ["C102","D150","D201","D251","D451","F300","H200","J304","J306","J307","J501","N451","N452"]
ensNOdisc   = ["F300","J306"]

@info("Ready")


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> 1D HVP computation [LO, NLOa, NLOb] <==========================##

diag = ""  # LO  NLOa  NLOb  NLOa&b

IMPR_SET = ["1","2"] # ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

STD_DERIV = false
RESC      = false

OVERWRITE = false  # Allows to erase data and overwrite it with new data, use carefully !!

path_bdio = path_bdio_dict["local"]


@info("STARTING HVP COMPUTATION [diag. $diag; wind. subtracted SD]")
STD_DERIV ? @info("STANDARD DERIVATIVE is being employed in the IMPROVEMENT") : nothing


# LO => (alpha/pi)^2 || NLO => (alpha/pi)^3
exp_diag = diag == "LO" ? 2 : 3

corr33tl_ll, corr33tl_lc = read_tree_level_v33(path_tl, cons=true)
corr33tl_v3s03_ll, corr33tl_v3s03_lc = read_tree_level_v3sig03(path_tl, cons=true, massless=true)

corr33tl_ll = uwreal.(corr33tl_ll); corr33tl_lc = uwreal.(corr33tl_lc)

@time begin
    for ens in ensInfo

        @info("Computing for ensemble $(ens.id)")
        ens.id ∈ ensNOcharm ? @info("  > NO CHARM DATA FOR $(ens.id)") : nothing
        ens.id ∈ ensNOdisc  ? @info("  > NO DISC. DATA FOR $(ens.id)") : nothing

        T = HVPobs.Data.get_T(ens.id)
        sym_points = Int64(T/2+1); t = collect(1:sym_points)
        if !RESC
            println("   - Reading t0...")
            t0 = BDIOread_t0(path_bdio,ens)
            aens = sqrtt0_ph / sqrt(t0); tfm = aens.*(t.-1)
            factor = hbarc * sqrt(t0) / sqrtt0_ph  
            if (ens.id ∉ ensNOcharm) || (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)
                a_ß = sqrtt0_ph / sqrt(t0sym(ens.beta)); tfm_ß = a_ß.*(t.-1)
                factor_ß = hbarc * sqrt(t0sym(ens.beta)) / sqrtt0_ph
            end
        else
            println("   - Reading fPi...")
            fPi = BDIOread_fPS(path_bdio,ens)["fPi"]
            aens = hbarc * fPi / fPi_ph; tfm = aens.*(t.-1)
            factor = fPi_ph / fPi  
            if (ens.id ∉ ensNOcharm) || (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)
                a_ß = hbarc * fPiph(ens.beta) / fPi_ph; tfm_ß = a_ß.*(t.-1)
                factor_ß = fPi_ph / fPiph(ens.beta)
            end
        end

        println("   - Reading TMR...")
        TMR  = BDIOread_TMR(path_bdio,ens,diag,resc=RESC,beta=false); TMR = TMR[t]
        if (ens.id ∉ ensNOcharm) || (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)
            TMRbeta = BDIOread_TMR(path_bdio,ens,diag,resc=RESC,beta=true); TMRbeta = TMRbeta[t]
        end

        TMRw = TMR .* Window("SD")(tfm)

        TMRb(Q::Float64) = ((16/(Q/factor)^2)^2 * π^2 * (massmu/factor)^2) .* C4[diag].((massmu/factor).*(t .- 1)) .* sin.((Q/factor/4) .* (t .- 1)).^4
        TMRsub(Q::Float64) = TMRw .- (Window("SD")(0) .* TMRb(Q))
        if (ens.id ∉ ensNOcharm) || (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)
            TMRwbeta = TMRbeta .* Window("SD")(tfm_ß)

            TMRb_SU3(Q::Float64)   = ((16/(Q/factor_ß)^2)^2 * π^2 * (massmu/factor_ß)^2) .* C4[diag].((massmu/factor_ß).*(t .- 1)) .* sin.((Q/factor_ß/4) .* (t .- 1)).^4
            TMRsub_SU3(Q::Float64) = TMRwbeta .- (Window("SD")(0) .* TMRb_SU3(Q))
        end
        

        for impr_set in IMPR_SET

            println("   - Starting set "*impr_set)

            println("      - Reading corr...")
            corr = BDIOread_corr(path_bdio,ens,impr_set,STD=STD_DERIV)

            if (ens.kappa_l != ens.kappa_s) && (ens.id ∉ ensNOdisc)
                corr["g88_ll"] = corr["g88conn_ll"] .+ corr["g88disc_ll"] .+ (2).*corr["g08conn_ll"] .+ corr["g08disc_ll"] .+ corr["g80disc_ll"]
                corr["g88_lc"] = corr["g88conn_lc"] .+ corr["g88disc_lc"] .+ corr["g08conn_lc"] .+ corr["g08disc_lc"]
            end

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

            println("      - Computing HVPs...\n        For Q = $Qlist")

            HVPQ = Dict{String, Array{uwreal}}()
            HVP = Dict{String, Array{uwreal}}()
            HVPsyst = Dict{String, Array{Float64}}()

            # HVPQ["g33_tl0"] = tl_cont

            for  key in ["g33_ll","g33_lc","g33tl_ll","g33tl_lc"]
                println("         - $key")

                obs   = corr[key][t]
                HVPQ[key] = []
                for Q in Qlist

                    int = obs .* TMRsub(Q)
                    amu = (alpha/pi)^exp_diag * sum( obs .* TMRsub(Q)) * 1e10

                    push!(HVPQ[key],amu)
                end
            end


            if ens.id ∉ ensNOcharm
                for key in ["gCCconn_ll_sim","gCCconn_lc_sim","gCCconn_ll_sim+","gCCconn_lc_sim+"]
                    println("         - $key  &  $(key[1:8])SU3$(key[8:end])  &  ∆lc_b$(key[8:end])")

                    obs  = corr[key][t]
                    ∆obs = obs .- 2 .* corr["g33_"*key[9:10]][t]

                    HVPQ[key] = []
                    HVPQ[key[1:8]*"SU3"*key[8:end]] = []
                    HVPQ["∆lc_b"*key[8:end]] = []
                    # HVPQ["∆lc_b_"*key[5:6]*"_beta_"*key[13:end]] = []
                    for Q in Qlist
                        int    = obs .* TMRsub(Q)
                        intSU3 = obs .* TMRsub_SU3(Q)
                        intb   = ∆obs .* TMRb(Q)
                        # intbSU3 = ∆obs .* TMRb_SU3(Q)

                        amu    = (alpha/pi)^exp_diag * sum(int) * 1e10
                        amuSU3 = (alpha/pi)^exp_diag * sum(intSU3) * 1e10
                        amub   = (alpha/pi)^exp_diag * sum(intb) * 1e10
                        # amubSU3 = (alpha/pi)^exp_diag * sum(intbSU3) * 1e10

                        push!(HVPQ[key], amu)
                        push!(HVPQ[key[1:8]*"SU3"*key[8:end]], amuSU3)
                        push!(HVPQ["∆lc_b"*key[8:end]], amub)
                        # push!(HVPQ["∆lc_b_"*key[5:6]*"_beta_"*key[13:end]], amubSU3)
                    end
                end

                println("      - Reading kappaC target...")

                mDs_ph_prime, mDs_beta, Ds_dict = BDIOread_mDs_kappaC(path_bdio,ens)
                kappaC_tar = Ds_dict["kappaC"]
                println("         - Interpolating charm contribution...")

                kappaC = [kcd_in[ens.id]["kappaC_sim"],kcd_in[ens.id]["kappaC_sim_plus"]]
                
                DkappaC = abs(kappaC[1]-kappaC[2])
                DkappaC_tar = minimum([abs(value(kappaC_tar)-kappaC[1]),abs(value(kappaC_tar)-kappaC[2])])

                for key in ["gCCconn_ll","gCCconn_lc","gCCconn_SU3_ll","gCCconn_SU3_lc","∆lc_b_ll","∆lc_b_lc"]
                    HVPQ[key] = []
                    HVPsyst[key] = []
                    for i=1:length(Qlist)
                        obs  = [HVPQ["$(key)_sim"][i],HVPQ["$(key)_sim+"][i]]
                    
                        par, _ = lin_fit(kappaC,obs,lineprint=false)
                        obs_tar = y_lin_fit(par,kappaC_tar)

                        if kappaC[1] <= value(kappaC_tar) <= kappaC[2]
                            push!(HVPQ[key], obs_tar)
                            push!(HVPsyst[key], 0.0)
                        else
                            push!(HVPQ[key], obs_tar)  # + 1.5 * (abs(DkappaC_tar) / DkappaC)^2 * uwreal([0.0,err(obs_tar)],"kappaC syst")
                            uwerr(obs_tar); push!(HVPsyst[key], 1.5 * (abs(DkappaC_tar) / DkappaC)^2 * err(obs_tar))
                        end
                    end
                end
            end

            if (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)
                println("      - Computing ∆ls(amu)...")
                for discr in ["ll","lc"]
                    key = "∆ls_amu_$discr"
                    println("         - $key")
                    obs = corr["g88_$discr"][t] .- corr["g33_$discr"][t]

                    int = obs .* TMRw
                    amu = (alpha/pi)^exp_diag * sum(int) * 1e10

                    HVP[key] = [amu]
                end
                println("      - Computing gCdisc...")
                for (k,key) in enumerate(["gCCdisc_cc","gC8disc_cc"])
                    println("         - $key")

                    obs = corr[key][t]

                    int = k==1 ? (obs .* TMRwbeta) : (obs .* TMRw)
                    amu = (alpha/pi)^exp_diag * sum(int) * 1e10

                    HVP[key] = [amu]
                end
            else
                println("         - 'G33 = G88 & GCCdisc = 0' for $(ens.id)")
            end

            println("      - Writing BDIO & JDL2...")

            DERstr = STD_DERIV ? "_std" : "" 
            RESstr = RESC ? "_resc" : ""
            pBDIO  = create_path(path_bdio,["HVP&FVC","SDsub",ens.id,"$(ens.id)_HVP$(diag)_set$(impr_set)$(RESstr)$(DERstr)"],OVERWRITE=OVERWRITE)
            pjdl2  = create_path(path_bdio,["HVP&FVC","SDsub",ens.id,"$(ens.id)_HVP$(diag)_info_set$(impr_set)$(RESstr)$(DERstr).jld2"],OVERWRITE=OVERWRITE)

            io = IOBuffer()
            write(io, "$(ens.id) HVP")
            fb = ALPHAdobs_create(pBDIO, io)

            extra = Dict{String, Any}("ens" => ens.id, "impr_set" => impr_set, "diag" => diag, "wind" => "SDsub", "Qlist" => Qlist)
            ALPHAdobs_write(fb, HVPQ, extra=extra)
            (ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc) ? ALPHAdobs_write(fb, HVP, extra=extra) : nothing

            ALPHAdobs_close(fb)

            HVPinfo = Dict{String,Dict}("HVPsyst" => HVPsyst) 
            save(pjdl2,"HVPinfo",HVPinfo)
        end # end impr_set loop
    end # end ens loop
end # end timer


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> 1D FVC computation [LO, NLOa, NLOb] <==========================##

diag = ""  # LO  NLOa  NLOb  NLOa&b

OVERWRITE = false  # Allows to erase data and overwrite it with new data, use carefully !!

RESC = false
VREF = false

path_bdio = path_bdio_dict["local"]

if !VREF
    @info("STARTING FVC COMPUTATION [diag. $diag; wind. SDsub]")
else
    @info("STARTING FVC COMPUTATION (to Vref) [diag. $diag; wind. SDsub]")
end


@time begin
    for ens in ensInfo

        @info("Computing for ensemble $(ens.id)")

        sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1); t = collect(1:sym_points)

        if !RESC
            println("   - Reading t0...")
            t0 = BDIOread_t0(path_bdio,ens)
            aens = sqrtt0_ph / sqrt(t0); tfm = aens.*(t.-1)
            factor = hbarc * sqrt(t0) / sqrtt0_ph  
            # if (ens.id ∉ ensNOcharm) || (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)
            #     a_ß = sqrtt0_ph / sqrt(t0sym(ens.beta)); tfm_ß = a_ß.*(t.-1)
            #     factor_ß = hbarc * sqrt(t0sym(ens.beta)) / sqrtt0_ph
            # end
        else
            println("   - Reading fPi...")
            fPi = BDIOread_fPS(path_bdio,ens)["fPi"]
            aens = hbarc * fPi / fPi_ph; tfm = aens.*(t.-1)
            factor = fPi_ph / fPi  
            # if (ens.id ∉ ensNOcharm) || (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)
            #     a_ß = hbarc * fPiph(ens.beta) / fPi_ph; tfm_ß = a_ß.*(t.-1)
            #     factor_ß = fPi_ph / fPiph(ens.beta)
            # end
        end

        println("   - Reading TMR...")
        TMR = BDIOread_TMR(path_bdio,ens,diag,beta=false,resc=RESC); TMR = TMR[t]
        # if ens.id ∉ ensNOcharm
        #     TMR_beta = BDIOread_TMR(path_bdio,ens,diag,beta=true)
        # end

        TMRw = TMR .* Window("SD")(tfm)

        TMRb(Q::Float64) = ((16/(Q/factor)^2)^2 * π^2 * C4[diag].((massmu/factor).*(t .- 1)) * (massmu/factor)^2) .* sin.((Q/factor/4) .* (t .- 1)).^4
        TMRsub(Q::Float64) = TMRw .- (Window("SD")(0) .* TMRb(Q))

        # if ens.id ∉ ensNOcharm
        #     TMRwbeta = TMR_beta .* Window("SD")(tfm)
        #     TMRb_SU3(Q::Float64)   = ((16/(Q/factor_ß)^2)^2 * C4[diag] * (massmu/factor_ß)^2) .* sin.((Q/factor_ß/4) .* (t .- 1)).^4
        # end

        println("   - Reading FVC...")
        FVC_corr = BDIOread_FVCcorr(path_bdio,ens,Vref=VREF)
        
        # LO => (alpha/pi)^2 || NLO => (alpha/pi)^3
        exp_diag = diag == "LO" ? 2 : 3

        println("   - Computing (TMRw x FVC)...")

        tFVC = 1:Int64(length(FVC_corr["FVCPi1"])/2)

        # for the substracted computation, we only compute the full FVC (i.e. n=6)
        fvcPi = vcat(TMRw[1],-FVC_corr["FVCPi6"][tFVC])
        fvcK  = vcat(TMRw[1],-FVC_corr["FVCK6" ][tFVC])

        FVCPi = (alpha/pi)^exp_diag * sum(fvcPi .* TMRw) * 1e10
        FVCK  = (alpha/pi)^exp_diag * sum(fvcK  .* TMRw) * 1e10

        FVC∆ls_amu = ens.kappa_l == ens.kappa_s ? (uwreal([0.0,0.0],"")) : ((-1.0) * (alpha/pi)^exp_diag * sum((-(1/6.).*fvcK .+ fvcPi) .* TMRw) * 1e10)

        dataFVC =  Dict{String, Array{uwreal}}(
            "FVCPi" => [FVCPi], 
            "FVCK"  => [FVCK/2],
            "FVC∆ls_amu" => [FVC∆ls_amu], 
        )

        println("   - Computing (TMRsub(Q) x FVC)...")

        FVCPisub = Vector{uwreal}(); FVCKsub = Vector{uwreal}()
        FVCPib   = Vector{uwreal}(); FVCKb   = Vector{uwreal}()
        FVC33    = Vector{uwreal}(); FVC88   = Vector{uwreal}()
        FVC∆lc_b = Vector{uwreal}()# ; FVC∆lc_bß = Vector{uwreal}()
        for Q in Qlist

            FVCPisub_ = (alpha/pi)^exp_diag * sum(fvcPi .*  TMRsub(Q)) * 1e10
            FVCKsub_  = (alpha/pi)^exp_diag * sum(fvcK  .*  TMRsub(Q)) * 1e10

            FVCPib_ = (alpha/pi)^exp_diag * sum(fvcPi .*  TMRb(Q)) * 1e10
            FVCKb_  = (alpha/pi)^exp_diag * sum(fvcK  .*  TMRb(Q)) * 1e10

            # FVCPibß_ = (alpha/pi)^exp_diag * sum(fvcPi .*  TMRb_SU3(Q)) * 1e10
            # FVCKbß_  = (alpha/pi)^exp_diag * sum(fvcK  .*  TMRb_SU3(Q)) * 1e10

            if ens.kappa_l == ens.kappa_s
                FVC33_ = (3/2)*FVCPisub_
                FVC88_ = (3/2)*FVCPisub_
                FVC∆lc_b_  = -(3)*FVCPib_ # this minus sign was not there before. !!??
                # FVC∆lc_bß_ = -(3)*FVCKbß_ # this minus sign was not there before. !!??
            else
                FVC33_ = FVCPisub_ + FVCKsub_/2
                FVC88_ = 2/3 * FVCKsub_
                FVC∆lc_b_  = (-2)*(FVCPib_ + FVCKb_/2) # there's a factor 2
                # FVC∆lc_bß_ = (-2)*(FVCPibß_ + FVCKbß_) # there's a factor 2
            end
                        
            push!(FVCPisub, FVCPisub_); push!(FVCKsub, FVCKsub_)
            push!(FVCPib, FVCPib_); push!(FVCKb, FVCKb_)
            push!(FVC33, FVC33_); push!(FVC88, FVC88_)
            push!(FVC∆lc_b, FVC∆lc_b_)
            # push!(FVC∆lc_bß, FVC∆lc_bß_)
        end

        dataFVCQ =  Dict{String, Array{uwreal}}(
            "FVCPisub" => FVCPisub, 
            "FVCKsub"  => FVCKsub/2,
            "FVCPib"   => FVCPib, 
            "FVCKb"    => FVCKb/2,
            "FVCg33"   => FVC33, 
            "FVCg88"   => FVC88,
            "FVC∆lc_b" => FVC∆lc_b
            # "FVC∆lc_b_beta" => FVC∆lc_bß
        )

        println("   - Writing BDIO...")

        REFstr = VREF ? "_Vref" : ""
        RESstr = RESC ? "_resc" : ""
        pBDIO = create_path(path_bdio,["HVP&FVC","SDsub",ens.id,"$(ens.id)_FVC$(diag)$(RESstr)$(REFstr)"],OVERWRITE=OVERWRITE)

        io = IOBuffer()
        write(io, "$(ens.id) FVC")
        fb = ALPHAdobs_create(pBDIO, io)

        extra = Dict{String, Any}("ens" => ens.id, "diag" => diag, "wind" => "SDsub", "Qlist" => Qlist)
        ALPHAdobs_write(fb, dataFVC, extra=extra)
        ALPHAdobs_write(fb, dataFVCQ, extra=extra)

        ALPHAdobs_close(fb)
    end
end # end timer

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##


##==========================> READING TEST <==========================##

diag  = ""
ensid = ""

impr_set = ""

HVP, info = BDIOread_HVPens(path_bdio,diag,"SDsub",ensid,impr_set,info=true)

FVC = BDIOread_FVCens(path_bdio,diag,"SDsub",ensid)
