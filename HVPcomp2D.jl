# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

using ProgressBars
using Suppressor
using TimerOutputs

# Include necessary functions

include("HVPtools/Utils.jl")

include("HVPtools/Reader.jl")

include("HVPtools/Writer.jl")

# BDIO path definition (set 'STD_DERIV = true' to use the standard sym. derivative in the impr.)

julia_script_directory = @__DIR__

path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")

# Blind analysis (Simon K.) safe ensembles: 
# SU(3) sym.  H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

# All considered ensembles are:
ensList = ["A654","B450","C101","C102","D150","D200","D201","D450","D451","D452","E250","E300","F300","H101","H102","H200","J303","J304","J306","J307","J500","J501","N101","N200","N202","N203","N300","N302","N451","N452","S400"]
# ensList = ["D251"] 

ensInfo = EnsInfo.(ensList)

# We do not have charm or disconnected data for some of the ensembles
ensNOcharm = ["J501","N451","D150","D451","J304","C102","D251","D201","J306","J307","F300","H200","N452"]
ensNOdisc  = ["D251","J306","J307","F300","D450"]


@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> 2D HVP computation (+ BM) [NLOc] <==========================##

diag = "NLOc"  # NLOc
wind = "ID"  # NW      to add: (SD  ID  LD  ILD)

IMPR_SET = ["1","2"]

@info("STARTING HVP COMPUTATION [diag. NLOc; wind. NW]")
STD_DERIV ? @info("SRANDART DERIVATIVE is being employed in the IMPROVEMENT") : nothing

@time begin
    for ens in ensInfo

        @info("Computing for ensemble $(ens.id)")

        println("   - Reading t0...")

        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_t0"),"r")
        t0 = uwreal(0.0)
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            t0 = ALPHAdobs_read_next(fb, size=sz, keys=ks)["t0"][1]
        end

        println("   - Reading TMR...")

        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_TMR_NLO"),"r") 
        TMRDict = Dict{String, Any}()
        partial_res = Vector{Dict}()
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            push!(partial_res,ALPHAdobs_read_next(fb, size=sz, keys=ks))
        end
        BDIO_close!(fb)
        for dict in partial_res
            merge!(TMRDict, dict)
        end
        TMR = TMRDict["TMRc"]

        aens = sqrtt0_ph / sqrt(t0)
            
        sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1)
        t = collect(1:sym_points)
        tfm = aens.*(t.-1)

        # define windowed kernel
        if wind == "NW" 
            TMRw = TMR
        else
            TMRw = TMR .* Window2D(wind)(tfm,tfm)
        end

        for impr_set in IMPR_SET

            println("   - Starting set "*impr_set)
            println("      - Reading corr...")

            fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_corr_set"*impr_set),"r")
            corr = Dict{String, Array{uwreal}}()
            while ALPHAdobs_next_p(fb)
                d = ALPHAdobs_read_parameters(fb)
                sz = tuple(d["size"]...)
                ks = collect(d["keys"])
                corr = ALPHAdobs_read_next(fb, size=sz, keys=ks)
            end
            BDIO_close!(fb)

            if ens.kappa_l == ens.kappa_s
                println("      - SU(3) flavour sym point!")

                corr["g88_ll"] = corr["g88_ll_conn"]
                corr["g88_lc"] = corr["g88_lc_conn"]
            else
                corr["g88_ll"] = corr["g88_ll_conn"] .+ corr["g88_ll_disc"] .+ (2).*corr["g08_ll_conn"] .+ corr["g08_ll_disc"] .+ corr["g80_ll_disc"]
                corr["g88_lc"] = corr["g88_lc_conn"] .+ corr["g88_lc_disc"] .+ corr["g08_lc_conn"] .+ corr["g08_lc_disc"]
            end

            println("      - Computing HVPs...")

            tcut0 = 10
            tstep = 1
            tcut_fix = 1.2  # fm

            mpi = m_ens[ens.id]["m_pi"] 
            mrho = m_ens[ens.id]["m_rho"]
            L = ens.L
            E2pi = 2*sqrt(mpi^2 + (2*pi/L)^2)
            E3pi = 2*sqrt(mpi^2 + (2*pi/L)^2) + sqrt(mpi^2 + 2(2*pi/L)^2)

            HVP     =  Dict{String, uwreal}()
            HVPsyst =  Dict{String, Float64}()
            plateau =  Dict{String, Vector{Float64}}()

            for discr in ["ll","lc"] # ["ll-ll","ll-lc","lc-lc"] or ["ll","lc","mixed"]
                println("         - discr. = $discr")

                obs33 = corr["g33_$discr"]
                obs88 = corr["g88_$discr"]

                int3333 = (obs33[t].*hcat(obs33[t]...)) .* TMRw
                int8888 = (obs88[t].*hcat(obs88[t]...)) .* TMRw
                int3388 = (obs33[t].*hcat(obs88[t]...)) .* TMRw

                ub3333 = Vector{uwreal}(); lb3333 = Vector{uwreal}()
                ub8888 = Vector{uwreal}(); lb8888 = Vector{uwreal}()
                ub3388 = Vector{uwreal}(); lb3388 = Vector{uwreal}()

                println("            - Computing bounded correlators...")

                Eeff33 = uwreal(0.0); Eeff88 = uwreal(0.0)
                if ens.kappa_l == ens.kappa_s
                    for tcut in ProgressBar(tcut0:tstep:t[end-1])

                        if aens*tcut < tcut_fix  # we fix the eff energy at some point
                            Eeff33=Eeff(tcut, obs33)
                        end
                        UB33 = mrho < E2pi ? corrBound(t, tcut, obs33, ens, mrho) : corrBound(t, tcut, obs33, ens, E2pi)
                        LB33 = corrBound(t, tcut, obs33, ens, Eeff33)

    
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
                    for tcut in ProgressBar(tcut0:tstep:t[end-1])

                        if aens*tcut < tcut_fix  # we fix the eff energy at some point
                            Eeff33=Eeff(tcut, obs33)
                            Eeff88=Eeff(tcut, obs88)
                        end
                        UB33 = mrho < E2pi ? corrBound(t, tcut, obs33, ens, mrho) : corrBound(t, tcut, obs33, ens, E2pi)
                        LB33 = corrBound(t, tcut, obs33, ens, Eeff33)
                        UB88 = mrho < E3pi ? corrBound(t, tcut, obs88, ens, mrho) : corrBound(t, tcut, obs88, ens, E3pi)
                        LB88 = corrBound(t, tcut, obs88, ens, Eeff88)

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

                        ub_amuNLO3333 = (alpha/pi)^3 * (sum(int3333[1:tcut,1:tcut]) + sum(UBInt3333_1) + 2*sum(UBInt3333_2)) * 1e10
                        ub_amuNLO8888 = (alpha/pi)^3 * (sum(int8888[1:tcut,1:tcut]) + sum(UBInt8888_1) + 2*sum(UBInt8888_2)) * 1e10
                        ub_amuNLO3388 = (alpha/pi)^3 * (sum(int3388[1:tcut,1:tcut]) + sum(UBInt3388_1) + sum(UBInt3388_2) + sum(UBInt8833_2)) * 1e10
            
                        lb_amuNLO3333 = (alpha/pi)^3 * (sum(int3333[1:tcut,1:tcut]) + sum(LBInt3333_1) + 2*sum(LBInt3333_2)) * 1e10
                        lb_amuNLO8888 = (alpha/pi)^3 * (sum(int8888[1:tcut,1:tcut]) + sum(LBInt8888_1) + 2*sum(LBInt8888_2)) * 1e10
                        lb_amuNLO3388 = (alpha/pi)^3 * (sum(int3388[1:tcut,1:tcut]) + sum(LBInt3388_1) + sum(LBInt3388_2) + sum(LBInt8833_2)) * 1e10
            
                        push!(ub3333, ub_amuNLO3333); push!(lb3333, lb_amuNLO3333)
                        push!(ub8888, ub_amuNLO8888); push!(lb8888, lb_amuNLO8888)
                        push!(ub3388, ub_amuNLO3388); push!(lb3388, lb_amuNLO3388)
                    end
                end

                ub = [ub3333,ub8888,ub3388]
                lb = [lb3333,lb8888,lb3388]

                println("            - Applying BM...")

                for (i,comp) in enumerate(["3333","8888","3388"])
                    averb = (ub[i].+lb[i])./2; uwerr.(averb)
                    x0    = findfirst(abs.(value.(ub[i]).-value.(lb[i])) .< 0.75.*err.(averb))
                    xend_x0 = findfirst(abs.(averb[x0:end].-averb[x0]) .> 0.5*err(averb[x0]))
                    xend_ = isnothing(xend_x0) ? x0+15  : x0-1 + xend_x0
                    if xend_-x0 > 5
                        xend = xend_>length(averb) ? length(averb) : xend_
                    else
                        xend=x0+3
                        x0 > 3 ? x0-=3 : x0=1
                    end
                    plateau_fm = value(aens).*(collect(x0:xend).+tcut0.-2)
                    amu = sum(averb[x0:xend])/length(averb[x0:xend])
            
                    aux1 = sum(averb[x0:xend].^2)/length(averb[x0:xend])
                    aux2 = amu^2
                    syst = sqrt(abs(value(aux1 - aux2)))

                    HVP["g$(comp)_$(discr)"]     = amu
                    HVPsyst["g$(comp)_$(discr)"] = syst
                    plateau["g$(comp)_$(discr)"] = [plateau_fm[1],plateau_fm[end]]
                end

                obsCC  = corr["gcc_$(discr)_conn_sim"]
                obsCCp = corr["gcc_$(discr)_conn_sim+"]

                HVP["gCCCC_$(discr)_sim"]  = (alpha/pi)^3 * sum((obsCC[t] .* hcat(obsCC[t]...)) .* TMRw) * 1e10
                HVP["gCCCC_$(discr)_sim+"] = (alpha/pi)^3 * sum((obsCCp[t] .* hcat(obsCCp[t]...)) .* TMRw) * 1e10
                HVP["g33CC_$(discr)_sim"]  = (alpha/pi)^3 * sum((obs33[t] .* hcat(obsCC[t]...)) .* TMRw) * 1e10
                HVP["g33CC_$(discr)_sim+"] = (alpha/pi)^3 * sum((obs33[t] .* hcat(obsCCp[t]...)) .* TMRw) * 1e10
                HVP["g88CC_$(discr)_sim"]  = (alpha/pi)^3 * sum((obs88[t] .* hcat(obsCC[t]...)) .* TMRw) * 1e10
                HVP["g88CC_$(discr)_sim+"] = (alpha/pi)^3 * sum((obs88[t] .* hcat(obsCCp[t]...)) .* TMRw) * 1e10
            end

            println("      - Writing BDIO & JDL2...")
            
            HVPinfo = Dict{String,Dict}(
                "HVPsyst" => HVPsyst,
                "plateau" => plateau,
            )

            pens = joinpath(path_bdio,"HVP&FVC",wind,ens.id)
            !ispath(pens) ? mkdir(pens) : nothing

            pHVP = joinpath(pens,"HVP")
            !ispath(pHVP) ? mkdir(pHVP) : nothing

            io = IOBuffer()
            write(io, "$(ens.id) HVP")
            fb = ALPHAdobs_create(joinpath(pHVP,"$(ens.id)_HVP$(diag)_set$(impr_set)"), io)

            extra = Dict{String, Any}("ens" => ens.id, "impr_set" => impr_set, "diag" => diag, "wind" => wind)
            ALPHAdobs_write(fb, HVP, extra=extra)

            ALPHAdobs_close(fb)

            pinfo = joinpath(pHVP,"$(ens.id)_HVP$(diag)_BMinfo_set$(impr_set).jld2")
            save(pinfo,"HVPinfo",HVPinfo)

            GC.gc() # call garvage collector

        end # end impr_set loop
    end # end ens loop
end # end timer


## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> 2D FVC computation [NLOc] <==========================##

diag = "NLOc"  # NLOc
wind = "NW"    # NW      to add: (SD  ID  LD  ILD)

OVERWRITE = false  # Allows to erase data and overwrite it with new data, use carefully !!

@info("STARTING FVC COMPUTATION [diag. $diag; wind. $wind]")
STD_DERIV ? @info("SRANDART DERIVATIVE is being employed in the IMPROVEMENT") : nothing

@time begin
    for ens in ensInfo

        @info("Computing for ensemble $(ens.id)")

        println("   - Reading TMR...")

        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_TMR_NLO"),"r") 
        TMR = Dict{String, Any}()
        partial_res = Vector{Dict}()
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            push!(partial_res,ALPHAdobs_read_next(fb, size=sz, keys=ks))
        end
        BDIO_close!(fb)
        for dict in partial_res
            merge!(TMR, dict)
        end
        TMR = TMR["TMRc"]

        aens = sqrtt0_ph / sqrt(t0)
            
        sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1)
        t = collect(1:sym_points)
        tfm = aens.*(t.-1)

        # define windowed kernel
        if wind == "NW" 
            TMRw = TMR
        else
            TMRw = TMR .* Window2D(wind)(tfm,tfm)
        end

        println("   - Reading FVC...")

        FVC_corr = Dict{String,Vector{uwreal}}()
        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_FVC"),"r") 
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            FVC_corr = ALPHAdobs_read_next(fb, size=sz, keys=ks)
        end
        BDIO_close!(fb)

        println("   - Computing FVC...")

        for impr_set in ["1","2"]

            println("      - Starting set "*impr_set)
            println("         - Reading corr...")

            fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_corr_set"*impr_set),"r")
            corr = Dict{String, Array{uwreal}}()
            while ALPHAdobs_next_p(fb)
                d = ALPHAdobs_read_parameters(fb)
                sz = tuple(d["size"]...)
                ks = collect(d["keys"])
                corr = ALPHAdobs_read_next(fb, size=sz, keys=ks)
            end
            BDIO_close!(fb)

            tFVC = 1:Int64(length(FVC_corr["FVCPi1"])/2)

            DeltaG33 = Vector{Vector{uwreal}}()
            DeltaG88 = Vector{Vector{uwreal}}()
            if ens.kappa_l == ens.kappa_s
                println("         - SU(3) flavour sym point!")

                corr["g88_ll"] = corr["g88_ll_conn"]
                corr["g88_lc"] = corr["g88_lc_conn"]

                for i in collect(1:6)
                    fullFVCPi = vcat(TMR[1,1],-FVC_corr["FVCPi$i"][tFVC])
                    fullFVCK  = vcat(TMR[1,1],-FVC_corr["FVCK$i" ][tFVC])

                    push!(DeltaG33, (3/2).*fullFVCPi)
                    push!(DeltaG88, (3/2).*fullFVCPi)
                end
            else
                corr["g88_ll"] = corr["g88_ll_conn"] .+ corr["g88_ll_disc"] .+ (2).*corr["g08_ll_conn"] .+ corr["g08_ll_disc"] .+ corr["g80_ll_disc"]
                corr["g88_lc"] = corr["g88_lc_conn"] .+ corr["g88_lc_disc"] .+ corr["g08_lc_conn"] .+ corr["g08_lc_disc"]

                for i in collect(1:6)
                    fullFVCPi = vcat(TMR[1,1],-FVC_corr["FVCPi$i"][tFVC])
                    fullFVCK  = vcat(TMR[1,1],-FVC_corr["FVCK$i" ][tFVC])

                    push!(DeltaG33, fullFVCPi .+ fullFVCK)
                    push!(DeltaG88, (2/3).*fullFVCK)
                end
            end

            dataΔG = Dict{String, Array{uwreal}}()

            println("         - Computing FVC...")


            for discr in ["ll","lc"]
                println("            - discr. $discr")

                dataΔG["FVC3333_$discr"] = []
                dataΔG["FVC8888_$discr"] = []
                dataΔG["FVC3388_$discr"] = []
                dataΔG["FVC33CC_$discr"] = []
                dataΔG["FVC88CC_$discr"] = []
                for i in ProgressBar(collect(1:6))
                    DeltaG8888 = TMR .* (corr["g88_$discr"][t].*hcat(DeltaG88[i]...) .+ DeltaG88[i].*hcat(corr["g88_$discr"][t]...) .+ DeltaG88[i].*hcat(DeltaG88[i]...))
                    DeltaG3388 = TMR .* (corr["g33_$discr"][t].*hcat(DeltaG88[i]...) .+ DeltaG88[i].*hcat(corr["g33_$discr"][t]...) .+ DeltaG33[i].*hcat(DeltaG88[i]...))
                    DeltaG3333 = TMR .* (corr["g33_$discr"][t].*hcat(DeltaG33[i]...) .+ DeltaG33[i].*hcat(corr["g33_$discr"][t]...) .+ DeltaG33[i].*hcat(DeltaG33[i]...))
                    DeltaG33CC = TMR .* (DeltaG33[i].*hcat(corr["gcc_$(discr)_conn"][t]...))
                    DeltaG88CC = TMR .* (DeltaG88[i].*hcat(corr["gcc_$(discr)_conn"][t]...))

                    push!(dataΔG["FVC3333_$discr"],(alpha/pi)^3 * sum(DeltaG3333) * 1e10)
                    push!(dataΔG["FVC8888_$discr"],(alpha/pi)^3 * sum(DeltaG8888) * 1e10)
                    push!(dataΔG["FVC3388_$discr"],(alpha/pi)^3 * sum(DeltaG3388) * 1e10)
                    push!(dataΔG["FVC33CC_$discr"],(alpha/pi)^3 * sum(DeltaG33CC) * 1e10)
                    push!(dataΔG["FVC88CC_$discr"],(alpha/pi)^3 * sum(DeltaG88CC) * 1e10)
                end
            end

            println("      - Writing BDIO...")

            pens = joinpath(path_bdio,"HVP&FVC",wind,ens.id)
            !ispath(pens) ? mkdir(pens) : nothing

            pFVC = joinpath(pens,"FVC")
            !ispath(pFVC) ? mkdir(pFVC) : nothing
            
            pBDIO = joinpath(pFVC,"$(ens.id)_FVC$(diag)_set$impr_set")
            if ispath(pBDIO)
                OVERWRITE ? rm(pBDIO, recursive=true) : error("This information already exist; set 'OVERWRITE = true' to overwrite")
            end

            io = IOBuffer()
            write(io, "$(ens.id) FVC")
            fb = ALPHAdobs_create(pBDIO, io)
        
            extra = Dict{String, Any}("Ens" => ens.id, "Diag" => diag)
            ALPHAdobs_write(fb, dataΔG, extra=extra)

            ALPHAdobs_close(fb)
        end
    end
end # end timer


## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##


##==========================> READING TEST <==========================##

data     = "HVP"   # "HVP"  "FVC"
wind     = "ID"    # "NW"  "SD"  "ID"  "LD"  "ILD"
ensid    = "H102"
diag     = "LO"    # "LO"  "NLOa"  "NLOb"  "NLOc"
impr_set = "1"


hvp, info = BDIOread_HVPens(path_bdio,diag,wind,ensid,impr_set,info=true)

fvc = BDIOread_FVCens(path_bdio,diag,wind,ensid)


##==> Result computation in the 'flavour basis'

ensid    = "H101"
impr_set = "1"
discr    = "ll"

HVP = BDIOread_HVPens(path_bdio,diag,wind,ensid,impr_set,info=false)

light   = 5/9*(2*HVP["g33_$(discr)"]); uwerr(light)
strange = 1/9*(3*HVP["g88_$(discr)_conn"]-HVP["g33_$(discr)"]); uwerr(strange)
charm   = 4/9*HVP["gcc_$(discr)_conn"]; uwerr(charm)

println("Ensemble: $ensid [impr. set $impr_set & discr. $discr]")
println("- amu(l) = $light")
println("- amu(s) = $strange")
println("- amu(c) = $charm")

