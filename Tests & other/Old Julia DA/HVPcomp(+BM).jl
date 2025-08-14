# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

using ProgressBars
using TimerOutputs

# Include necessary functions

# include("amuNLO.jl")
# export amuHVPNLO, amu∆G

include("data_management.jl")

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302"]
ensInfo = EnsInfo.(ensList)

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> HVP computation + BM [a & b] <==========================##

tcut0 = 10
tstep = 1

tcut_fix = 1.2  # fm

function windowfunc(t::Union{Int64,Float64,uwreal},tstar::Union{Int64,Float64},Delta::Union{Int64,Float64})
    return 0.5*(1 + tanh((t-tstar)/Delta))
end
Delta = 1.5

@time begin
    for ens in ensInfo

        @info("Computing HVP for diag. 'a' & 'b'; ensemble: $(ens.id)")

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

        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_TMR"),"r") 
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
        TMRa = TMR["TMRa"]
        TMRb = TMR["TMRb"]

        for impr_set in ["1","2"]

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

            sym_points = Int64(length(corr["g33_ll"])/2+1)
            t = collect(1:sym_points)
            tcut0 = 10

            aens = t0_ph / sqrt(t0)

            mpi = m_ens[ens.id]["m_pi"] 
            mrho = m_ens[ens.id]["m_rho"]
            L = ens.L
            E2pi = 2*sqrt(mpi^2 + (2*pi/L)^2)
            E3pi = 2*sqrt(mpi^2 + (2*pi/L)^2) + sqrt(mpi^2 + 2(2*pi/L)^2)

            HVP     =  Dict{String, Array{uwreal}}()
            HVPsyst =  Dict{String, Vector{Float64}}()
            plateau =  Dict{String, Vector{Vector{Float64}}}()

            for key in ["g33_ll","g33_lc","g88_ll","g88_lc"]
                println("         - $key")
                
                ub_a = Vector{uwreal}()
                ub_b = Vector{uwreal}()

                lb_a = Vector{uwreal}()
                lb_b = Vector{uwreal}()

                obs = corr[key]
                inta = obs[t] .* TMRa
                intb = obs[t] .* TMRb

                Eeff_ = uwreal(0.0)
                for tcut in tcut0:t[end-1]

                    if aens*tcut < tcut_fix  # we fix the eff energy at some point
                        Eeff_=Eeff(tcut, obs)
                    end
                    if key[2:3] == "33"
                        UB = mrho < E2pi ? Gbounding(t, tcut, obs, ens, mrho) : Gbounding(t, tcut, obs, ens, E2pi)
                    elseif key[2:3] == "88"
                        UB = Gbounding(t, tcut, obs, ens, mrho)
                    end
                    # print("$(value(Eeff_))\n")
                    LB = Gbounding(t, tcut, obs, ens, Eeff_)
            
                    UBInta = UB .* TMRa[tcut+1:end]
                    UBIntb = UB .* TMRb[tcut+1:end]

                    LBinta = LB .* TMRa[tcut+1:end]
                    LBintb = LB .* TMRb[tcut+1:end]
            
                    ub_amuNLOa = (alpha/pi)^3 * (sum(inta[1:tcut])+sum(UBInta)) * 1e10
                    ub_amuNLOb = (alpha/pi)^3 * (sum(intb[1:tcut])+sum(UBIntb)) * 1e10

                    lb_amuNLOa = (alpha/pi)^3 * (sum(inta[1:tcut])+sum(LBinta)) * 1e10
                    lb_amuNLOb = (alpha/pi)^3 * (sum(intb[1:tcut])+sum(LBintb)) * 1e10

                    push!(ub_a, ub_amuNLOa)
                    push!(ub_b, ub_amuNLOb)
        
                    push!(lb_a, lb_amuNLOa)
                    push!(lb_b, lb_amuNLOb)
                end

                ub = [ub_a,ub_b]
                lb = [lb_a,lb_b]

                HVP[key] = Vector{uwreal}(); HVPsyst[key] = Vector{Float64}(); plateau[key] = Vector{Vector{Float64}}()
                for i in [1,2]
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

                    push!(HVP[key], amu)
                    push!(HVPsyst[key], syst)
                    push!(plateau[key], [plateau_fm[1],plateau_fm[end]])
                end
            end

            for key in ["gcc_ll_conn","gcc_lc_conn"]
                println("         - $key")

                obs = corr[key][t]

                inta = obs .* TMRa
                intb = obs .* TMRb
                
                amua = (alpha/pi)^3 * sum(inta) * 1e10
                amub = (alpha/pi)^3 * sum(intb) * 1e10

                HVP[key] = [amua,amub]
            end


            if ens.kappa_l != ens.kappa_s
                for key in ["gcc_cc_disc","gc8_cc_disc"]
                    println("         - $key")

                    obs = corr[key][t]; uwerr.(obs)

                    tstar = findfirst(abs.(err.(obs[2:end])./value.(obs[2:end])) .> 0.5) + 1
                    windvec = 1 .- windowfunc.(t,tstar,Delta)

                    intWa = obs .* windvec .* TMRa
                    intWb = obs .* windvec .* TMRb

                    amua = (alpha/pi)^3 * sum(intWa) * 1e10
                    amub = (alpha/pi)^3 * sum(intWb) * 1e10

                    HVP[key]     = [amua,amub]
                    HVPsyst[key] = 0.5.*value.([amua,amub])
                end
            end

            println("      - Writing BDIO & JDL2...")
            
            HVPinfo = Dict{String,Dict}(
                "HVPsyst" => HVPsyst,
                "plateau" => plateau,
            )

            pens = joinpath(path_bdio,"HVP&FVC",ens.id)
            !ispath(pens) ? mkdir(pens) : nothing

            pHVP = joinpath(pens,"HVP")
            !ispath(pHVP) ? mkdir(pHVP) : nothing

            io = IOBuffer()
            write(io, "$(ens.id) HVP")
            fb = ALPHAdobs_create(joinpath(pHVP,"$(ens.id)_HVP_set$(impr_set)"), io)

            extra = Dict{String, Any}("Ens" => ens.id, "Set" => impr_set)
            ALPHAdobs_write(fb, HVP, extra=extra)

            ALPHAdobs_close(fb)

            pinfo = joinpath(pHVP,"$(ens.id)_HVPinfo_set$(impr_set).jld2")
            save(pinfo,"HVPinfo",HVPinfo)
        end
    end
end # end timer

##==========================> HVP computation + BM [c] <==========================##

tcut0 = 10
tstep = 1
tcut_fix = 1.2  # fm

IMPR_SET = ["1","2"]

@time begin
    for ens in ensInfo

        @info("Computing HVP for diag. 'c' (ensemble $(ens.id))")

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

        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_TMR"),"r") 
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
        TMRc = TMR["TMRc"]

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

            sym_points = Int64(length(corr["g33_ll"])/2+1)
            t = collect(1:sym_points)

            aens = t0_ph / sqrt(t0)

            mpi = m_ens[ens.id]["m_pi"] 
            mrho = m_ens[ens.id]["m_rho"]
            L = ens.L
            E2pi = 2*sqrt(mpi^2 + (2*pi/L)^2)
            E3pi = 2*sqrt(mpi^2 + (2*pi/L)^2) + sqrt(mpi^2 + 2(2*pi/L)^2)

            HVP     =  Dict{String, uwreal}()
            HVPsyst =  Dict{String, Float64}()
            plateau =  Dict{String, Vector{Float64}}()

            for discr in ["ll" ,"lc"]
                println("         - discr. = $discr")

                obs33 = corr["g33_$discr"]
                obs88 = corr["g88_$discr"]

                int3333 = (obs33[t].*hcat(obs33[t]...)) .* TMRc
                int8888 = (obs88[t].*hcat(obs88[t]...)) .* TMRc
                int3388 = (obs33[t].*hcat(obs88[t]...)) .* TMRc

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
                        UB33 = mrho < E2pi ? Gbounding(t, tcut, obs33, ens, mrho) : Gbounding(t, tcut, obs33, ens, E2pi)
                        LB33 = Gbounding(t, tcut, obs33, ens, Eeff33)

    
                        UBInt3333_1 = (UB33 .* hcat(UB33...)) .* TMRc[tcut+1:end,tcut+1:end]
                        UBInt3333_2 = (obs33[1:tcut] .* hcat(UB33...)) .* TMRc[1:tcut,tcut+1:end]
    
                        LBInt3333_1 = (LB33 .* hcat(LB33...)) .* TMRc[tcut+1:end,tcut+1:end]
                        LBInt3333_2 = (obs33[1:tcut] .* hcat(LB33...)) .* TMRc[1:tcut,tcut+1:end]
    
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
                        UB33 = mrho < E2pi ? Gbounding(t, tcut, obs33, ens, mrho) : Gbounding(t, tcut, obs33, ens, E2pi)
                        LB33 = Gbounding(t, tcut, obs33, ens, Eeff33)
                        UB88 = mrho < E3pi ? Gbounding(t, tcut, obs88, ens, mrho) : Gbounding(t, tcut, obs88, ens, E3pi)
                        LB88 = Gbounding(t, tcut, obs88, ens, Eeff88)

                        UBInt3333_1 = (UB33 .* hcat(UB33...)) .* TMRc[tcut+1:end,tcut+1:end]
                        UBInt8888_1 = (UB88 .* hcat(UB88...)) .* TMRc[tcut+1:end,tcut+1:end]
                        UBInt3388_1 = (UB33 .* hcat(UB88...)) .* TMRc[tcut+1:end,tcut+1:end]
                        UBInt3333_2 = (obs33[1:tcut] .* hcat(UB33...)) .* TMRc[1:tcut,tcut+1:end]
                        UBInt8888_2 = (obs88[1:tcut] .* hcat(UB88...)) .* TMRc[1:tcut,tcut+1:end]
                        UBInt3388_2 = (obs33[1:tcut] .* hcat(UB88...)) .* TMRc[1:tcut,tcut+1:end]
                        UBInt8833_2 = (obs88[1:tcut] .* hcat(UB33...)) .* TMRc[1:tcut,tcut+1:end]

                        LBInt3333_1 = (LB33 .* hcat(LB33...)) .* TMRc[tcut+1:end,tcut+1:end]
                        LBInt8888_1 = (LB88 .* hcat(LB88...)) .* TMRc[tcut+1:end,tcut+1:end]
                        LBInt3388_1 = (LB33 .* hcat(LB88...)) .* TMRc[tcut+1:end,tcut+1:end]
                        LBInt3333_2 = (obs33[1:tcut] .* hcat(LB33...)) .* TMRc[1:tcut,tcut+1:end]
                        LBInt8888_2 = (obs88[1:tcut] .* hcat(LB88...)) .* TMRc[1:tcut,tcut+1:end]
                        LBInt3388_2 = (obs33[1:tcut] .* hcat(LB88...)) .* TMRc[1:tcut,tcut+1:end]
                        LBInt8833_2 = (obs88[1:tcut] .* hcat(LB33...)) .* TMRc[1:tcut,tcut+1:end]

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

                obsCC = corr["gcc_$(discr)_conn"]

                HVP["gCCCC_$(discr)"] = (alpha/pi)^3 * sum((obsCC[t] .* hcat(obsCC[t]...)) .* TMRc) * 1e10
                HVP["g33CC_$(discr)"] = (alpha/pi)^3 * sum((obs33[t] .* hcat(obsCC[t]...)) .* TMRc) * 1e10
                HVP["g88CC_$(discr)"] = (alpha/pi)^3 * sum((obs88[t] .* hcat(obsCC[t]...)) .* TMRc) * 1e10

            end

            println("      - Writing BDIO & JDL2...")
            
            HVPinfo = Dict{String,Dict}(
                "HVPsyst" => HVPsyst,
                "plateau" => plateau,
            )

            pens = joinpath(path_bdio,"HVP&FVC",ens.id)
            !ispath(pens) ? mkdir(pens) : nothing

            pHVP = joinpath(pens,"HVP")
            !ispath(pHVP) ? mkdir(pHVP) : nothing

            io = IOBuffer()
            write(io, "$(ens.id) HVP")
            fb = ALPHAdobs_create(joinpath(pHVP,"$(ens.id)_HVPC_set$(impr_set)"), io)

            extra = Dict{String, Any}("Ens" => ens.id, "Set" => impr_set)
            ALPHAdobs_write(fb, HVP, extra=extra)

            ALPHAdobs_close(fb)

            pinfo = joinpath(pHVP,"$(ens.id)_HVPinfoC_set$(impr_set).jld2")
            save(pinfo,"HVPinfo",HVPinfo)
        end
    end
end # end timer

##==========================> HVP computation + BM [LO] <==========================##

tcut0 = 10
tstep = 1

tcut_fix = 1.2  # fm

function windowfunc(t::Union{Int64,Float64,uwreal},tstar::Union{Int64,Float64},Delta::Union{Int64,Float64})
    return 0.5*(1 + tanh((t-tstar)/Delta))
end
Delta = 1.5

@time begin
    for ens in ensInfo

        @info("Computing LO HVP; ensemble: $(ens.id)")

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

        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_TMR_LO"),"r") 
        TMRLO = Dict{String, Any}()
        partial_res = Vector{Dict}()
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            TMRLO = ALPHAdobs_read_next(fb, size=sz, keys=ks)
        end
        BDIO_close!(fb)
        TMR = TMRLO["TMR"]

        for impr_set in ["1","2"]

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

                # corr["g88_ll"] = corr["g88_ll_conn"]
                # corr["g88_lc"] = corr["g88_lc_conn"]

                light_keys = ["g33_ll","g33_lc"]
            else
                corr["g88_ll"] = corr["g88_ll_conn"] .+ corr["g88_ll_disc"] .+ (2).*corr["g08_ll_conn"] .+ corr["g08_ll_disc"] .+ corr["g80_ll_disc"]
                corr["g88_lc"] = corr["g88_lc_conn"] .+ corr["g88_lc_disc"] .+ corr["g08_lc_conn"] .+ corr["g08_lc_disc"]

                light_keys = ["g33_ll","g33_lc","g88_ll_conn","g88_lc_conn","g88_ll","g88_lc"]
            end

            println("      - Computing HVPs...")

            sym_points = Int64(length(corr["g33_ll"])/2+1)
            t = collect(1:sym_points)
            tcut0 = 10

            aens = t0_ph / sqrt(t0)

            mpi = m_ens[ens.id]["m_pi"] 
            mrho = m_ens[ens.id]["m_rho"]
            L = ens.L
            E2pi = 2*sqrt(mpi^2 + (2*pi/L)^2)
            E3pi = 2*sqrt(mpi^2 + (2*pi/L)^2) + sqrt(mpi^2 + 2(2*pi/L)^2)

            HVP     =  Dict{String, uwreal}()
            HVPsyst =  Dict{String, Float64}()
            plateau =  Dict{String, Vector{Float64}}()

            for key in light_keys
                println("         - $key")
                
                ub = Vector{uwreal}()
                lb = Vector{uwreal}()

                obs = corr[key]
                int = obs[t] .* TMR

                Eeff_ = uwreal(0.0)
                for tcut in tcut0:t[end-1]

                    if aens*tcut < tcut_fix  # we fix the eff energy at some point
                        Eeff_=Eeff(tcut, obs)
                    end
                    if key[2:3] == "33"
                        UB = mrho < E2pi ? Gbounding(t, tcut, obs, ens, mrho) : Gbounding(t, tcut, obs, ens, E2pi)
                    elseif key[2:3] == "88"
                        UB = Gbounding(t, tcut, obs, ens, mrho)
                    end
                    # print("$(value(Eeff_))\n")
                    LB = Gbounding(t, tcut, obs, ens, Eeff_)
            
                    UBInt = UB .* TMR[tcut+1:end]
                    LBint = LB .* TMR[tcut+1:end]
            
                    ub_amuNLO = (alpha/pi)^2 * (sum(int[1:tcut])+sum(UBInt)) * 1e10
                    lb_amuNLO = (alpha/pi)^2 * (sum(int[1:tcut])+sum(LBint)) * 1e10

                    push!(ub, ub_amuNLO)
                    push!(lb, lb_amuNLO)
                end

                averb = (ub.+lb)./2; uwerr.(averb)
                x0    = findfirst(abs.(value.(ub).-value.(lb)) .< 0.75.*err.(averb))
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

                HVP[key] = amu
                HVPsyst[key] = syst
                plateau[key] = [plateau_fm[1],plateau_fm[end]]
            end

            if ens.kappa_l == ens.kappa_s
                HVP["g88_ll_conn"] = HVP["g88_ll"] = HVP["g33_ll"]
                HVP["g88_lc_conn"] = HVP["g88_lc"] = HVP["g33_lc"]

                HVPsyst["g88_ll_conn"] = HVPsyst["g88_ll"] = HVPsyst["g33_ll"]
                HVPsyst["g88_lc_conn"] = HVPsyst["g88_lc"] = HVPsyst["g33_lc"]

                plateau["g88_ll_conn"] = plateau["g88_ll"] = plateau["g33_ll"]
                plateau["g88_lc_conn"] = plateau["g88_lc"] = plateau["g33_lc"]
            end


            for key in ["gcc_ll_conn","gcc_lc_conn"]
                println("         - $key")

                obs = corr[key][t]

                int = obs .* TMR
                
                amu = (alpha/pi)^2 * sum(int) * 1e10

                HVP[key] = amu
            end


            if ens.kappa_l != ens.kappa_s
                for key in ["gcc_cc_disc","gc8_cc_disc"]
                    println("         - $key")

                    obs = corr[key][t]; uwerr.(obs)

                    tstar = findfirst(abs.(err.(obs[2:end])./value.(obs[2:end])) .> 0.5) + 1
                    windvec = 1 .- windowfunc.(t,tstar,Delta)

                    intW = obs .* windvec .* TMR

                    amu = (alpha/pi)^2 * sum(intW) * 1e10

                    HVP[key]     = amu
                    HVPsyst[key] = 0.5.*value(amu)
                end
            end

            println("      - Writing BDIO & JDL2...")
            
            HVPinfo = Dict{String,Dict}(
                "HVPsyst" => HVPsyst,
                "plateau" => plateau,
            )

            pens = joinpath(path_bdio,"HVP&FVC",ens.id)
            !ispath(pens) ? mkdir(pens) : nothing

            pHVP = joinpath(pens,"HVP")
            !ispath(pHVP) ? mkdir(pHVP) : nothing

            io = IOBuffer()
            write(io, "$(ens.id) HVP")
            fb = ALPHAdobs_create(joinpath(pHVP,"$(ens.id)_HVPLO_set$(impr_set)"), io)

            extra = Dict{String, Any}("Ens" => ens.id, "Set" => impr_set)
            ALPHAdobs_write(fb, HVP, extra=extra)

            ALPHAdobs_close(fb)

            pinfo = joinpath(pHVP,"$(ens.id)_HVPLOinfo_set$(impr_set).jld2")
            save(pinfo,"HVPinfo",HVPinfo)
        end
    end
end # end timer


## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> FVC computation [NLO] <==========================##

FVCab = false
FVCc  = true

IMPR_SET = ["1","2"]

if FVCab && FVCc
    DIAGstring = "'a', 'b' & 'c'"
elseif FVCab
    DIAGstring = "'a' & 'b'"
elseif FVCc
    DIAGstring = "'c'"
end

@time begin
    for ens in ensInfo

        @info("Computing FVCs (ensemble $(ens.id); NLO diag. $DIAGstring; impr. set $IMPR_SET)")

        println("   - Reading TMR...")

        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_TMR"),"r") 
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
        TMRa = TMR["TMRa"]
        TMRb = TMR["TMRb"]
        TMRc = TMR["TMRc"]

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

        if FVCab
            println("   - Computing FVC for 'a' & 'b'...")

            mult = ens.kappa_l == ens.kappa_s ? 1.5 : 1.0

            FVCPi_a = Vector{uwreal}()
            FVCPi_b = Vector{uwreal}()
            FVCK_a  = Vector{uwreal}()
            FVCK_b  = Vector{uwreal}()

            tFVC = 1:Int64(length(FVC_corr["FVCPi1"])/2)
            for i in collect(1:6)
                fullFVCPi = vcat(TMRa[1],-FVC_corr["FVCPi$i"][tFVC])
                fullFVCK  = vcat(TMRa[1],-FVC_corr["FVCK$i" ][tFVC])
                # print(length(fullFVCPi))
                # print(length(TMRa))

                FVCPiinta = mult .* fullFVCPi .* TMRa
                FVCPiintb = mult .* fullFVCPi .* TMRb

                FVCKinta  = mult .* fullFVCK .* TMRa
                FVCKintb  = mult .* fullFVCK .* TMRb

                push!(FVCPi_a,(alpha/pi)^3 * sum(FVCPiinta) * 1e10)
                push!(FVCPi_b,(alpha/pi)^3 * sum(FVCPiintb) * 1e10)
                push!(FVCK_a ,(alpha/pi)^3 * sum(FVCKinta ) * 1e10)
                push!(FVCK_b ,(alpha/pi)^3 * sum(FVCKintb ) * 1e10)
            end

            matrixPi = hcat([FVCPi_a,FVCPi_b]...)
            matrixK  = hcat([FVCK_a ,FVCK_b ]...)

            dataFVC_ab =  Dict{String, Array{uwreal}}("FVCPi" => matrixPi, "FVCK" => matrixK)

            println("      - Writing BDIO ('a' & 'b')...")

            pens = joinpath(path_bdio,"HVP&FVC",ens.id)
            !ispath(pens) ? mkdir(pens) : nothing

            pFVC = joinpath(pens,"FVC")
            !ispath(pFVC) ? mkdir(pFVC) : nothing

            io = IOBuffer()
            write(io, "$(ens.id) FVC a & b")
            fb = ALPHAdobs_create(joinpath(pFVC,"$(ens.id)_FVC_ab"), io)

            extra = Dict{String, Any}("Ens" => ens.id, "Diag" => "a&b")
            ALPHAdobs_write(fb, dataFVC_ab, extra=extra)

            ALPHAdobs_close(fb)
        end
        if FVCc
            println("   - Computing FVC for 'c'...")

            for impr_set in IMPR_SET

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

                sym_points = Int64(length(corr["g33_ll"])/2+1)
                t = collect(1:sym_points)
                tFVC = 1:Int64(length(FVC_corr["FVCPi1"])/2)

                DeltaG33 = Vector{Vector{uwreal}}()
                DeltaG88 = Vector{Vector{uwreal}}()
                if ens.kappa_l == ens.kappa_s
                    println("         - SU(3) flavour sym point!")
    
                    corr["g88_ll"] = corr["g88_ll_conn"]
                    corr["g88_lc"] = corr["g88_lc_conn"]

                    for i in collect(1:6)
                        fullFVCPi = vcat(TMRa[1],-FVC_corr["FVCPi$i"][tFVC])
                        fullFVCK  = vcat(TMRa[1],-FVC_corr["FVCK$i" ][tFVC])
    
                        push!(DeltaG33, (3/2).*fullFVCPi)
                        push!(DeltaG88, (3/2).*fullFVCPi)
                    end
                else
                    corr["g88_ll"] = corr["g88_ll_conn"] .+ corr["g88_ll_disc"] .+ (2).*corr["g08_ll_conn"] .+ corr["g08_ll_disc"] .+ corr["g80_ll_disc"]
                    corr["g88_lc"] = corr["g88_lc_conn"] .+ corr["g88_lc_disc"] .+ corr["g08_lc_conn"] .+ corr["g08_lc_disc"]

                    for i in collect(1:6)
                        fullFVCPi = vcat(TMRa[1],-FVC_corr["FVCPi$i"][tFVC])
                        fullFVCK  = vcat(TMRa[1],-FVC_corr["FVCK$i" ][tFVC])
    
                        push!(DeltaG33, fullFVCPi .+ fullFVCK)
                        push!(DeltaG88, (2/9).*fullFVCK)
                    end
                end

                dataΔG = Dict{String, Array{uwreal}}()

                println("         - Computing...")


                for discr in ["ll","lc"]
                    println("            - discr. $discr")

                    dataΔG["g3333_$discr"] = []
                    dataΔG["g8888_$discr"] = []
                    dataΔG["g3388_$discr"] = []
                    dataΔG["g33CC_$discr"] = []
                    dataΔG["g88CC_$discr"] = []
                    for i in ProgressBar(collect(1:6))
                        DeltaG8888 = TMRc .* (corr["g88_$discr"][t].*hcat(DeltaG88[i]...) .+ DeltaG88[i].*hcat(corr["g88_$discr"][t]...) .+ DeltaG88[i].*hcat(DeltaG88[i]...))
                        DeltaG3388 = TMRc .* (corr["g33_$discr"][t].*hcat(DeltaG88[i]...) .+ DeltaG88[i].*hcat(corr["g33_$discr"][t]...) .+ DeltaG33[i].*hcat(DeltaG88[i]...))
                        DeltaG3333 = TMRc .* (corr["g33_$discr"][t].*hcat(DeltaG33[i]...) .+ DeltaG33[i].*hcat(corr["g33_$discr"][t]...) .+ DeltaG33[i].*hcat(DeltaG33[i]...))
                        DeltaG33CC = TMRc .* (DeltaG33[i].*hcat(corr["gcc_$(discr)_conn"][t]...))
                        DeltaG88CC = TMRc .* (DeltaG88[i].*hcat(corr["gcc_$(discr)_conn"][t]...))

                        push!(dataΔG["g3333_$discr"],(alpha/pi)^3 * sum(DeltaG3333) * 1e10)
                        push!(dataΔG["g8888_$discr"],(alpha/pi)^3 * sum(DeltaG8888) * 1e10)
                        push!(dataΔG["g3388_$discr"],(alpha/pi)^3 * sum(DeltaG3388) * 1e10)
                        push!(dataΔG["g33CC_$discr"],(alpha/pi)^3 * sum(DeltaG33CC) * 1e10)
                        push!(dataΔG["g88CC_$discr"],(alpha/pi)^3 * sum(DeltaG88CC) * 1e10)
                    end
                end

                println("      - Writing BDIO ('c')...")

                pens = joinpath(path_bdio,"HVP&FVC",ens.id)
                !ispath(pens) ? mkdir(pens) : nothing
    
                pFVC = joinpath(pens,"FVC")
                !ispath(pFVC) ? mkdir(pFVC) : nothing
            
                io = IOBuffer()
                write(io, "$(ens.id) FVC c")
                fb = ALPHAdobs_create(joinpath(pFVC,"$(ens.id)_FVC_c_set$impr_set"), io)
            
                extra = Dict{String, Any}("Ens" => ens.id, "Diag" => "c")
                ALPHAdobs_write(fb, dataΔG, extra=extra)

                ALPHAdobs_close(fb)
            end
        end
    end
end # end timer


##==========================> FVC computation [LO] <==========================##

IMPR_SET = ["1","2"]

@time begin
    for ens in ensInfo

        @info("Computing FVCs (ensemble $(ens.id); LO; impr. set $IMPR_SET)")

        println("   - Reading TMR...")

        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_TMR_LO"),"r") 
        TMRLO = Dict{String, Any}()
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            TMRLO = ALPHAdobs_read_next(fb, size=sz, keys=ks)
        end
        BDIO_close!(fb)
        TMR = TMRLO["TMR"]

        println("   - Reading FVC (corr)...")

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

        mult = ens.kappa_l == ens.kappa_s ? 1.5 : 1.0

        FVCPi = Vector{uwreal}()
        FVCK  = Vector{uwreal}()

        tFVC = 1:Int64(length(FVC_corr["FVCPi1"])/2)
        for i in collect(1:6)
            fullFVCPi = vcat(TMR[1],-FVC_corr["FVCPi$i"][tFVC])
            fullFVCK  = vcat(TMR[1],-FVC_corr["FVCK$i" ][tFVC])

            FVCPiint = mult .* fullFVCPi .* TMR
            FVCKint  = mult .* fullFVCK  .* TMR

            push!(FVCPi,(alpha/pi)^2 * sum(FVCPiint) * 1e10)
            push!(FVCK ,(alpha/pi)^2 * sum(FVCKint ) * 1e10)
        end

        dataFVC =  Dict{String, Array{uwreal}}("FVCPi" => FVCPi, "FVCK" => FVCK)

        println("      - Writing BDIO [LO]...")

        pens = joinpath(path_bdio,"HVP&FVC",ens.id)
        !ispath(pens) ? mkdir(pens) : nothing

        pFVC = joinpath(pens,"FVC")
        !ispath(pFVC) ? mkdir(pFVC) : nothing

        io = IOBuffer()
        write(io, "$(ens.id) FVC LO")
        fb = ALPHAdobs_create(joinpath(pFVC,"$(ens.id)_FVC_LO"), io)

        extra = Dict{String, Any}("Ens" => ens.id, "Diag" => "LO")
        ALPHAdobs_write(fb, dataFVC, extra=extra)

        ALPHAdobs_close(fb)
    end
end # end timer

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##



##==========================> READING TEST <==========================##

ensid = "N300"
impr_set = "1"
extract_data = "HVP LO"   # "HVP LO"  "HVPab"  "HVPc"  "FVC LO"  "FVCab"  "FVCc"

if extract_data == "HVP LO"
    fb = BDIO_open(joinpath(path_bdio,"HVP&FVC",ensid,"HVP","$(ensid)_HVPLO_set$(impr_set)"),"r")
    val = Dict{String, Dict{String, uwreal}}()
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        nobs = d["nobs"]
        dims = d["dimensions"]
        ks = collect(d["keys"])
        val["HVP"] =  ALPHAdobs_read_next(fb, keys=ks)
    end
    BDIO_close!(fb)
    info = load(joinpath(path_bdio,"HVP&FVC",ensid,"HVP","$(ensid)_HVPLOinfo_set$(impr_set).jld2"), "HVPinfo")
    res = merge(val,info)
elseif extract_data == "HVPab"
    fb = BDIO_open(joinpath(path_bdio,"HVP&FVC",ensid,"HVP","$(ensid)_HVP_set$(impr_set)"),"r")
    val = Dict{String, Dict{String, Vector{uwreal}}}()
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        val["HVP"] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
    end
    BDIO_close!(fb)
    info = load(joinpath(path_bdio,"HVP&FVC",ensid,"HVP","$(ensid)_HVPinfo_set$(impr_set).jld2"), "HVPinfo")
    res = merge(val,info)
elseif extract_data == "HVPc"
    fb = BDIO_open(joinpath(path_bdio,"HVP&FVC",ensid,"HVP","$(ensid)_HVPC_set$(impr_set)"),"r")
    val = Dict{String, Dict{String, uwreal}}()
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        nobs = d["nobs"]
        dims = d["dimensions"]
        ks = collect(d["keys"])
        val["HVP"] =  ALPHAdobs_read_next(fb, keys=ks)
    end
    BDIO_close!(fb)
    info = load(joinpath(path_bdio,"HVP&FVC",ensid,"HVP","$(ensid)_HVPinfoC_set$(impr_set).jld2"), "HVPinfo")
    res = merge(val,info)
elseif extract_data == "FVC LO"
    res = Dict{String, Any}()
    fb_ab = BDIO_open(joinpath(path_bdio,"HVP&FVC",ensid,"FVC","$(ensid)_FVC_ab"),"r")
    val_ab = Dict{String, Matrix{uwreal}}()
    while ALPHAdobs_next_p(fb_ab)
        d = ALPHAdobs_read_parameters(fb_ab)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        res = ALPHAdobs_read_next(fb_ab, size=sz, keys=ks)
    end
    BDIO_close!(fb_ab)
elseif extract_data == "FVCab"
    res = Dict{String, Any}()
    fb_ab = BDIO_open(joinpath(path_bdio,"HVP&FVC",ensid,"FVC","$(ensid)_FVC_ab"),"r")
    val_ab = Dict{String, Matrix{uwreal}}()
    while ALPHAdobs_next_p(fb_ab)
        d = ALPHAdobs_read_parameters(fb_ab)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        val_ab = ALPHAdobs_read_next(fb_ab, size=sz, keys=ks)
    end
    BDIO_close!(fb_ab)
    for FVC in ["FVCPi","FVCK"]
        res[FVC] = [val_ab[FVC][:,1],val_ab[FVC][:,2]]
    end
elseif extract_data == "FVCc"
    fb_c = BDIO_open(joinpath(path_bdio,"HVP&FVC",ensid,"FVC","$(ensid)_FVC_c_set$(impr_set)"),"r")
    res = Dict{String, Vector{uwreal}}()
    while ALPHAdobs_next_p(fb_c)
        d = ALPHAdobs_read_parameters(fb_c)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        res = ALPHAdobs_read_next(fb_c, size=sz, keys=ks)
    end
    BDIO_close!(fb_c)
end

res

##==> Result computation in the 'flavour basis'

ensid = "H101"
impr_set = "1"
discr = "lc"

fb = BDIO_open(joinpath(path_bdio,"HVP&FVC",ensid,"HVP","$(ensid)_HVPLO_set$(impr_set)"),"r")
val = Dict{String, Dict{String, uwreal}}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    ks = collect(d["keys"])
    val["HVP"] =  ALPHAdobs_read_next(fb, keys=ks)
end
BDIO_close!(fb)
info = load(joinpath(path_bdio,"HVP&FVC",ensid,"HVP","$(ensid)_HVPLOinfo_set$(impr_set).jld2"), "HVPinfo")
res = merge(val,info)
HVP = res["HVP"]

light   = 5/9*(2*HVP["g33_$(discr)"]); uwerr(light)
strange = 1/9*(3*HVP["g88_$(discr)_conn"]-HVP["g33_$(discr)"]); uwerr(strange)
charm   = 4/9*HVP["gcc_$(discr)_conn"]; uwerr(charm)

println("Ensemble: $ensid & impr. set $impr_set")
println("- amu(l)[$(discr)] = $light")
println("- amu(s)[$(discr)] = $strange")
println("- amu(c)[$(discr)] = $charm")

