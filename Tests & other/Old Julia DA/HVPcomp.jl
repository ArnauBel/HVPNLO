# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO


# Include necessary functions

# include("amuNLO.jl")
# export amuHVPNLO, amu∆G

include("const.jl")

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302"]
ensInfo = EnsInfo.(ensList)

## <<----------------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> HVP computation <==========================##

for ens in ensInfo

    @info("Reading data ensemble: $(ens.id)")

    println("   - Reading TMR")

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

    for impr_set in ["1","2"]

        println("   - Starting set "*impr_set)
        println("      - Reading corr")

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
            keys = ["g33_ll","g33_lc","g88_ll","g88_lc","gcc_ll_conn","gcc_lc_conn"]

            corr["g88_ll"] = corr["g88_ll_conn"]
            corr["g88_lc"] = corr["g88_lc_conn"]
        else
            keys = ["g33_ll","g33_lc","g88_ll","g88_lc","gcc_ll_conn","gcc_lc_conn","gcc_cc_disc","gc8_cc_disc"]

            corr["g88_ll"] = corr["g88_ll_conn"] .+ corr["g88_ll_disc"]
            corr["g88_lc"] = corr["g88_lc_conn"] .+ corr["g88_lc_disc"]

            corr["g08_ll"] = corr["g08_ll_conn"] .+ corr["g08_ll_disc"] .+ corr["g80_ll_disc"]
            corr["g08_lc"] = corr["g08_lc_conn"] .+ corr["g08_lc_disc"]
        end

        sym_points = Int64(length(corr["g33_ll"])/2+1)

        println("      - Computing HVPs")

        data_HVP =  Dict{String, Array{uwreal}}()
        for key in keys
            obs = corr[key][1:sym_points]

            inta = obs .* TMRa
            intb = obs .* TMRb
            intc = (obs .* hcat(obs...)) .* TMRc
            
            amua = (alpha/pi)^3 * sum(inta) * 1e10
            amub = (alpha/pi)^3 * sum(intb) * 1e10
            amuc = (alpha/pi)^3 * sum(intc) * 1e10

            data_HVP[key] = [amua,amub,amuc]
        end

        println("      - Writing BDIO")

        # Create path if it does not exist
        pens = joinpath(path_bdio,"HVP",ens.id)
        !ispath(pens) ? mkdir(pens) : nothing
            
        jo = IOBuffer()
        write(jo, "$(ens.id)  amu HVP, improvement set $(impr_set)")
        fb = ALPHAdobs_create(joinpath(pens,"$(ens.id)_HVP_set$(impr_set)"), jo)

        extra = Dict{String, Any}("Ens" => ens.id, "Impr_Set" => impr_set)
        ALPHAdobs_write(fb, data_HVP, extra=extra)
            
        ALPHAdobs_close(fb)
    end
end

##==========================> FVC computation <==========================##

for ens in ensInfo

    @info("Reading data ensemble: $(ens.id)")

    println("   - Reading TMR")

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

    println("   - Computing FVC for 'a' and 'b'...")

    mult = ens.kappa_l == ens.kappa_s ? 1.5 : 1.0

    FVCPi_a = Vector{uwreal}()
    FVCPi_b = Vector{uwreal}()
    FVCK_a = Vector{uwreal}()
    FVCK_b = Vector{uwreal}()
    for i in 1:length(FVC_corr)
        fullFVCPi = vcat(TMRa[1],-FVC_corr["FVCPi$i"])
        fullFVCK  = vcat(TMRa[1],-FVC_corr["FVCK$i"])

        FVCPiinta = mult .* fullFVCPi .* TMRa
        FVCPiintb = mult .* fullFVCPi .* TMRb
        FVCKinta  = mult .* fullFVCK .* TMRa
        FVCKintb  = mult .* fullFVCK .* TMRb

        push!(FVCPi_a,(alpha/pi)^3 * sum(FVCPiinta) * 1e10)
        push!(FVCPi_b,(alpha/pi)^3 * sum(FVCPiintb) * 1e10)
        push!(FVCK_a,(alpha/pi)^3 * sum(FVCKinta) * 1e10)
        push!(FVCK_b,(alpha/pi)^3 * sum(FVCKintb) * 1e10)
    end

    data_FVC =  Dict{String, Array{uwreal}}()

    data_FVC["FVCPi"] = [FVCPi_a,FVCPi_b]

    data_FVC["FVCK"]  = [FVCK_a,FVCK_b]

    for impr_set in ["1","2"]

        println("   - Starting set "*impr_set)
        println("      - Reading corr")

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
            corr["g_ll"] = corr["g33_ll"] .+ (1/3).*corr["g88_ll_conn"] .+ (4/9).*corr["gcc_ll_conn"]
            corr["g_lc"] = corr["g33_lc"] .+ (1/3).* corr["g88_lc_conn"] .+ (4/9).*corr["gcc_lc_conn"]
        else
            corr["g88_ll"] = corr["g88_ll_conn"] .+ corr["g08_ll_conn"] .+ corr["g88_ll_disc"] .+ corr["g08_ll_disc"]
            corr["g88_lc"] = corr["g88_lc_conn"] .+ corr["g08_lc_conn"] .+ corr["g88_lc_disc"] .+ corr["g08_lc_disc"]

            corr["g_ll"] = corr["g33_ll"] .+ (1/3).*corr["g88_ll"] .+ (4/9).*corr["gcc_ll_conn"] .+ (4/9).*corr["gcc_cc_disc"] .+ (2/(3*sqrt(3))).*corr["gc8_cc_disc"]
            corr["g_lc"] = corr["g33_lc"] .+ (1/3).*corr["g88_lc"] .+ (4/9).*corr["gcc_lc_conn"] .+ (4/9).*corr["gcc_cc_disc"] .+ (2/(3*sqrt(3))).*corr["gc8_cc_disc"] # .+ (4/(3*sqrt(3))).*corr["gc8_cc_disc"]
        end

        sym_points = Int64(length(corr["g33_ll"])/2+1)

        println("      - Computing FVC for 'c'")

        FVC_c_ll = Vector{uwreal}()
        FVC_c_lc = Vector{uwreal}()
        for i in 1:length(FVC_corr)
            fullFVC = vcat(TMRa[1],-FVC_corr["FVC_HP$i"])
            obsll = corr["g_ll"][1:sym_points]
            obslc = corr["g_lc"][1:sym_points]

            FVCintc_ll = (mult^2 .* (fullFVC .* hcat(fullFVC...)) .+ mult .* (obsll .* hcat(fullFVC...) .+ fullFVC .* hcat(obsll...))) .* TMRc
            FVCintc_lc = (mult^2 .* (fullFVC .* hcat(fullFVC...)) .+ mult .* (obslc .* hcat(fullFVC...) .+ fullFVC .* hcat(obslc...))) .* TMRc

            push!(FVC_c_ll,(alpha/pi)^3 * sum(FVCintc_ll) * 1e10)
            push!(FVC_c_lc,(alpha/pi)^3 * sum(FVCintc_lc) * 1e10)
        end

        data_FVC["FVCc_ll_Set"*impr_set] = FVC_c_ll
        data_FVC["FVCc_lc_Set"*impr_set] = FVC_c_lc
    end

    println("   - Writing BDIO")

    # Create path if it does not exist
    pens = joinpath(path_bdio,"HVP",ens.id)
    !ispath(pens) ? mkdir(pens) : nothing

    io = IOBuffer()
    write(io, "$(ens.id)  amu FVC")
    fb = ALPHAdobs_create(joinpath(pens,"$(ens.id)_FVC"), io)

    extra = Dict{String, Any}("Ens" => ens.id)
    ALPHAdobs_write(fb, data_FVC, extra=extra)
        
    ALPHAdobs_close(fb)
end

## <<----------------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> READING TEST <==========================##

ensid = "H102"
impr_set = "1"
extract_data = "HVP"   # "HVP"  "FVC"

extract_data == "HVP" ? fb = BDIO_open(joinpath(path_bdio,"HVP",ensid,ensid*"_HVP_set"*impr_set),"r") : nothing
extract_data == "FVC" ? fb = BDIO_open(joinpath(path_bdio,"HVP",ensid,ensid*"_FVC"),"r") : nothing

res = Dict{String, Any}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    res = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end

res
