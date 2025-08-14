# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
# using OrderedCollections

using TimerOutputs
using Suppressor

# Include TMR functions and extra analysis functioins

include("TMRKernel.jl")
export Tildef4aInner, Tildef4a, Tildef4bInner, Tildef4b, Tildef4cInner, Tildef4c, Tildef2Inner, Tildef2

include("tools/extra_func.jl")              #already includes "const.jl"
export get_Z3, get_Z8, get_Z08, corr33, corr88_conn, corr08_conn, corrR, ∆GHP, ensCheck

# include("amuNLO.jl")
# export amuHVPNLO, amu∆G

# Data treatment 

IMPR      = true
RENORM    = true
STD_DERIV = false

# Path definition

julia_script_directory = @__DIR__

path_HVP   = joinpath(julia_script_directory, "..", "LMEData", "HVP_data")
path_rw    = joinpath(julia_script_directory, "..", "LMEData", "rwf_deflated")
path_ms    = joinpath(julia_script_directory, "..", "LMEData", "ms_t0_data")
path_fvcPI = joinpath(julia_script_directory, "..", "LMEData", "FSE", "JKMPI")
path_fvcK  = joinpath(julia_script_directory, "..", "LMEData", "FSE", "JKMK")

path_coef = joinpath(julia_script_directory, "..", "Coefficients")

if STD_DERIV
    path_bdio = joinpath(julia_script_directory, "..", "ObsBDIOstd")
else
    path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")
end

# Ensamble choice

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302"]
myensList, badensList = ensCheck(ensList, path_HVP, path_rw, path_ms, path_fvcPI, showbad=true)

ensInfo = EnsInfo.(myensList)

isempty(badensList) ? @info("Enough information has been found for all ensembles") : @info("Not enough information has been found concerning ensables $(join(badensList, ", "))\n")


##======================= Production Dict change of names


for ens in ensInfo[end-1:end]

    for impr_set in ["1","2"]

        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_corr_set"*impr_set),"r")

        corr = Dict{String,Vector{uwreal}}()
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            corr = ALPHAdobs_read_next(fb, size=sz, keys=ks)
        end

        println(keys(corr))
    end
end

##

for ens in ensInfo

    println("- Ens $(ens.id)")

    for impr_set in ["1","2"]
        prinln("   - Impr. set $(impr_set)")

        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_corr_set"*impr_set),"r")

        corr = Dict{String,Vector{uwreal}}()
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            corr = ALPHAdobs_read_next(fb, size=sz, keys=ks)
        end

        newcorr = Dict{String,Array{uwreal}}()

        if ens.kappa_l == ens.kappa_s
            mykeys = ["g33_ll", "g88_ll_conn", "g33_lc", "g88_lc_conn"]
        else
            mykeys = ["g33_ll", "gc8_cc_disc", "gcc_cc_disc", "g8c_lc_disc", "g88_ll_disc", "g08_lc_conn", "g08_lc_disc", "gc8_ll_disc", "g8c_ll_disc", "g88_lc_conn", "g88_lc_disc", "gc8_lc_disc", "g08_ll_disc", "g80_lc_disc", "gcc_lc_disc", "g88_ll_conn", "g08_ll_conn", "g8c_cc_disc", "gcc_ll_disc", "g33_lc", "g80_ll_disc"]
        end

        for key in mykeys
            newcorr[key] = corr[key]
        end

        newcorr["gcc_ll_conn_sim"] = corr["gcc_ll_conn"]
        newcorr["gcc_lc_conn_sim"] = corr["gcc_lc_conn"]
        newcorr["gcc_ll_conn_sim+"] = corr["gcc_ll_conn_p"]
        newcorr["gcc_lc_conn_sim+"] = corr["gcc_lc_conn_p"]

        # Create path if it does not exist
        pens = joinpath(path_bdio,"Corr&Kernel&t0",ens.id)
        !ispath(pens) ? mkdir(pens) : nothing

        pBDIO = joinpath(pens,"$(ens.id)_corr_set$(impr_set)")
        if ispath(pBDIO)
            rm(pBDIO, recursive=true)
        else
            error("ens $(ens.id) impr. set $impr_set not found in data storage!!")
        end
        
        io = IOBuffer()
        write(io, "$(ens.id)  HVP correlators, improvement set $(impr_set)")
        fb = ALPHAdobs_create(pBDIO, io)

        extra = Dict{String, Any}("Ens" => ens.id, "Impr_Set" => impr_set)
        ALPHAdobs_write(fb, newcorr, extra=extra)
        
        ALPHAdobs_close(fb)
    end
end