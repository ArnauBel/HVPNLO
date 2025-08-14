##====================================================================================================##
## This code can be used to check wether two BDIO folders contain the same information
##====================================================================================================##

using ADerrors
using HVPobs

using BDIO

include("IO_BDIO.jl")
export read_BDIO

julia_script_directory = @__DIR__

path_bdio    = joinpath(julia_script_directory, "..", "obsBDIO")
path_bdiopre = joinpath(julia_script_directory, "..", "obsBDIOpre")

ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302"]
ensInfo = EnsInfo.(ensList)


#=============> Set of diagrams to compute <=============##

setDiagrams = ["a","b"]

#=============> BDIO; HVP extraction <=============##

# Initialize the Dictionaries (new data)
t0 = Dict{String, uwreal}()
HVP_33ll = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_33ll[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_33lc = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_33lc[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_88ll = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_88ll[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_88lc = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_88lc[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_08ll = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_08ll[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_08lc = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_08lc[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP∆G    = Dict{String, Dict{String, Dict{String, Vector{uwreal}}}}()
[HVP∆G[ensid] = Dict{String, Dict}("a" => Dict{String, Vector{uwreal}}(), "b" => Dict{String, Vector{uwreal}}(), "c" => Dict{String, Vector{uwreal}}()) for ensid in getfield.(ensInfo,:id)]

HVP_Rll  = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_Rll[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_Rlc  = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_Rlc[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]

# Initialize the Dictionaries (old data data)
t0pre = Dict{String, uwreal}()
HVP_33llpre = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_33llpre[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_33lcpre = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_33lcpre[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_88llpre = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_88llpre[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_88lcpre = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_88lcpre[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_08llpre = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_08llpre[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_08lcpre = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_08lcpre[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP∆Gpre    = Dict{String, Dict{String, Dict{String, Vector{uwreal}}}}()
[HVP∆Gpre[ensid] = Dict{String, Dict}("a" => Dict{String, Vector{uwreal}}(), "b" => Dict{String, Vector{uwreal}}(), "c" => Dict{String, Vector{uwreal}}()) for ensid in getfield.(ensInfo,:id)]

HVP_Rllpre  = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_Rllpre[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]
HVP_Rlcpre  = Dict{String, Dict{String, Dict{String, uwreal}}}()
[HVP_Rlcpre[ensid] = Dict{String, Dict}("a" => Dict{String, uwreal}(), "b" => Dict{String, uwreal}(), "c" => Dict{String, uwreal}()) for ensid in getfield.(ensInfo,:id)]

# Read BDIO (new data)
for ensid in getfield.(ensInfo,:id)
    p = joinpath(path_bdio, ensid, string(ensid,"a","_amuObs_Set1.bdio"))
    t0[ensid] = read_BDIO(p, "HVP", "t0")[1]

    for diagram in setDiagrams
        for IMPR_SET in ["1","2"]
            p = joinpath(path_bdio, ensid, string(ensid,diagram,"_amuObs_Set$(IMPR_SET).bdio"))
            HVP_33ll[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "33_ll")[1]; uwerr(HVP_33ll[ensid][diagram][IMPR_SET])
            HVP_33lc[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "33_lc")[1]; uwerr(HVP_33lc[ensid][diagram][IMPR_SET])
            HVP_88ll[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "88_ll")[1]; uwerr(HVP_88ll[ensid][diagram][IMPR_SET])
            HVP_88lc[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "88_lc")[1]; uwerr(HVP_88lc[ensid][diagram][IMPR_SET])
            HVP_08ll[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "08_ll")[1]; uwerr(HVP_08ll[ensid][diagram][IMPR_SET])
            HVP_08lc[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "08_lc")[1]; uwerr(HVP_08lc[ensid][diagram][IMPR_SET])
            HVP∆G[ensid][diagram][IMPR_SET]    = read_BDIO(p, "HVP", "FVC_HP");   uwerr.(HVP∆G[ensid][diagram][IMPR_SET])

            HVP_Rll[ensid][diagram][IMPR_SET]  = HVP_33ll[ensid][diagram][IMPR_SET] + HVP_88ll[ensid][diagram][IMPR_SET] + HVP_08ll[ensid][diagram][IMPR_SET]; uwerr(HVP_Rll[ensid][diagram][IMPR_SET])
            HVP_Rlc[ensid][diagram][IMPR_SET]  = HVP_33lc[ensid][diagram][IMPR_SET] + HVP_88lc[ensid][diagram][IMPR_SET] + HVP_08lc[ensid][diagram][IMPR_SET]; uwerr(HVP_Rlc[ensid][diagram][IMPR_SET])
        end
    end
end

# Read BDIO (old data)
for ensid in getfield.(ensInfo,:id)
    p = joinpath(path_bdiopre, ensid, string(ensid,"a","_amuObs_Set1.bdio"))
    t0pre[ensid] = read_BDIO(p, "HVP", "t0")[1]

    for diagram in setDiagrams
        for IMPR_SET in ["1","2"]
            p = joinpath(path_bdio, ensid, string(ensid,diagram,"_amuObs_Set$(IMPR_SET).bdio"))
            HVP_33llpre[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "33_ll")[1]; uwerr(HVP_33llpre[ensid][diagram][IMPR_SET])
            HVP_33lcpre[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "33_lc")[1]; uwerr(HVP_33lcpre[ensid][diagram][IMPR_SET])
            HVP_88llpre[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "88_ll")[1]; uwerr(HVP_88llpre[ensid][diagram][IMPR_SET])
            HVP_88lcpre[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "88_lc")[1]; uwerr(HVP_88lcpre[ensid][diagram][IMPR_SET])
            HVP_08llpre[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "08_ll")[1]; uwerr(HVP_08llpre[ensid][diagram][IMPR_SET])
            HVP_08lcpre[ensid][diagram][IMPR_SET] = read_BDIO(p, "HVP", "08_lc")[1]; uwerr(HVP_08lcpre[ensid][diagram][IMPR_SET])
            HVP∆Gpre[ensid][diagram][IMPR_SET]    = read_BDIO(p, "HVP", "FVC_HP");   uwerr.(HVP∆Gpre[ensid][diagram][IMPR_SET]   )

            HVP_Rllpre[ensid][diagram][IMPR_SET]  = HVP_33ll[ensid][diagram][IMPR_SET] + HVP_88ll[ensid][diagram][IMPR_SET] + HVP_08ll[ensid][diagram][IMPR_SET]; uwerr(HVP_Rllpre[ensid][diagram][IMPR_SET])
            HVP_Rlcpre[ensid][diagram][IMPR_SET]  = HVP_33lc[ensid][diagram][IMPR_SET] + HVP_88lc[ensid][diagram][IMPR_SET] + HVP_08lc[ensid][diagram][IMPR_SET]; uwerr(HVP_Rlcpre[ensid][diagram][IMPR_SET])
        end
    end
end

##==

dis = 0
for ensid in getfield.(ensInfo,:id)
    @info("=======> ens: $ensid <=======")
    for diagram in setDiagrams
        @info("Diagram $diagram")
        for IMPR_SET in ["1","2"]
            @info("Set $IMPR_SET")
            if (value(HVP_Rll[ensid][diagram][IMPR_SET]) != value(HVP_Rllpre[ensid][diagram][IMPR_SET])) || (ADerrors.err(HVP_Rll[ensid][diagram][IMPR_SET]) != ADerrors.err(HVP_Rllpre[ensid][diagram][IMPR_SET]))
                @info("DISCREPANCY !!!!!!!!!!")
                dis = dis+1
            end
            print("HVP_33ll|| new: $(HVP_33ll[ensid][diagram][IMPR_SET]), old: $(HVP_33llpre[ensid][diagram][IMPR_SET]) \n")
            print("HVP_33lc|| new: $(HVP_33lc[ensid][diagram][IMPR_SET]), old: $(HVP_33lcpre[ensid][diagram][IMPR_SET]) \n")
            print("HVP_88ll|| new: $(HVP_88ll[ensid][diagram][IMPR_SET]), old: $(HVP_88llpre[ensid][diagram][IMPR_SET]) \n")
            print("HVP_88lc|| new: $(HVP_88lc[ensid][diagram][IMPR_SET]), old: $(HVP_88lcpre[ensid][diagram][IMPR_SET]) \n")
            print("HVP_08ll|| new: $(HVP_08ll[ensid][diagram][IMPR_SET]), old: $(HVP_08llpre[ensid][diagram][IMPR_SET]) \n")
            print("HVP_08lc|| new: $(HVP_08lc[ensid][diagram][IMPR_SET]), old: $(HVP_08lcpre[ensid][diagram][IMPR_SET]) \n")
            print("HVP∆G|| new: $(HVP∆G[ensid][diagram][IMPR_SET][end]), old: $(HVP∆Gpre[ensid][diagram][IMPR_SET][end]) \n")

            print("HVP_Rll|| new: $(HVP_Rll[ensid][diagram][IMPR_SET]), old: $(HVP_Rllpre[ensid][diagram][IMPR_SET]) \n")
            print("HVP_Rlc|| new: $(HVP_Rlc[ensid][diagram][IMPR_SET]), old: $(HVP_Rlcpre[ensid][diagram][IMPR_SET]) \n")
        end
    end
end
print("\n The number of discrepancies is: $dis")

##
HVP_Rll["N202"]["a"]["1"]
HVP_Rllpre["N202"]["a"]["1"]
##
(value(HVP_Rll["N202"]["a"]["1"]) != value(HVP_Rllpre["N202"]["a"]["1"])) || (ADerrors.err(HVP_Rll["N202"]["a"]["1"]) != 4)