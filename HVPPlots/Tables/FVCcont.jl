# Import packages

using Revise

include("../HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs

using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2
using DelimitedFiles

# include uwreal constants

# include("HVPtool/uwConst.jl")

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "..", "ObsBDIO"),
    # "local" => joinpath(julia_script_directory, "..", "..", "ObsCrosschecks", "Obs4LD(t0=Regensburg)", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "..", "ObsExternal", "PortableSSD", "ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_bPert   = joinpath(julia_script_directory, "..", "..", "PertSD")
path_FVCcont = joinpath(julia_script_directory, "..", "..", "FVCcont")

path_plot    = joinpath(julia_script_directory, "..", "..", "Slides & Plots", "Plots")

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

for diag in ["NLOa","NLOb","NLOa&b"]
    println("\n $diag :")
    for wind in ["SD","ID","LD","NW"]
        println(" - $wind : $(JDL2read_FVC_ChPT(path_FVCcont,diag,wind)*10)")
    end
end