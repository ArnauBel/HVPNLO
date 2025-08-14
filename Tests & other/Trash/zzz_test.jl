
using Revise


include("HVPtool/HVPtool.jl"); using .HVPtool
include("HVPtool/uwConst.jl")

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "ClusterMount", "mogon_mount", "ObsBDIO")
)


BDIOread_TMR(path_bdio_dict["local"],"A654","all")


path_HVP   = joinpath(julia_script_directory, "..", "ClusterMount", "mogon_mount", "HVPdata", "HVP_data")
path_rw_   = joinpath(julia_script_directory, "..", "ClusterMount", "mogon_mount", "HVPdata", "reweight")
path_ms    = joinpath(julia_script_directory, "..", "ClusterMount", "mogon_mount", "HVPdata", "ms_t0_dat")
path_fvcPI = joinpath(julia_script_directory, "..", "ClusterMount", "mogon_mount", "HVPdata", "FSE", "JKMPI")
path_fvcK  = joinpath(julia_script_directory, "..", "ClusterMount", "mogon_mount", "HVPdata", "FSE", "JKMK")

using HVPobs

ens = EnsInfo("E250")

path_rw = !isempty(filter(x->occursin(ens.id, x), readdir(joinpath(path_rw_,"reweight_deflated"), join=true))) ? joinpath(path_rw_,"reweight_deflated") : path_rw_

g33_ll, g33_lc = corr33(path_HVP, ens, path_rw=path_rw, impr=true, impr_set="1", cons=true, frw_bcwd=true, std=false)


get_w_from_fitres(fitRes["gCCconn_SU3_ll"])

path_bdio = path_bdio_dict["clust"]

using ADerrors

t0 = uwreal[]

t0 = BDIOread_t0(path_bdio,ens)


factor = hbarc * sqrt(t0)/sqrtt0_ph

sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1)
t_hat = (massmu/factor) .* collect(0:sym_points-1)

path_coef  = joinpath(julia_script_directory, "..", "KernelCoeff")

TMR = factor^2 .* Tildef2(t_hat,path_coef)
