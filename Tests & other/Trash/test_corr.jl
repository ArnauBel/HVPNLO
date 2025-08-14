# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
# using OrderedCollections

using TimerOutputs

# Include data management and TMR functions

include("../data_management.jl")              #already includes "const.jl"
export get_Z3, get_Z8, get_Z08, corr33, corr88_conn, corr08_conn, corrR, ∆GHP, ensCheck

# Path definition

julia_script_directory = @__DIR__

path_HVP  = joinpath(julia_script_directory, "..", "..", "LatticeData", "HVP_data")
path_rw   = joinpath(julia_script_directory, "..", "..", "LatticeData", "rwf_deflated")
path_ms   = joinpath(julia_script_directory, "..", "..", "LatticeData", "ms_t0_dat")
path_fvc  = joinpath(julia_script_directory, "..", "..", "LatticeData", "JKMPI")


# Ensemble choice

ens = EnsInfo("E250")

##==========================> Connected <==========================##

## G33

g33ll, g33lc = corr33(path_HVP, ens, sector="light", path_rw=path_rw, impr=false, impr_set="1", cons=true, frw_bcwd=true, std=false)

g33ll_std, g33lc_std  = corr33(path_HVP, ens, sector="light", path_rw=path_rw, impr=false, impr_set="1", cons=true, frw_bcwd=true, std=true)

g33ll_impr, g33lc_impr = corr33(path_HVP, ens, sector="light", path_rw=path_rw, impr=false, impr_set="1", cons=true, frw_bcwd=true, std=false)

Z3 = get_Z3(ens, impr_set="1")
g33ll_R = g33ll_impr; g33lc_R = g33lc_impr
renormalize!(g33ll_R, Z3^2)
renormalize!(g33lc_R, Z3)



##==========================> Disconnected <==========================##

#------- G88

# From data management

g88lldisc1, g88lcdisc1 = corrfl_disc(path_HVP, ens, "88", path_rw=path_rw, impr=false, impr_set="1", frw_bcwd=true, discr=["ll","lc"], std=false); uwerr.(g88lldisc1.obs); uwerr.(g88lcdisc1.obs)

g88lldisc_std1, g88lcdisc_std1 = corrfl_disc(path_HVP, ens, "88", path_rw=path_rw, impr=false, impr_set="1", frw_bcwd=true, discr=["ll","lc"], std=true); uwerr.(g88lldisc_std1.obs); uwerr.(g88lcdisc_std1.obs)

g88lldisc_impr1, g88lcdisc_impr1 = corrfl_disc(path_HVP, ens, "88", path_rw=path_rw, impr=true, impr_set="1", frw_bcwd=true, discr=["ll","lc"], std=false); uwerr.(g88lldisc_impr1.obs); uwerr.(g88lcdisc_impr1.obs)

g88lldisc_R1, g88lcdisc_R1 = corrfl_disc(path_HVP, ens, "88", path_rw=path_rw, impr=true, impr_set="1", frw_bcwd=true, discr=["ll","lc"], std=false); uwerr.(g88lldisc_R1.obs); uwerr.(g88lcdisc_R1.obs)

Z8 = get_Z8(ens, impr_set="1")
renormalize!(g88lldisc_R1, Z8^2); uwerr.(g88lldisc_R1.obs)
renormalize!(g88lcdisc_R1, Z8); uwerr.(g88lcdisc_R1.obs)

# From HVPobs

g88lldisc2, g88lcdisc2, g88ccdisc2 = corrDisconnected(path_HVP, ens, "88"; path_rw=path_rw, impr=false, impr_set="1", std=false); uwerr.(g88lldisc2.obs); uwerr.(g88lcdisc2.obs)

g88lldisc_std2, g88lcdisc_std2, g88ccdisc_std2 = corrDisconnected(path_HVP, ens, "88"; path_rw=path_rw, impr=false, impr_set="1", std=true); uwerr.(g88lldisc_std2.obs); uwerr.(g88lcdisc_std2.obs)

g88lldisc_impr2, g88lcdisc_impr2, g88ccdisc_impr2 = corrDisconnected(path_HVP, ens, "88"; path_rw=path_rw, impr=true, impr_set="1", std=false); uwerr.(g88lldisc_impr2.obs); uwerr.(g88lcdisc_impr2.obs)

g88lldisc_R2, g88lcdisc_R2, g88ccdisc_R2 = corrDisconnected(path_HVP, ens, "88"; path_rw=path_rw, impr=true, impr_set="1", std=false); uwerr.(g88lldisc_R2.obs); uwerr.(g88lcdisc_R2.obs)

Z8 = get_Z8(ens, impr_set="1")
renormalize!(g88lldisc_R2, Z8^2); uwerr.(g88lldisc_R2.obs)
renormalize!(g88lcdisc_R2, Z8); uwerr.(g88lcdisc_R2.obs)

##

g88lldisc_impr2.obs[2:10]
g88lcdisc_impr2.obs[2:10]
g88ccdisc_impr2.obs[2:10]

g88lldisc1.obs[2:5]
g88lldisc2.obs[2:5]
g88lldisc_std1.obs[2:5]
g88lldisc_std2.obs[2:5]
g88lldisc_impr1.obs[2:5]
g88lldisc_impr2.obs[2:5]
g88lldisc_R1.obs[2:10]./3
g88lldisc_R2.obs[2:10]./3

#------- G08

# From data management

g08lldisc1, g08lcdisc1 = corrfl_disc(path_HVP, ens, "08", path_rw=path_rw, impr=false, impr_set="1", frw_bcwd=true, discr=["ll","lc"], std=false); uwerr.(g08lldisc1.obs); uwerr.(g08lcdisc1.obs)

g08lldisc_std1, g08lcdisc_std1 = corrfl_disc(path_HVP, ens, "08", path_rw=path_rw, impr=false, impr_set="1", frw_bcwd=true, discr=["ll","lc"], std=true); uwerr.(g08lldisc_std1.obs); uwerr.(g08lcdisc_std1.obs)

g08lldisc_impr1, g08lcdisc_impr1 = corrfl_disc(path_HVP, ens, "08", path_rw=path_rw, impr=true, impr_set="1", frw_bcwd=true, discr=["ll","lc"], std=false); uwerr.(g08lldisc_impr1.obs); uwerr.(g08lcdisc_impr1.obs)

g08lldisc_R1, g08lcdisc_R1 = corrfl_disc(path_HVP, ens, "08", path_rw=path_rw, impr=true, impr_set="1", frw_bcwd=true, discr=["ll","lc"], std=false); uwerr.(g08lldisc_R1.obs); uwerr.(g08lcdisc_R1.obs)

Z8 = get_Z8(ens, impr_set="1")
Z08 = get_Z08(ens, impr_set="1")
renormalize!(g08lldisc_R1, Z8*Z08); uwerr.(g08lldisc_R1.obs)
renormalize!(g08lcdisc_R1, Z08); uwerr.(g08lcdisc_R1.obs)

# From HVPobs

g08lldisc2, g08lcdisc2, g08ccdisc2 = corrDisconnected(path_HVP, ens, "08"; path_rw=path_rw, impr=false, impr_set="1", std=false); uwerr.(g08lldisc2.obs); uwerr.(g08lcdisc2.obs)

g08lldisc_std2, g08lcdisc_std2, g08ccdisc_std2 = corrDisconnected(path_HVP, ens, "08"; path_rw=path_rw, impr=false, impr_set="1", std=true); uwerr.(g08lldisc_std2.obs); uwerr.(g08lcdisc_std2.obs)

g08lldisc_impr2, g08lcdisc_impr2, g08ccdisc_impr2 = corrDisconnected(path_HVP, ens, "08"; path_rw=path_rw, impr=true, impr_set="1", std=false); uwerr.(g08lldisc_impr2.obs); uwerr.(g08lcdisc_impr2.obs)

g08lldisc_R2, g08lcdisc_R2, g08ccdisc_R2 = corrDisconnected(path_HVP, ens, "08"; path_rw=path_rw, impr=true, impr_set="1", std=false); uwerr.(g08lldisc_R2.obs); uwerr.(g08lcdisc_R2.obs)

Z8 = get_Z8(ens, impr_set="1")
Z08 = get_Z08(ens, impr_set="1")
renormalize!(g08lldisc_R2, Z8*Z08); uwerr.(g08lldisc_R2.obs)
renormalize!(g08lcdisc_R2, Z08); uwerr.(g08lcdisc_R2.obs)

##

g08lldisc_impr1.obs[2:10]
g08lldisc_impr2.obs[2:10]

g08lldisc_R2.obs

#------- GC8

# From HVPobs

gC8lldisc2, gC8lcdisc2, gC8ccdisc2 = corrDisconnected(path_HVP, ens, "c8"; path_rw=path_rw, impr=false, impr_set="1", std=false); uwerr.(gC8lldisc2.obs); uwerr.(gC8lcdisc2.obs)

gC8lldisc_std2, gC8lcdisc_std2, gC8ccdisc_std2 = corrDisconnected(path_HVP, ens, "c8"; path_rw=path_rw, impr=false, impr_set="1", std=true); uwerr.(gC8lldisc_std2.obs); uwerr.(gC8lcdisc_std2.obs)

gC8lldisc_impr2, gC8lcdisc_impr2, gC8ccdisc_impr2 = corrDisconnected(path_HVP, ens, "c8"; path_rw=path_rw, impr=true, impr_set="1", std=false); uwerr.(gC8ccdisc_impr2.obs)
g8Clldisc_impr2, g8Clcdisc_impr2, g8Cccdisc_impr2 = corrDisconnected(path_HVP, ens, "8c"; path_rw=path_rw, impr=true, impr_set="1", std=false); uwerr.(g8Cccdisc_impr2.obs)

##

gC8ccdisc_impr2.obs[2:10]
g8Cccdisc_impr2.obs[2:10]

