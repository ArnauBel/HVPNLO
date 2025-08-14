# this code takes a bdio file with obs saved for each ensembles and replaces one 
# specific ensemble for which statistics has been updated

using Revise
using HVPobs
using ADerrors
using BDIO
using ALPHAio
import ADerrors:err

include("const.jl")
include("types.jl")
# include("plot_utils.jl")
include("IO_BDIO.jl")
# include("tools.jl")


const path_main ="/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/low_q_kernel/scale_error_artificial" 
path_new_ens = joinpath(path_main,"tmp/PI_33_spectroscopy.bdio")
path_old_bdio = joinpath(path_main, "old_data/update_data_no_spectroscopy/PI_33.bdio")
path_save_updated_bdio = path_main


# ensembles = ["E300", "D251", "F300", "J306", "J307", "D450"] # ensembles to replace
ensembles = [ "E250", "D200"] # ensembles to replace

## read old BDIO

fb = BDIO_open(path_old_bdio, "r")
res_old = Dict()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    if extra["Ens"] ∈ ensembles
        continue
    end
    res_old[extra["Ens"]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

# change_id(from="D450", to="D450_old")
## read updated statistics ensembles from bdio
fb = BDIO_open(path_new_ens, "r")
res_updated = Dict()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    res_updated[extra["Ens"]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)


## merge the new ensembles with the older one dictionaries 
res_new = merge(res_old, res_updated)

## save the newly created BDIO with updated statistics
@info("Saving PI CC results in BDIO")
io = IOBuffer()
write(io, "PI 33 low q. ")
fname= basename(path_old_bdio)
fb = ALPHAdobs_create(joinpath(path_save_updated_bdio, fname), io)

ens_tot = collect(keys(res_new))
for (k,ens) in enumerate(ens_tot)
    extra = Dict{String,Any}("Ens"=> ens)
    data = res_new[ens]
    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)
println("# Saving complete!")


## test reading new saved file
fb = BDIO_open(joinpath(path_save_updated_bdio, fname), "r")
res = Dict()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    sz = tuple(d["size"]...)
    extra = d["extra"]
    ks = collect(d["keys"])
    res[extra["Ens"]] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

