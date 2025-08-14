using Revise, HVPobs
using PyPlot, LaTeXStrings
using BDIO, ADerrors, ALPHAio
import ADerrors: err

const path_pi_to_combine = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/low_q_kernel/scale_error_artificial/tmp/"
const path_pi_final = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/PIdata/impr_deriv/low_q_kernel/scale_error_artificial/"

##
fname = filter(x-> occursin(".bdio", x), readdir(path_pi_to_combine))

res_dict = Dict()
for f in fname
    fb = BDIO_open(joinpath(path_pi_to_combine, f), "r")
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        extra = d["extra"]
        ks = collect(d["keys"])
        res = ALPHAdobs_read_next(fb, size=sz, keys=ks)
        res_dict[extra["Ens"]] = res
    end

    BDIO_close!(fb)

end

## pop out new ens
# pop!(res_dict, "E250")
# pop!(res_dict, "D150")
# pop!(res_dict, "D452")
# pop!(res_dict, "D251")
# pop!(res_dict, "D200")

## save file in in single BDIO

ensembles = sort(collect(keys(res_dict)))

io = IOBuffer()
write(io, "PI 08  low q")
fb = ALPHAdobs_create(joinpath(path_pi_final, "PI_08.bdio"), io)
for (k, ens) in enumerate(ensembles)
    extra = Dict{String, Any}("Ens"=> ens)
    data = Dict{String, Array{uwreal}}(
        # "pi33_ll_s1" => res_dict[ens]["pi33_ll_s1"],
        # "pi33_ll_s2" => res_dict[ens]["pi33_ll_s2"],
        "pi08_lc_s1" => res_dict[ens]["pi08_lc_s1"],
        "pi08_lc_s2" => res_dict[ens]["pi08_lc_s2"]

        # "pi33_ll_ILD_s1" => res_dict[ens]["pi33_ll_ILD_s1"],
        # "pi33_lc_ILD_s1" => res_dict[ens]["pi33_lc_ILD_s1"],
        # "pi33_ll_ILD_s2" => res_dict[ens]["pi33_ll_ILD_s2"],
        # "pi33_lc_ILD_s2" => res_dict[ens]["pi33_lc_ILD_s2"]
    )
    ALPHAdobs_write(fb, data, extra=extra)
end
ALPHAdobs_close(fb)

## test reading
fb = BDIO_open(joinpath(path_pi_final, "PI_88_ILD.bdio"), "r")
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