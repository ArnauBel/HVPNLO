function read_BDIO(path::String, uinfo::Int64)
    r = Vector{uwreal}(undef, 0)
    fb = BDIO_open(path, "r")
    BDIO_seek!(fb)
    
    if BDIO_get_uinfo(fb) == uinfo
        push!(r, read_uwreal(fb))
    end

    while BDIO_seek!(fb, 2) 
        if BDIO_get_uinfo(fb) == uinfo
            push!(r, read_uwreal(fb))
        end
    end

    BDIO_close!(fb)
    return r
end

function read_BDIO(path::String, type::String, obs::String)
    dict_HVP = Dict(
        "t0"     => 0,
        "33_ll"  => 1,
        "33_lc"  => 2,
        "88_ll"  => 3,
        "88_lc"  => 4,
        "08_ll"  => 5,
        "08_lc"  => 6,
        "FVC_HP" => 7
    )
    dict_Integrand = Dict(
        "33_ll"  => 0,
        "33_lc"  => 1,
        "88_ll"  => 2,
        "88_lc"  => 3,
        "08_ll"  => 4,
        "08_lc"  => 5,
        "FVC_HP1" => 6,
        "FVC_HP2" => 7,
        "FVC_HP3" => 8,
        "FVC_HP4" => 9,
        "FVC_HP5" => 10,
        "FVC_HP6" => 11
    )
    dict_Corr = Dict(
        "33_ll"     => 0,
        "33_lc"     => 1,
        "88_ll_con" => 2,
        "88_lc_con" => 3,
        "08_ll_con" => 4,
        "08_lc_con" => 5
    )

    dict2dict = Dict("HVP" => dict_HVP, "Integrand" => dict_Integrand, "Corr" => dict_Corr)
    
    if !(type in keys(dict2dict))
        error("Incorrect type.\ntype = $(keys(dict2dict))")
    end
    if !(obs in keys(dict2dict[type]))
        error("Incorrect obs.\nobs = $(keys(dict2dict[type]))")
    end
    return read_BDIO(path, dict2dict[type][obs])
end
