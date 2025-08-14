module WriteRead

using HVPobs
using ADerrors
import ADerrors: err

using BDIO
using JLD2
using DelimitedFiles

using ALPHAio

include("../isovModel.jl")

# needed functions

DictComptoKey = Dict{String,Vector{String}}(
    "g33"      => ["g33_ll","g33_lc"],
    "g88"      => ["g88_ll","g88_lc"],
    # only interested in the lc (local-conserved) discr. for the cc conn
    "gCCconn"  => ["gCCconn_SU3_lc"],
    # "gCCconn"  => ["gCCconn_SU3_ll","gCCconn_SU3_lc"],

    "∆ls_amu" => ["∆ls_amu_ll","∆ls_amu_lc"],
    "∆lc_b"   => ["∆lc_b_ll","∆lc_b_lc"],

    # only interested in the cc (conserved-conserved) discr. for the cc disc  & c8 disc
    "gCCdisc" => ["gCCdisc_cc"],
    "gC8disc" => ["gC8disc_cc"],

    "g3333"    => ["g3333_ll","g3333_lc"],
    "g8888"    => ["g8888_ll","g8888_lc"],
    "gCCCC"    => ["gCCCC_ll","gCCCC_lc"],
    "g3388"    => ["g3388_ll","g3388_lc"],
    "g33CC"    => ["g33CC_ll","g33CC_lc"],
    "g88CC"    => ["g88CC_ll","g88CC_lc"]
)

function paste_str(str_vec::Vector{String})::String
    str = str_vec[1]
    for i=collect(2:length(str_vec))
        str *= ","
        str *= string(str_vec[i])
    end
    return str
end

funcOrder_conv = [
    "a2","a2loga","a3","a4",
    "a2phi2","a2phi4","a3phi2","phi2","phi2sqr","phi2log","phi2inv","logphi2","phi4","phi4sqr","phi4log","phi4inv","logphi4",
    "a2y","a2z","a3y","y","ysqr","ylog","yinv","logy","z","zsqr","zlog","zinv","logz"
]
order_map = Dict(s => i for (i, s) in enumerate(funcOrder_conv))
function func_str(func_vec::Vector{Function};Order::Bool=true)::String
    str_vec = string.(func_vec)
    if Order
        str_vec = sort(str_vec, by = x -> order_map[x])
    end
    return paste_str(str_vec)
end

function create_path(path_bdio::String,path_create::Vector{String};OVERWRITE::Bool=false)
    p_ = path_bdio
    if length(path_create) > 1
        for p in path_create[1:end-1]
            p_ = joinpath(p_,p)
            !ispath(p_) ? mkdir(p_) : nothing
        end
    end
    p_ = joinpath(p_,path_create[end])
    if ispath(p_)
        if OVERWRITE
            rm(p_, recursive=true)
            @info("File will be overwriten!")
        else
            error("This information already exist; set 'OVERWRITE = true' to overwrite")
        end
    end
    return p_
end

export paste_str, func_str, create_path


# include("Writer.jl")
# export ...

# General BDIO read structures

function BDIOread_simple(path::String; extra::Bool=false)::Dict
    fb = BDIO_open(path,"r")
    mydict = Dict(); xtra = 0.0
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        xtra = d["extra"]
        mydict = ALPHAdobs_read_next(fb, size=sz, keys=ks)
    end
    extra ? (return mydict, xtra) : (return mydict)
end

function BDIOread_general(path::String; merge::Bool=true, extra::Bool=false)
    fb = BDIO_open(path,"r")
    mydict = Dict(); xtra = 0.0
    mydictvec = Vector{Dict}()
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        # xtra = d["extra"]
        push!(mydictvec,ALPHAdobs_read_next(fb, size=sz, keys=ks))
    end
    BDIO_close!(fb)
    if merge
        if length(mydictvec)>1
            for dict in mydictvec
                merge!(mydict, dict)
            end
        else
            mydict = mydictvec[1]
        end
    else
        mydict = mydictvec
    end
    extra ? (return mydict, xtra) : (return mydict)
end
function BDIOread_general_2first(path::String; merge::Bool=true, extra::Bool=false)
    fb = BDIO_open(path,"r")
    mydict = Dict(); xtra = 0.0
    mydictvec = Vector{Dict}()
    i = 0
    while ALPHAdobs_next_p(fb) && i < 2
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        xtra = d["extra"]
        push!(mydictvec,ALPHAdobs_read_next(fb, size=sz, keys=ks))
        i+=1
    end
    BDIO_close!(fb)
    if merge
        if length(mydictvec)>1
            for dict in mydictvec
                merge!(mydict, dict)
            end
        else
            mydict = mydictvec[1]
        end
    else
        mydict = mydictvec
    end
    extra ? (return mydict, xtra) : (return mydict)
end

function BDIOread_scalar(path::String; extra::Bool=false)
    fb = BDIO_open(path,"r")
    mydict = Dict(); xtra = 0.0
    while ALPHAdobs_next_p(fb)
        if extra
            d = ALPHAdobs_read_parameters(fb)
            xtra = d["extra"]
        end
        mydict = ALPHAdobs_read_next(fb)
    end
    return mydict
    extra ? (return mydict, xtra) : (return mydict)
end

function BDIOread_dim0(path::String; extra::Bool=false)
    fb = BDIO_open(path,"r")
    mydict = Dict(); xtra = 0.0
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        ks = collect(d["keys"])
        xtra = d["extra"]
        mydict = ALPHAdobs_read_next(fb, keys=ks)
    end
    extra ? (return mydict, xtra) : (return mydict)
end

include("Reader.jl")
export BDIOread_t0, BDIOread_TMR, BDIOread_corr, BDIOread_FVCcorr
export BDIOread_HVPens, BDIOread_FVCens
export BDIOread_XYdata, JDL2read_FitRes
export JDL2read_ModelInfo, BDIOread_res, BDIOread_MA, BDIOread_MAtot
export BDIOread_mPP, BDIOread_fPS, BDIOread_mDs_kappaC, BDIOread_mDs
export TXTread_FVCcorr_GS, JDL2read_FVC_ChPT, TXTread_bQ

end