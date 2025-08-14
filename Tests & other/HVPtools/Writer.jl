using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2


function paste_str(str_vec::Vector{String})::String
    str = str_vec[1]
    for i=collect(2:length(str_vec))
        str *= ","
        str *= string(str_vec[i])
    end
    return str
end

funcOrder_conv = ["a2","a2loga","a3","a4","a2phi2","a2phi4","a3phi2","phi2","phi2sqr","phi2log","phi2inv","logphi2","phi4","phi4sqr","phi4log","phi4inv","logphi4"]
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

