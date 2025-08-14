module ModAver

using HVPobs
using ADerrors
import ADerrors: err

mutable struct FitCat
    xdata::Array{uwreal}
    ydata::Vector{uwreal}
    info::String
    fit::Vector{FitRes}
    function FitCat(xx, yy, ii)
        a = new()
        a.xdata = xx
        a.ydata = yy
        a.info  = ii
        a.fit   = Vector{FitRes}(undef, 0)
        return a
    end 
end
function Base.show(io::IO, a::FitCat)
    print(io, "FitCat info: ", a.info)
end
export FitCat

include("GetWeight.jl")
export get_w_from_fitres, get_w_from_fitcat

include("ModelAverage.jl")
export model_average, renorm_w!

end