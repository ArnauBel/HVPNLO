module TMR

using ADerrors
import ADerrors: err

using QuadGK

using SpecialFunctions

include("../Const/Const.jl")

include("../WriteRead/WriteRead.jl")
using .WriteRead

# We deffine a new parse function to extract rational numbers from ratioinal-like strings:

function Baseparse(x::Union{String,SubString})#::Rational{BigInt}
    try
        if occursin('/', x)
            ms, ns = split(x, '/')
            m = parse(BigInt, ms)
            n = parse(BigInt, ns)
            return m // n
        else
            return parse(BigInt, x) // 1
        end
    catch
        return parse(Float64, x)
    end
end

# Function to extract Mathematica lists from a "contents" string

function extract_coef(contents::String, list_name::String)
    start_index = findfirst(x -> contains(x, "# $(list_name)"), split(contents, '\n'))
    if start_index === nothing
        error("List of coefficients $(an) has not been found")
    end

    start_index += 1
    end_index = findfirst(x -> isempty(x) || x[1] == '#', split(contents, '\n')[start_index:end])
    end_index = end_index === nothing ? length(contents) : start_index + end_index - 2

    raw_list = string.(split(contents, '\n')[start_index:end_index])
    
    return Baseparse.(string.(raw_list))
end

# Custom log regulator

function custom_log(t::Union{Int64,Float64,uwreal})
    if typeof(t) == uwreal
        value(t)  == 0 ? (0) : (log(t))
    else
        t == 0 ? (0) : (log(abs(t)))
    end
end

# g(w) function for numerical transformation

function g(x)
    return x^2 - 4*sin(x/2)^2
end


include("TMRLO.jl")
export f2
export Tildef2_num
export Tildef2, Tildef2Inner

include("TMRNLOa.jl")
export f4a
export Tildef4a_num
export Tildef4a, Tildef4aInner

include("TMRNLOb.jl")
export f4b
export f4btau
export Tildef4b_num
export Tildef4btau_num
export Tildef4b, Tildef4bInner

include("TMRNLOc.jl")
export Tildef4c, Tildef4cInner

end