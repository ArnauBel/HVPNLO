using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using PyPlot
using Plots
using SpecialFunctions
using LinearAlgebra

## Some constants

include("Const.jl")

include("Reader.jl")

#============= Coefficient extraction functions =============#

## We deffine a new parse function to extract rational numbers from ratioinal-like strings:

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

## Function to extract Mathematica lists from a "contents" string

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

## Custom log regulator

function custom_log(t::Union{Int64,Float64,uwreal})
    if typeof(t) == uwreal
        value(t)  == 0 ? (0) : (log(t))
    else
        t == 0 ? (0) : (log(abs(t)))
    end
end

#============= NLO TRM Kernels =============#

## Tilde f[4a]

function Tildef4aa(t::Union{Int64,Float64,uwreal})
    return t^2 * (197/144 + pi^2/12 - pi^2*log(2)/2 + 3*zeta(3)/4)/2
end

function Tildef4ab13(t::Union{Int64,Float64,uwreal})
    return - pi * t /8 + (log(t) + GammaEuler) * (1 - 5 / (12 * t^2)) + 653/216 - 127 * pi^2/144 - 7 * zeta(3)/4 + 7 * pi^2 * log(2)/6
end

function Tildef4aInner(t::Union{Int64,Float64,uwreal}, an::Vector{Rational{BigInt}}, bn::Vector{Rational{BigInt}}, cn::Vector{Rational{BigInt}}, dn::Vector{Rational{BigInt}}, anb11::Vector{Float64}, anb12::Vector{Float64}, anb21::Vector{Float64}, anb22::Vector{Float64}, anb23::Vector{Float64})
    idcsSmall = collect(4:2:Integer(length(an)*2+2))
    factorials = factorial.(big.(idcsSmall))
    idcsLarge = collect(0:length(anb11)-1)

    tstar = 3.82
    t0 = 5

    if t <= tstar
        return (16 * pi^2 / massmu^2) * sum((an .+  bn .* pi^2 .+ cn .* (custom_log(t)+GammaEuler) .+ dn .* ((custom_log(t)+GammaEuler)^2)) .* ((t .^ idcsSmall) ./ factorials))
    else
        partial_result = Tildef4aa(t) + Tildef4ab13(t)  + sum(anb11 .* ((t0^2/t^2 -  1) .^idcsLarge))/t  + sum(anb12 .* ((t0^2/t^2 -  1) .^idcsLarge))/t^2  + exp(-2*t) * sum(anb21 .* ((t0/t -  1) .^idcsLarge)) + exp(-2*t) * custom_log(t) * sum(anb22 .* ((t0/t -  1) .^idcsLarge))/sqrt(t) + exp(-2*t) * sum(anb23 .* ((t0/t -  1) .^idcsLarge))/sqrt(t)
        return (16 * pi^2 / massmu^2) * partial_result
    end
end

function Tildef4aInner(t::Union{Int64,Float64,uwreal}, an::Union{Vector{Rational{BigInt}},Vector{Float64}}, bn::Union{Vector{Rational{BigInt}},Vector{Float64}}, cn::Union{Vector{Rational{BigInt}},Vector{Float64}}, dn::Union{Vector{Rational{BigInt}},Vector{Float64}}, anb11::Vector{Float64}, anb12::Vector{Float64}, anb21::Vector{Float64}, anb22::Vector{Float64}, anb23::Vector{Float64}, idcsSmall::Vector{Int64}, factorials::Vector{BigInt}, idcsLarge::Vector{Int64})

    tstar = 3.82
    t0 = 5

    if t <= tstar
        return (16 * pi^2 / massmu^2) * sum((an .+  bn .* pi^2 .+ cn .* (custom_log(t)+GammaEuler) .+ dn .* ((custom_log(t)+GammaEuler)^2)) .* ((t .^ idcsSmall) ./ Float64.(factorials,RoundUp)))
    else
        partial_result = Tildef4aa(t) + Tildef4ab13(t)  + sum(anb11 .* ((t0^2/t^2 -  1) .^idcsLarge))/t  + sum(anb12 .* ((t0^2/t^2 -  1) .^idcsLarge))/t^2  + exp(-2*t) * sum(anb21 .* ((t0/t -  1) .^idcsLarge)) + exp(-2*t) * custom_log(t) * sum(anb22 .* ((t0/t -  1) .^idcsLarge))/sqrt(t) + exp(-2*t) * sum(anb23 .* ((t0/t -  1) .^idcsLarge))/sqrt(t)
        return (16 * pi^2 / massmu^2) * partial_result
    end
end

function Tildef4a(t::Union{Vector{Int64},Vector{Float64},Vector{uwreal}}, path_to_coef::String)
    contents4a = read(joinpath(path_to_coef, "NLO_diagram4a.txt"), String)
    an, bn, cn, dn, anb11, anb12, anb21, anb22, anb23 = [Float64.(vec) for vec in extract_coef.(contents4a, ["an","bn","cn","dn","anb11","anb12","anb21","anb22","anb23"])]

    idcsSmall = collect(4:2:Integer(length(an)*2+2))
    factorials = factorial.(big.(idcsSmall))
    idcsLarge = collect(0:length(anb11)-1)
    
    if length(an) != length(bn) != length(cn) != length(dn)
        error("Lengths of SMALL range coeff. not compatible")
    #elseif length(anb11) != length(anb12) != length(anb21) != length(anb22) != length(anb23)       # this is not required for the func. to work
        #error("Lengths of LONG range coeff. not compatible")
    else
        result = Tildef4aInner.(t, Ref(Float64.(an,RoundUp)), Ref(Float64.(bn,RoundUp)), Ref(Float64.(cn,RoundUp)), Ref(Float64.(dn,RoundUp)), Ref(anb11), Ref(anb12), Ref(anb21), Ref(anb22), Ref(anb23), Ref(idcsSmall), Ref(factorials), Ref(idcsLarge))
        #return replace(value.(result), NaN => 0.0)
        return result
    end
end


function Tildef4a(t::Union{Int64,Float64,uwreal}, path_to_coef::String)
    contents4a = read(joinpath(path_to_coef, "NLO_diagram4a.txt"), String)
    an, bn, cn, dn, anb11, anb12, anb21, anb22, anb23 = [Float64.(vec) for vec in extract_coef.(contents4a, ["an","bn","cn","dn","anb11","anb12","anb21","anb22","anb23"])]

    idcsSmall = collect(4:2:Integer(length(an)*2+2))
    factorials = factorial.(big.(idcsSmall))
    idcsLarge = collect(0:length(anb11)-1)
    
    if length(an) != length(bn) != length(cn) != length(dn)
        error("Lengths of SMALL range coeff. not compatible")
    #elseif length(anb11) != length(anb12) != length(anb21) != length(anb22) != length(anb23)       # this is not required for the func. to work
        #error("Lengths of LONG range coeff. not compatible")
    else
        result = Tildef4aInner(t, Float64.(an,RoundUp), Float64.(bn,RoundUp), Float64.(cn,RoundUp), Float64.(dn,RoundUp), anb11, anb12, anb21, anb22, anb23, idcsSmall, factorials, idcsLarge)
        #return replace(value.(result), NaN => 0.0)
        return result
    end
end


## Tilde f[4b]

function Tildef4bInner(t::Union{Int64,Float64,uwreal}, an::Vector{Float64}, bn::Vector{Float64}, cn::Vector{Float64}, dn::Vector{Float64}; OverErr::Bool=false)
    idcs = collect(4:2:Integer(length(an)*2+2))
    factorials = factorial.(big.(idcs))

    tstar = 4

    if t <= tstar || OverErr
        return (16 * pi^2 / massmu^2) * sum((an .+  bn .* pi^2 .+ cn .* (custom_log(t)+GammaEuler) .+ dn .* ((custom_log(t)+GammaEuler)^2)) .* ((t .^ idcs) ./ Float64.(factorials,RoundUp)))
    else
        error("mµ t > 4: No numerical precision <1e-8 is guaranteed beyond this point, set OverErr = true to overrun this error.")
    end
end

function Tildef4bInner(t::Union{Int64,Float64,uwreal}, path_to_coef::String; OverErr::Bool=false)
    contents4b = read(joinpath(path_to_coef, "NLO_diagram4b.txt"), String)
    an, bn, cn, dn= extract_coef.(contents4b, ["an","bn","cn","dn"])

    idcs = collect(4:2:Integer(length(an)*2+2))
    factorials = factorial.(big.(idcs))

    tstar = 4

    if t <= tstar || OverErr
        return (16 * pi^2 / massmu^2) * sum((an .+  bn .* pi^2 .+ cn .* (custom_log(t)+GammaEuler) .+ dn .* ((custom_log(t)+GammaEuler)^2)) .* ((t .^ idcs) ./ Float64.(factorials,RoundUp)))
    else
        error("mµ t > 4: No numerical precision <1e-8 is guaranteed beyond this point, set OverErr = true to overrun this error.")
    end
end

function Tildef4bInner(t::Union{Int64,Float64,uwreal}, an::Vector{Float64}, bn::Vector{Float64}, cn::Vector{Float64}, dn::Vector{Float64}, idcs::Vector{Int64}, factorials::Vector{BigInt}; OverErr::Bool=false)
    first_term = an .+  bn .* pi^2 .+ cn .* (custom_log(t)+GammaEuler) .+ dn .* ((custom_log(t)+GammaEuler)^2)
    second_term = (t .^ idcs) ./ Float64.(factorials,RoundUp)

    tstar = 4

    if t <= tstar || OverErr
        return sum(first_term .* second_term)
    else
        error("mµ t > 4: No numerical precision <1e-8 is guaranteed beyond this point, set OverErr = true to overrun this error.")
    end
end

function Tildef4b(t::Union{Vector{Int64},Vector{Float64},Vector{uwreal}}, path_to_coef::String; OverErr::Bool=false)
    contents4b = read(joinpath(path_to_coef, "NLO_diagram4b.txt"), String)
    an, bn, cn, dn = extract_coef.(contents4b, ["an","bn","cn","dn"])

    idcs = collect(4:2:Integer(length(an)*2+2))
    factorials = factorial.(big.(idcs))
    
    length(an) != length(bn) != length(cn) != length(dn) ? error("Lengths of coeff. not compatible") : nothing
    
    result = (16 * pi^2 / massmu^2) .* Tildef4bInner.(t, Ref(an), Ref(bn), Ref(cn), Ref(dn), Ref(idcs), Ref(factorials), OverErr=OverErr)
    #return replace(value.(result), NaN => 0.0)
    return result
end

function Tildef4b(t::Union{Int64,Float64,uwreal}, path_to_coef::String; OverErr::Bool=false)
    contents4b = read(joinpath(path_to_coef, "NLO_diagram4b.txt"), String)
    an, bn, cn, dn = extract_coef.(contents4b, ["an","bn","cn","dn"])

    idcs = collect(4:2:Integer(length(an)*2+2))
    factorials = factorial.(big.(idcs))
    
    length(an) != length(bn) != length(cn) != length(dn) ? error("Lengths of coeff. not compatible") : nothing
    
    result = (16 * pi^2 / massmu^2) * Tildef4bInner(t, an, bn, cn, dn, idcs, factorials, OverErr=OverErr)
    #return replace(value.(result), NaN => 0.0)
    return result
end

## Tilde f[4c]

function Tildef4cInner(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal}, a2n::Vector{Rational{BigInt}}, b2n::Vector{Rational{BigInt}}, a3n::Vector{Rational{BigInt}}, b3n::Vector{Rational{BigInt}})

    aidcs = collect(4:2:Integer(length(a2n)*2-2))
    bidcs = collect(6:2:Integer(length(a3n)*2-2))
        
    aterm1 = sum(a2n[3:end] .* (t^2 .* (tau .^aidcs) .+ tau^2 .* (t .^aidcs)))
    aterm2 = sum(b2n[3:end] .* (t^2 .* (custom_log(tau)+GammaEuler) .* (tau .^aidcs) .+ tau^2 .* (custom_log(t)+GammaEuler) .* (t .^aidcs)))

    bterm1 = 2 * sum(a3n[4:end] .* ((t + tau) .^bidcs .+ (t - tau) .^bidcs .- 2 .* t .^bidcs .- 2 .* tau .^bidcs)) 
    bterm2 = 2 * sum(b3n[4:end] .* ((custom_log(t+tau)+GammaEuler) .* (t + tau) .^bidcs .+ (custom_log(t-tau)+GammaEuler) .* (t - tau) .^bidcs .- 2 .* (custom_log(t)+GammaEuler) .* t .^bidcs .- 2 .* (custom_log(tau)+GammaEuler) .* tau .^bidcs))

    return (32 * pi^4 / massmu^4) * (aterm1 + aterm2 + bterm1 + bterm2)  
end

function Tildef4cInner(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal}, path_to_coef::String)

    contents4c = read(joinpath(path_to_coef, "NLO_diagram4c.txt"), String)
    a2n, b2n, a3n, b3n = [Float64.(vec) for vec in extract_coef.(contents4c, ["a2n","b2n","a3n","b3n"])]
    
    aidcs = collect(4:2:Integer(length(a2n)*2-2))
    bidcs = collect(6:2:Integer(length(a3n)*2-2))
        
    aterm1 = sum(a2n[3:end] .* (t^2 .* (tau .^aidcs) .+ tau^2 .* (t .^aidcs)))
    aterm2 = sum(b2n[3:end] .* (t^2 .* (custom_log(tau)+GammaEuler) .* (tau .^aidcs) .+ tau^2 .* (custom_log(t)+GammaEuler) .* (t .^aidcs)))
    
    bterm1 = 2 * sum(a3n[4:end] .* ((t + tau) .^bidcs .+ (t - tau) .^bidcs .- 2 .* t .^bidcs .- 2 .* tau .^bidcs)) 
    bterm2 = 2 * sum(b3n[4:end] .* ((custom_log(t+tau)+GammaEuler) .* (t + tau) .^bidcs .+ (custom_log(t-tau)+GammaEuler) .* (t - tau) .^bidcs .- 2 .* (custom_log(t)+GammaEuler) .* t .^bidcs .- 2 .* (custom_log(tau)+GammaEuler) .* tau .^bidcs))
    
    return (32 * pi^4 / massmu^4) * (aterm1 + aterm2 + bterm1 + bterm2)
end

function Tildef4cInner(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal}, a2n::Union{Vector{Rational{BigInt}},Vector{Float64}}, b2n::Union{Vector{Rational{BigInt}},Vector{Float64}}, a3n::Union{Vector{Rational{BigInt}},Vector{Float64}}, b3n::Union{Vector{Rational{BigInt}},Vector{Float64}}, aidcs::Vector{Int64}, bidcs::Vector{Int64})
        
    aterm1 = sum(a2n[3:end] .* (t^2 .* (tau .^aidcs) .+ tau^2 .* (t .^aidcs)))
    aterm2 = sum(b2n[3:end] .* (t^2 .* (custom_log(tau)+GammaEuler) .* (tau .^aidcs) .+ tau^2 .* (custom_log(t)+GammaEuler) .* (t .^aidcs)))
    
    bterm1 = 2 * sum(a3n[4:end] .* ((t + tau) .^bidcs .+ (t - tau) .^bidcs .- 2 .* t .^bidcs .- 2 .* tau .^bidcs)) 
    bterm2 = 2 * sum(b3n[4:end] .* ((custom_log(t+tau)+GammaEuler) .* (t + tau) .^bidcs .+ (custom_log(t-tau)+GammaEuler) .* (t - tau) .^bidcs .- 2 .* (custom_log(t)+GammaEuler) .* t .^bidcs .- 2 .* (custom_log(tau)+GammaEuler) .* tau .^bidcs))
    
    #return (32 * pi^4 / massmu^4) * (aterm1 + aterm2 + bterm1 + bterm2)
    return aterm1 + aterm2 + bterm1 + bterm2
end

function Tildef4cSUPa(t::Union{Int64,Float64,uwreal}, cn::Vector{Float64}, cidcs::Vector{Int64})
    t0 = 2.2
    mySUP = exp(-2*t)/sqrt(t) * sum(cn .* (t0/t -  1) .^cidcs)
    #return (32 * pi^4 / massmu^4)*sqrt(pi) * mySUP
    return sqrt(pi) * mySUP
end

function Tildef4cDOM_blue(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal})
    myDOM = -11/6 - t^2/2 + t^2/tau^2 - tau^2/2 + tau^2/t^2 + (t^2*tau^2)/4 - 1/48*pi*(t + tau)*(-45 + 4*(t+tau)^2) + 1/48*pi*(t-tau)*(-45 + 4*(t-tau)^2) + 2*(t^2 + tau^2 - 1)*(GammaEuler + custom_log(t*tau)) + (1 - (t+tau)^2)*custom_log(t+tau) + (1 - (t-tau)^2)*custom_log(t-tau)
    #return (32 * pi^4 / massmu^4) * myDOM
    return myDOM
end
function Tildef4cSUPP_blue(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal}, c2n::Vector{Float64}, c3n::Vector{Float64}, cidcs::Vector{Int64})
    mySUPP = 0.25*(tau^2*Tildef4cSUPa(t,c2n,cidcs) + t^2*Tildef4cSUPa(tau,c2n,cidcs)) + 2*(Tildef4cSUPa(t+tau,c3n,cidcs) + Tildef4cSUPa(t-tau,c3n,cidcs) - 2*Tildef4cSUPa(t,c3n,cidcs) - 2*Tildef4cSUPa(tau,c3n,cidcs))
    #return (32 * pi^4 / massmu^4) * mySUPP
    return mySUPP
end

function Tildef4cDOM_red(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal})
    myDOM = -11/4 - t^2/2 + t^2/tau^2 - tau^2/2 + tau^2/t^2 + (t^2*tau^2)/4 - 1/48*pi*(t + tau)*(-45 + 4*(t+tau)^2) + 2*(t^2 + tau^2 - 1)*(GammaEuler + custom_log(t*tau)) + (1 - (t+tau)^2)*custom_log(t+tau) - GammaEuler
    #return (32 * pi^4 / massmu^4) * myDOM
    return myDOM
end
function Tildef4cSUPP_red(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal}, c2n::Vector{Float64}, c3n::Vector{Float64}, cidcs::Vector{Int64})
    mySUPP = 0.25*(tau^2*Tildef4cSUPa(t,c2n,cidcs) + t^2*Tildef4cSUPa(tau,c2n,cidcs)) + 2*(Tildef4cSUPa(t+tau,c3n,cidcs) - 2*Tildef4cSUPa(t,c3n,cidcs) - 2*Tildef4cSUPa(tau,c3n,cidcs))
    #return (32 * pi^4 / massmu^4) * mySUPP
    return mySUPP
end

function Tildef4cDOM_yellow(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal})
    myDOM =  tau^2/t^2+ 2*(t^2 + tau^2 - 1)custom_log(t) + (1 - (t+tau)^2)*custom_log(t+tau) + (1 - (t-tau)^2)*custom_log(t-tau)    #return (32 * pi^4 / massmu^4) * myDOM
    return myDOM
end
function Tildef4cSUPP_yellow(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal}, c2n::Vector{Float64}, c3n::Vector{Float64}, cidcs::Vector{Int64})
    mySUPP = 0.25*tau^2*Tildef4cSUPa(t,c2n,cidcs) + 2*(Tildef4cSUPa(t+tau,c3n,cidcs) + Tildef4cSUPa(t-tau,c3n,cidcs) - 2*Tildef4cSUPa(t,c3n,cidcs))
    #return (32 * pi^4 / massmu^4) * mySUPP
    return mySUPP
end

function Tildef4c(t::Union{Vector{Int64},Vector{Float64},Vector{uwreal}}, path_to_coef::String)
    
    contents4c = read(joinpath(path_to_coef, "NLO_diagram4c.txt"), String)
    a2n, b2n, a3n, b3n, c2n, c3n = [Float64.(vec) for vec in extract_coef.(contents4c, ["a2n","b2n","a3n","b3n","c2n","c3n"])]
    
    length(a2n) != length(b2n) ? error("The lengths of a2n and b2n should be the same") : nothing
    length(a3n) != length(b3n) ? error("The lengths of a3n and b3n should be the same") : nothing
    length(c2n) != length(c3n) ? error("The lengths of c2n and c3n should be the same") : nothing

    aidcs  = collect(4:2:Int64(length(a2n)*2-2))
    bidcs  = collect(6:2:Int64(length(a3n)*2-2))
    cidcs  = collect(0:1:Int64(length(c2n)-1))

    n = length(t)
    result_matrix = Array{Any}(undef, n, n)

    for i in 1:n
        for j in 1:i
            ti = t[i]; tj = t[j]
            if ti + tj <= 3.8
                result_matrix[i, j] = Tildef4cInner(ti,tj,a2n,b2n,a3n,b3n,aidcs,bidcs)
            elseif ti-tj <= 3 && tj >= 1.3
                result_matrix[i, j] = Tildef4cDOM_red(ti, tj) + Tildef4cSUPP_red(ti, tj, c2n, c3n, cidcs) + (GammaEuler-7/4)*(ti-tj)^2 + 2*sum((a3n[3:end].+b3n[3:end].*(GammaEuler+custom_log(ti-tj))).*(ti-tj).^aidcs)
            elseif tj <= 2.4
                result_matrix[i, j] = Tildef4cDOM_yellow(ti, tj) + Tildef4cSUPP_yellow(ti, tj, c2n, c3n, cidcs) + 3*tj^2 + ti^2*sum((a2n[3:end].+b2n[3:end].*(GammaEuler+custom_log(tj))).* tj.^aidcs) - 4*sum((a3n[3:end].+b3n[3:end].*(GammaEuler+custom_log(tj))).* tj.^aidcs)
            else
                result_matrix[i, j] = Tildef4cDOM_blue(ti, tj) + Tildef4cSUPP_blue(ti, tj, c2n, c3n, cidcs)
            end
            result_matrix[j,i] = result_matrix[i,j]
        end
    end
    return (32 * pi^4 / massmu^4) .* result_matrix
    #return result_matrix
end


#============= LO TRM Kernels =============#

## Tilde f[2]

function Tildef2DOM(t::Union{Int64,Float64,uwreal})
    return 1/(2*t^2) - (pi*t)/4 + t^2/8 + 1/4*(-1 + 4*GammaEuler + 4*custom_log(t))
end


function Tildef2Inner(t::Union{Int64,Float64,uwreal}, an::Vector{Rational{BigInt}}, cn::Vector{Rational{BigInt}}, anb23::Vector{Float64})
    idcsSmall = collect(4:2:Integer(length(an)*2+2))
    factorials = factorial.(big.(idcsSmall))
    idcsLarge = collect(0:length(anb23)-1)

    tstar = 1.05
    v0 = 0.7

    if t <= tstar
        return (16 * pi^2 / massmu^2) * sum((an .+ cn .* (custom_log(t)+GammaEuler)) .* ((t .^ idcsSmall) ./ factorials))
    else
        partial_result = Tildef2DOM(t) - (sqrt(pi)/(8*sqrt(t)))*exp(-2*t) * sum(anb23 .* ((1/t -  v0) .^idcsLarge))
        return (16 * pi^2 / massmu^2) * partial_result
    end
end

function Tildef2Inner(t::Union{Int64,Float64,uwreal}, an::Union{Vector{Rational{BigInt}},Vector{Float64}}, cn::Union{Vector{Rational{BigInt}},Vector{Float64}}, anb23::Vector{Float64}, idcsSmall::Vector{Int64}, factorials::Vector{BigInt}, idcsLarge::Vector{Int64})

    tstar = 1.05
    v0 = 0.7

    if t <= tstar
        return (16 * pi^2 / massmu^2) * sum((an .+ cn .* (custom_log(t)+GammaEuler)) .* ((t .^ idcsSmall) ./ Float64.(factorials,RoundUp)))
    else
        partial_result = Tildef2DOM(t) - (sqrt(pi)/(8*sqrt(t))) * exp(-2*t) * sum(anb23 .* ((1/t -  v0) .^idcsLarge))
        return (16 * pi^2 / massmu^2) * partial_result
    end
end

function Tildef2(t::Union{Vector{Int64},Vector{Float64},Vector{uwreal}}, path_to_coef::String)
    contents4a = read(joinpath(path_to_coef, "LO.txt"), String)
    an, cn, anb23 = [Float64.(vec) for vec in extract_coef.(contents4a, ["an","cn","anb23"])]

    idcsSmall = collect(4:2:Integer(length(an)*2+2))
    factorials = factorial.(big.(idcsSmall))
    idcsLarge = collect(0:length(anb23)-1)
    
    if length(an) != length(cn)
        error("Lengths of SMALL range coeff. not compatible")
    else
        result = Tildef2Inner.(t, Ref(Float64.(an,RoundUp)), Ref(Float64.(cn,RoundUp)), Ref(anb23), Ref(idcsSmall), Ref(factorials), Ref(idcsLarge))
        # result = Tildef2Inner.(t, Ref(an), Ref(cn), Ref(anb23), Ref(idcsSmall), Ref(factorials), Ref(idcsLarge))
        #return replace(value.(result), NaN => 0.0)
        return result
    end
end

function Tildef2(t::Union{Int64,Float64,uwreal}, path_to_coef::String)
    contents4a = read(joinpath(path_to_coef, "LO.txt"), String)
    an, cn, anb23 = [Float64.(vec) for vec in extract_coef.(contents4a, ["an","cn","anb23"])]

    idcsSmall = collect(4:2:Integer(length(an)*2+2))
    factorials = factorial.(big.(idcsSmall))
    idcsLarge = collect(0:length(anb23)-1)
    
    if length(an) != length(cn)
        error("Lengths of SMALL range coeff. not compatible")
    else
        result = Tildef2Inner(t, Float64.(an,RoundUp), Float64.(cn,RoundUp), anb23, idcsSmall, factorials, idcsLarge)
        # result = Tildef2Inner.(t, Ref(an), Ref(cn), Ref(anb23), Ref(idcsSmall), Ref(factorials), Ref(idcsLarge))
        #return replace(value.(result), NaN => 0.0)
        return result
    end
end
