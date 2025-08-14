
## f[2]

function f2(w::Union{Int64,Float64,uwreal})
    s = w^2
    Z(s) = - (s-sqrt(s^2+4*s))/(2*s)
    return (1/massmu^2) * s*Z(s)^3 * (1-s*Z(s))/(1+s*Z(s)^2)
end


## Tilde f[2] num.

 function Tildef2_num(t::Union{Int64,Float64,uwreal};rtol::Float64=1e-8)

    integrand(w) = 8π^2 * f2(w)/w * g(w*t)
    result, error = quadgk(integrand, 0.0, Inf, rtol=rtol)
    return result
end

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
