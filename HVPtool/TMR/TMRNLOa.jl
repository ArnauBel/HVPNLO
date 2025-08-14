
## f[4a]

function f4a(w::Union{Int64,Float64,uwreal})

    function dilog(z)
        -quadgk(t -> log(1 - z*t)/t, 0, 1, rtol=1e-8)[1]
    end

    function F4(u)
        term1 = (78 + 99*u - 57*u^2 - 86*u^3 + 124*u^4 - 37*u^5 + 23*u^6) / (72 * (u - 1)^2 * u * (1 + u))

        term2 = (7 + 8*u - 8*u^3 - 7*u^4) * log(1 - u) / (12 * u^2)

        term3 = (6 + 13*u - 15*u^3 + 4*u^4 + 21*u^5 - 78*u^6 - 11*u^7 + 12*u^8) * log(-u) / (12 * (u - 1)^3 * u * (1 + u)^2)

        term4 = (1 + u) * (6 + 8*u + 7*u^2 - u^3) * log(1 + u) / (12 * u^2)

        log_term = log(-u) * log((1 - u)^2 * (1 + u)) + 2 * dilog(-u) + 4 * dilog(u)
        term5 = (1 / (6 * u^2)) * (-3 - 5*u - 7*u^2 - 5*u^3 - 3*u^4) * log_term

        return term1 + term2 + term3 + term4 + term5
    end

    y(z) = (z - sqrt(z * (z - 4))) / (z + sqrt(z * (z - 4)))

    return (2 * F4(1 / y(-w^2))) / -w^2
end

## Tilde f[4a] num.

 function Tildef4a_num(t::Union{Int64,Float64,uwreal};rtol::Float64=1e-8)

    integrand(w) = 8π^2 * f4a(w)/w * g(w*t)
    result, error = quadgk(integrand, 0.0, 800, rtol=rtol)
    return result
end

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
