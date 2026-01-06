
## f[4b]

function f4b(w::Union{Int64,Float64,uwreal})

    function f2(w::Union{Int64,Float64,uwreal})
        s = w^2
        Z(s) = - (s-sqrt(s^2+4*s))/(2*s)
        return (1/massmu^2) * s*Z(s)^3 * (1-s*Z(s))/(1+s*Z(s)^2)
    end

    Fl(w,ml) = -5/9 + (4 * ml^2) / (3 * w^2) - (1/3) * sqrt(1 + (4 * ml^2) / w^2) * (1 - (2 * ml^2) / w^2) * log((sqrt(1 + (4 * ml^2) / w^2) - 1) / (sqrt(1 + (4 * ml^2) / w^2) + 1))

    return 2 * f2(w) * Fl(w,masse/massmu)
end


## Tilde f[4b] num.

function Tildef4b_num(t::Union{Int64,Float64,uwreal};rtol::Float64=1e-8)

    integrand(w) = 8π^2 * f4b(w)/w * g(w*t)
    result, error = quadgk(integrand, 0.0, Inf, rtol=rtol)
    return result
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

# we add the numerical kernel for the tau loop

function f4btau(w::Union{Int64,Float64,uwreal})

    function f2(w::Union{Int64,Float64,uwreal})
        s = w^2
        Z(s) = - (s-sqrt(s^2+4*s))/(2*s)
        return (1/massmu^2) * s*Z(s)^3 * (1-s*Z(s))/(1+s*Z(s)^2)
    end

    Fl(w,ml) = -5/9 + (4 * ml^2) / (3 * w^2) - (1/3) * sqrt(1 + (4 * ml^2) / w^2) * (1 - (2 * ml^2) / w^2) * log((sqrt(1 + (4 * ml^2) / w^2) - 1) / (sqrt(1 + (4 * ml^2) / w^2) + 1))

    return 2 * f2(w) * Fl(w,masstau/massmu)
end

function Tildef4btau_num(t::Union{Int64,Float64,uwreal};rtol::Float64=1e-8)

    integrand(w) = 8π^2 * f4btau(w)/w * g(w*t)
    result, error = quadgk(integrand, 0.0, Inf, rtol=rtol)
    return result
end