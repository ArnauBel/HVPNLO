
function get_TIC(fitresvec::Vector{FitRes})
    chi2tot = getfield.(fitresvec, :chi2)
    chi2exptot = getfield.(fitresvec, :chi2exp)
    ic = chi2tot .- 2 .* chi2exptot
    return ic
end
get_TIC(fitcat::ModAver.FitCat) = get_TIC(fitcat.fit)

function get_chi2(fitresvec::Vector{FitRes})
    chi2tot = getfield.(fitresvec, :chi2)
    chi2exptot = getfield.(fitresvec, :chi2exp)
    ic = chi2tot./chi2exptot
    return ic
end
get_chi2(fitcat::ModAver.FitCat) = get_chi2(fitcat.fit)

function get_bma_weight(ic::Vector{Float64})
    w = exp.(-0.5 .* ic)
    w = w ./ sum(w)
    return w
end

function get_w_from_fitres(resvec::Vector{FitRes}; norm::Bool=false, AIC::Bool=true)
    ic = AIC ? get_TIC(resvec) : get_chi2(resvec)
    ic = vcat(ic...)
    w = get_bma_weight(ic)
    if norm
        w_min = minimum(w)
        w_max = maximum(w)
        w = (w .- w_min) ./ (w_max - w_min)
    end
    return w
end

function get_w_from_fitcat(catvec::Vector{ModAver.FitCat}; norm::Bool=false)
    ic = get_TIC.(catvec)
    ic = vcat(ic...)
    w = get_bma_weight(ic)
    if norm
        w_min = minimum(w)
        w_max = maximum(w)
        w = (w .- w_min) ./ (w_max - w_min)
    end
    return w
end