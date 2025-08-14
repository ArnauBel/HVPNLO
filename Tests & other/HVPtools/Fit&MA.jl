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

function get_TIC(fitresvec::Vector{FitRes})
    chi2tot = getfield.(fitresvec, :chi2)
    chi2exptot = getfield.(fitresvec, :chi2exp)
    ic = chi2tot .- 2 .* chi2exptot
    return ic
end
get_TIC(fitcat::FitCat) = get_TIC(fitcat.fit)

function get_chi2(fitresvec::Vector{FitRes})
    chi2tot = getfield.(fitresvec, :chi2)
    chi2exptot = getfield.(fitresvec, :chi2exp)
    ic = chi2tot./chi2exptot
    return ic
end
get_chi2(fitcat::FitCat) = get_chi2(fitcat.fit)

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

function get_w_from_fitcat(catvec::Vector{FitCat}; norm::Bool=false)
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

function renorm_w!(w::Vector{Float64},nFit::Int64)
    length(w) % nFit != 0 ? error("Length of weights and models does not agree") : nothing
    nsteps = Int64(length(w)/nFit)
    for i in collect(0:nsteps-1)
        w[i*nFit+1:(i+1)*nFit] = w[i*nFit+1:(i+1)*nFit]./(nsteps*sum(w[i*nFit+1:(i+1)*nFit]))
    end
end

function model_average(results::Vector{uwreal}, ww)
    fin_res = sum(ww .* results)
    aux1 = sum( ww .* results.^2)
    aux2 = sum( ww .* results)^2
    syst = sqrt(abs(value(aux1 - aux2)))
    return [fin_res], syst
end


function set_fluc_to_zero!(a::uwreal, id_str::String)
    ws = ADerrors.wsg
    id_int = ws.str2id[id_str]
    idx = ADerrors.find_mcid(a, id_int)
    if isnothing(idx)
        error("No error available... maybe run uwerr")
    else
        nd = ws.fluc[ws.map_ids[a.ids[idx]]].nd
        for j in 1:length(a.prop)
            if (a.prop[j] && ((ws.map_nob[j] == a.ids[idx])))
                ws.fluc[j].delta[:] .= 0.0
            end
        end
        return nothing
    end
end

##


function set_err!(a::uwreal, id_str::String, arterr::Float64)
    ws = ADerrors.wsg
    id_int = ws.str2id[id_str]
    idx = ADerrors.find_mcid(a, id_int)
    if isnothing(idx)
        error("No error available... maybe run uwerr")
    else
        nd = ws.fluc[ws.map_ids[a.ids[idx]]].nd
        for j in 1:length(a.prop)
            if (a.prop[j] && ((ws.map_nob[j] == a.ids[idx])))
                ws.fluc[j].delta[:] .= arterr
                println("wow!")
            end
        end
        return nothing
    end
end

