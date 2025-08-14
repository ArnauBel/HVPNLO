
function model_average(results::Vector{uwreal}, ww)
    fin_res = sum(ww .* results)
    aux1 = sum( ww .* results.^2)
    aux2 = sum( ww .* results)^2
    syst = sqrt(abs(value(aux1 - aux2)))
    return [fin_res], syst
end

function renorm_w!(w::Vector{Float64},nFit::Int64)
    length(w) % nFit != 0 ? error("Length of weights and models does not agree") : nothing
    nsteps = Int64(length(w)/nFit)
    for i in collect(0:nsteps-1)
        w[i*nFit+1:(i+1)*nFit] = w[i*nFit+1:(i+1)*nFit]./(nsteps*sum(w[i*nFit+1:(i+1)*nFit]))
    end
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
                # println("wow!")
            end
        end
        return nothing
    end
end

