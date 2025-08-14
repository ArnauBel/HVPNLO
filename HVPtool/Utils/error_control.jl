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

function add_t0_err!(obs::Vector{uwreal}, t0phys::uwreal)
    uwerr(t0phys); uwerr.(obs)
    for k in eachindex(obs)
        err_t0 = abs.(mchist(obs[k], "sqrtt0 [fm]") * err(t0phys) / artificial_err)
        obs[k] = obs[k] + uwreal([0.0, err_t0[1]], "sqrtt0 [fm]")
    end
    return nothing
end

function get_t0err(obs::Vector{uwreal}, SCALEphys::uwreal; resc::Bool=false)
    uwerr(SCALEphys); uwerr.(obs)
    SCALEstr = !resc ? "sqrtt0 [fm]" : "fPi [GeV]"
    err_t0 = []
    for k in eachindex(obs)
        push!(err_t0, abs.(mchist(obs[k], SCALEstr)[1] * err(SCALEphys) / artificial_err))
    end
    return err_t0
end