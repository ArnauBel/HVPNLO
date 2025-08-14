function apply_syst_HVP(hvp::Dict,syst::Dict,diag::String,wind::String,ensid::String;systname::String="")::Dict
    hvpkeys  = keys(hvp)
    systkeys = keys(syst)

    HVP = Dict()
    for key in hvpkeys
        if key in systkeys
            if typeof(hvp[key]) == uwreal
                HVP[key] = hvp[key] + uwreal([0.0,syst[key]],"HVP syst. $systname [$ensid-$diag,$wind,$key]")
            elseif typeof(hvp[key]) == Vector{uwreal}
                HVP[key] = hvp[key] .+ [uwreal([0.0,syst[key][i]],"HVP syst. $systname [$ensid-$diag,$wind,$key]") for i=1:length(syst[key])]
            end
        else
            HVP[key] = hvp[key]
        end
    end
    return HVP
end

function apply_syst_HVP!(HVP::Dict,syst::Dict,diag::String,wind::String,ensid::String;systname::String="")
    systkeys = keys(syst)

    for key in systkeys
        if typeof(HVP[key]) == uwreal
            HVP[key] += uwreal([0.0,syst[key]],"HVP syst. $systname [$ensid-$diag,$wind,$key]")
        elseif typeof(HVP[key]) == Vector{uwreal}
            HVP[key] .+= [uwreal([0.0,syst[key][i]],"HVP syst. $systname [$ensid-$diag,$wind,$key]") for i=1:length(syst[key])]
        end
    end
end

function apply_syst_FVC(fvc::Dict,diag::String,wind::String,ensid::String;factor::Float64=0.15,IMPR_SET::Vector{String}=["1","2"])::Dict
    fvckeys  = keys(fvc)

    FVC = Dict()
    if diag != "NLOc"
        for key in fvckeys
            if typeof(FVC[key]) == uwreal
                FVC[key] = fvc[key] + uwreal([0.0,factor*value(FVC[key])],"FVC syst.  [$ensid-$diag,$wind,$key]")
            elseif typeof(FVC[key]) == Vector{uwreal}
                FVC[key] = fvc[key] .+ [uwreal([0.0,factor*value(FVC[key][i])],"FVC syst. $systname [$ensid-$diag,$wind,$key]") for i=1:length(FVC[key])]
            end
        end
    else
        for impr_set = IMPR_SET
            for key in fvckeys
                FVC[impr_set][key] = fvc[impr_set][key][end] + uwreal([0.0,factor*value(fvc[impr_set][key])],"FVC syst.  [$ensid-$diag,$wind,$key]")
            end
        end
    end
    return FVC
end

function apply_syst_FVC!(FVC::Dict,diag::String,wind::String,ensid::String;factor::Float64=0.15,IMPR_SET::Vector{String}=["1","2"])
    fvckeys  = keys(FVC)

    if diag != "NLOc"
        for key in fvckeys
            if typeof(FVC[key]) == uwreal
                FVC[key] += uwreal([0.0,factor*value(FVC[key])],"FVC syst.  [$ensid-$diag,$wind,$key]")
            elseif typeof(FVC[key]) == Vector{uwreal}
                FVC[key] .+= [uwreal([0.0,factor*value(FVC[key][i])],"FVC syst. [$ensid-$diag,$wind,$key]") for i=1:length(FVC[key])]
            end
        end
    else
        for impr_set = IMPR_SET
            for key in fvckeys
                FVC[impr_set][key] += uwreal([0.0,factor*value(FVC[impr_set][key][i])],"FVC syst.  [$ensid-$diag,$wind,$key]")
            end
        end
    end
    return FVC
end

function HVP_VolCorrect(HVP::Dict,FVC::Dict,diag::String;IMPR_SET::Vector{String}=["1","2"])::Dict
    hvpkeys = keys(HVP[IMPR_SET[1]])
    DISCR = ["ll","lc"]
    if diag != "NLOc"
        COMP = ["g33"]
        "g88_ll" in hvpkeys ? push!(COMP,"g88") : nothing
        "∆ls_amu_ll" in hvpkeys ? push!(COMP,"∆ls_amu") : nothing
        "∆lc_b_ll" in hvpkeys ? push!(COMP,"∆lc_b") : nothing
    else
        COMP = ["3333","8888","3388","33CC","88CC"]
    end

    amu_ens = deepcopy(HVP)

    for impr_set in IMPR_SET
        for comp in COMP
            for discr in DISCR
                if diag != "NLOc"
                    if typeof(amu_ens[impr_set]["$(comp)_$(discr)"]) == uwreal
                        amu_ens[impr_set]["$(comp)_$(discr)"] = HVP[impr_set]["$(comp)_$(discr)"] + FVC["FVC$(comp)"]
                    elseif typeof(amu_ens[impr_set]["$(comp)_$(discr)"]) == Vector{uwreal}
                        amu_ens[impr_set]["$(comp)_$(discr)"] = HVP[impr_set]["$(comp)_$(discr)"] .+ FVC["FVC$(comp)"]
                    else
                        @warn("FVC could not be applied to comp $comp impr. set $impr_set")
                    end
                else
                    amu_ens[impr_set]["g$(comp)_$(discr)"] = HVP[impr_set]["g$(comp)_$(discr)"] + FVC[impr_set]["FVC$(comp)_$(discr)"]
                end
            end
        end
    end
    return amu_ens
end

function HVP_VolCorrect!(HVP::Dict,FVC::Dict,diag::String;IMPR_SET::Vector{String}=["1","2"])
    hvpkeys = keys(HVP[IMPR_SET[1]])
    DISCR = ["ll","lc"]
    if diag != "NLOc"
        COMP = ["g33"]
        "g88_ll" in hvpkeys ? push!(COMP,"g88") : nothing
        "∆ls_amu_ll" in hvpkeys ? push!(COMP,"∆ls_amu") : nothing
        "∆lc_b_ll" in hvpkeys ? push!(COMP,"∆lc_b") : nothing
    else
        COMP = ["3333","8888","3388","33CC","88CC"]
    end

    for impr_set in IMPR_SET
        for comp in COMP
            for discr in DISCR
                if diag != "NLOc"
                    if typeof(HVP[impr_set]["$(comp)_$(discr)"]) == uwreal
                        HVP[impr_set]["$(comp)_$(discr)"] += FVC["FVC$(comp)"]
                    elseif typeof(HVP[impr_set]["$(comp)_$(discr)"]) == Vector{uwreal}
                        HVP[impr_set]["$(comp)_$(discr)"] .+= FVC["FVC$(comp)"]
                    else
                        @warn("FVC could not be applied to comp $comp impr. set $impr_set")
                    end
                else
                    HVP[impr_set]["g$(comp)_$(discr)"] += FVC[impr_set]["FVC$(comp)_$(discr)"]
                end
            end
        end
    end
end

function HVP_3limpr!(HVP::Dict,HVP3l0::Union{Vector{Float64},Float64};IMPR_SET::Vector{String}=["1","2"],meth::String="prod")
    for impr_set in IMPR_SET
        for discr in ["ll","lc"]
            if meth == "prod"
                if typeof(HVP[impr_set]["g33_$discr"]) == Vector{uwreal} && typeof(HVP3l0) == Vector{Float64}
                    HVP[impr_set]["g33_$discr"] .*= HVP3l0 ./ HVP[impr_set]["g33tl_$discr"]
                elseif typeof(HVP[impr_set]["g33_$discr"]) == uwreal && typeof(HVP3l0) == Float64
                    HVP[impr_set]["g33_$discr"]  *= HVP3l0  / HVP[impr_set]["g33tl_$discr"]
                else
                    @warn("3l impr. could not be applied")
                end
            elseif meth == "sum"
                if typeof(HVP[impr_set]["g33_$discr"]) == Vector{uwreal} && typeof(HVP3l0) == Vector{Float64}
                    HVP[impr_set]["g33_$discr"] .+= HVP3l0 .- HVP[impr_set]["g33tl_$discr"]
                elseif typeof(HVP[impr_set]["g33_$discr"]) == uwreal && typeof(HVP3l0) == Float64
                    HVP[impr_set]["g33_$discr"]  += HVP3l0  - HVP[impr_set]["g33tl_$discr"]
                else
                    @warn("3l impr. could not be applied")
                end
            else
                error("Method not recognised, please choose between 'mult' or 'sum'")
            end
        end
    end
end