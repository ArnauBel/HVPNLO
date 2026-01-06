using HVPobs


##-- corr functions

@doc raw"""
    corr33(path_data::String, ens::EnsInfo; path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, cons::Bool=true, std::Bool=false, lma::Bool=true)

This function return the G33 correlator given path_data and the EnsInfo data type of the corresponding ensemble.  

Optional flags:    

    - path_rw  : if provided, correlators are reweighted.  
    - L        : correlators are normalised with the volume L^3. L=1 by default.    
    - frw_bcwd : if true, the forward backward symmetrization is performed. 
    - impr     : if true, the correlator is improved.  
    - impr_set : either "1" or "2", to select set1 or set2 of improvement coefficients accordingly.  
    - cons     : if true, the function returns both the local-local and the local-conserved G88 correlators, else only the local-local is returned.  
    - std      : if true, standard symmetric derivatives are used for the vector-tensor correlator. If false, improved derivatives are used. The latter are thought to improve the short distance cutoff effects.
    -lma       : if true, allows for LMA data to be used

Examples:
```@example
g33_ll, g33_lc = corr33(path, ensinfo, path_rw=path_rw, frw_bcwd=true, impr=true, cons=true, std=false) 
```
"""
function corr33(path_data::String, ens::EnsInfo; path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, impr_set::String="1", cons::Bool=true, std::Bool=false, lma::Bool=true)
    
    Gamma_l = ["V1V1", "V2V2", "V3V3", "V1T10", "V2T20", "V3T30"]
    Gamma_c = ["V1V1c", "V2V2c", "V3V3c", "V1cT10", "V2cT20", "V3cT30"]

    v1v1 = get_corr(path_data, ens, "light", Gamma_l[1], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)
    v2v2 = get_corr(path_data, ens, "light", Gamma_l[2], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)
    v3v3 = get_corr(path_data, ens, "light", Gamma_l[3], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)

    vvobs = -(v1v1.obs .+ v2v2.obs .+ v3v3.obs)/3.

    if impr
        v1t10 = get_corr(path_data, ens, "light", Gamma_l[4], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)
        v2t20 = get_corr(path_data, ens, "light", Gamma_l[5], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)
        v3t30 = get_corr(path_data, ens, "light", Gamma_l[6], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)

        vtobs = -(v1t10.obs .+ v2t20.obs .+ v3t30.obs)/3.

        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta)
        elseif impr_set == "1old"
            cv_l = cv_loc_old(beta)
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
        end
        improve_corr_vkvk!(vvobs, vtobs, 2*cv_l, std=std)
    end
    g33 = Corr(0.5 .* vvobs, v1v1.id, "G33ll")
    if frw_bcwd
        frwd_bckwrd_symm!(g33)
    end

    if cons
        v1v1_c = get_corr(path_data, ens, "light", Gamma_c[1], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)
        v2v2_c = get_corr(path_data, ens, "light", Gamma_c[2], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)
        v3v3_c = get_corr(path_data, ens, "light", Gamma_c[3], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)

        vvobs_c = -(v1v1_c.obs .+ v2v2_c.obs .+ v3v3_c.obs)/3.

        if impr
            v1t10_c = get_corr(path_data, ens, "light", Gamma_c[4], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)
            v2t20_c = get_corr(path_data, ens, "light", Gamma_c[5], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)    
            v3t30_c = get_corr(path_data, ens, "light", Gamma_c[6], path_rw=path_rw, frw_bcwd=false, L=L, lma=lma)

            vtobs_c = -(v1t10_c.obs .+ v2t20_c.obs .+ v3t30_c.obs)/3.

            beta = ens.beta
            if impr_set == "1"
                cv_l = cv_loc(beta) 
                cv_c = cv_cons(beta)
            elseif impr_set == "1old"
                cv_l = cv_loc_old(beta)
                cv_c = cv_cons_old(beta)
            elseif impr_set == "2"
                cv_l = cv_loc_set2(beta)
                cv_c = cv_cons_set2(beta)
            end
            improve_corr_vkvk_cons!(vvobs_c, vtobs, vtobs_c, cv_l, cv_c, std=std)
        end
        g33_c = Corr(0.5 .* vvobs_c, v1v1_c.id, "G33lc")
        if frw_bcwd
            frwd_bckwrd_symm!(g33_c)
        end
    end

    !cons ? (return g33) : (return g33, g33_c)
end

@doc raw"""
    corr88_conn(path_data::String, ens::EnsInfo, g33_ll::Corr; g33_lc::Union{Nothing, Corr}=nothing, path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, cons::Bool=true, std::Bool=true)

This function return the G88 connected correlator given path_data and the EnsInfo data type for the ensemble of interest.  

Optional flags:  

    - g33_lc   : either nothing or a local-conserved G33 correlator. This is used to build the local-conserevd G88. It is required only if cons is set to true.
    - path_rw  : if provided, correlators are reweighted.  
    - L        : correlators are normalised with the volume L^3. L=1 by default.    
    - frw_bcwd : if true, the forward backward symmetrization is performed. 
    - impr     : if true, the correlator is improved.  
    - impr_set : either "1" or "2", to select set1 or set2 of improvement coefficients accordingly.  
    - cons     : if true, the function returns both the local-local and the local-conserved G88 correlators, else only the local-local is returned.  
    - std      : if true, standard symmetric derivatives are used for the vector-tensor correlator. If false, improved derivatives are used. The latter are thought to improve the short distance cutoff effects.

Examples:
```@example
corr88_conn(path, ensinfo, path_rw=path_rw, frw_bcwd=true, impr=true, std=true) 
```
"""
function corr88_conn(path_data::String, ens::EnsInfo, g33_ll::Corr; g33_lc::Union{Nothing, Corr}=nothing, path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, impr_set::String="1", cons::Bool=true, std::Bool=true)

    Gamma_l = ["V1V1", "V2V2", "V3V3", "V1T10", "V2T20", "V3T30"]
    Gamma_c = ["V1V1c", "V2V2c", "V3V3c", "V1cT10", "V2cT20", "V3cT30"]

    v1v1 = get_corr(path_data, ens, "strange", Gamma_l[1], path_rw=path_rw, frw_bcwd=false, L=L)
    v2v2 = get_corr(path_data, ens, "strange", Gamma_l[2], path_rw=path_rw, frw_bcwd=false, L=L)
    v3v3 = get_corr(path_data, ens, "strange", Gamma_l[3], path_rw=path_rw, frw_bcwd=false, L=L)

    vvobs = -(v1v1.obs .+ v2v2.obs .+ v3v3.obs)/3.

    if impr
        v1t10 = get_corr(path_data, ens, "strange", Gamma_l[4], path_rw=path_rw, frw_bcwd=false, L=L)
        v2t20 = get_corr(path_data, ens, "strange", Gamma_l[5], path_rw=path_rw, frw_bcwd=false, L=L)
        v3t30 = get_corr(path_data, ens, "strange", Gamma_l[6], path_rw=path_rw, frw_bcwd=false, L=L)

        vtobs = -(v1t10.obs .+ v2t20.obs .+ v3t30.obs)/3.

        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta)
        elseif impr_set == "1old"
            cv_l = cv_loc_old(beta)
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
        end 
        improve_corr_vkvk!(vvobs, vtobs, 2*cv_l, std=std)
    end

    g88_ll = Corr(1/6 .* (2 * g33_ll.obs + 2 * vvobs ), v1v1.id, "G88ll") 
    if frw_bcwd
        frwd_bckwrd_symm!(g88_ll)
    end

    if cons
        if isnothing(g33_lc)
            error("G33_lc is required to compute G88_lc")
        end
        v1v1_c = get_corr(path_data, ens, "strange", Gamma_c[1], path_rw=path_rw, frw_bcwd=false, L=L)
        v2v2_c = get_corr(path_data, ens, "strange", Gamma_c[2], path_rw=path_rw, frw_bcwd=false, L=L)
        v3v3_c = get_corr(path_data, ens, "strange", Gamma_c[3], path_rw=path_rw, frw_bcwd=false, L=L)
    
        vvobs_c = -(v1v1_c.obs .+ v2v2_c.obs .+ v3v3_c.obs)/3.
        
        if impr
            v1t10_c = get_corr(path_data, ens, "strange", Gamma_c[4], path_rw=path_rw, frw_bcwd=false, L=L)
            v2t20_c = get_corr(path_data, ens, "strange", Gamma_c[5], path_rw=path_rw, frw_bcwd=false, L=L)
            v3t30_c = get_corr(path_data, ens, "strange", Gamma_c[6], path_rw=path_rw, frw_bcwd=false, L=L)

            vtobs_c = -(v1t10_c.obs .+ v2t20_c.obs .+ v3t30_c.obs)/3.

            beta = ens.beta
            if impr_set == "1"
                cv_l = cv_loc(beta) 
                cv_c = cv_cons(beta)
            elseif impr_set == "1old"
                cv_l = cv_loc_old(beta)
                cv_c = cv_cons_old(beta)
            elseif impr_set == "2"
                cv_l = cv_loc_set2(beta)
                cv_c = cv_cons_set2(beta)
            end
            improve_corr_vkvk_cons!(vvobs_c, vtobs, vtobs_c, cv_l, cv_c, std=std)
        end
        g88_lc = Corr(1/6 .* (2 * g33_lc.obs + 2 * vvobs_c ), v1v1.id, "G88lc") 
        if frw_bcwd
            frwd_bckwrd_symm!(g88_lc)
        end
    end 

    !cons ? (return g88_ll) : (return g88_ll, g88_lc) 
end

@doc raw"""
corr08_conn(g33_ll::Corr, g88_ll::Corr; g33_lc::Union{Corr, Nothing}=nothing, g88_lc::Union{Corr, Nothing}=nothing)

This function return the local-local G08 connected correlator from the G33 and G88 connected according to:

```math
G_{08}^{conn} = \frac{\sqrt{3}}{2}(G_{33} - G_{88}^{conn})
```

If the local-conserved G33_lc and G88_lc correlators are passed, this functions also returns the local-consereved G08 connected correlator.

Examples:
```@example
g08_ll, g08_lc = corr08_conn(g33, g88, g33_lc=g33_lc, g88_lc=g88_lc) 
```
"""
function corr08_conn(g33_ll::Corr, g88_ll::Corr; g33_lc::Union{Corr, Nothing}=nothing, g88_lc::Union{Corr, Nothing}=nothing)
    
    g08_ll =  Corr(sqrt(3)/2 * (g33_ll.obs - g88_ll.obs), g33_ll.id, "G08ll")
    if !isnothing(g33_lc) && !isnothing(g88_lc)
        g08_lc =  Corr(sqrt(3)/2 * (g33_lc.obs - g88_lc.obs), g33_lc.id, "G08lc")
        return g08_ll, g08_lc
    else
        return g08_ll
    end
end

@doc raw"""
    corrC_conn(path_data::String, ens::EnsInfo; path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, impr_set::String="1", cons::Bool=true, std::Bool=true)

This function return the G charmed correlator given path_data and the EnsInfo data type of the corresponding ensemble.  

Optional flags:    

    - path_rw  : if provided, correlators are reweighted.  
    - L        : correlators are normalised with the volume L^3. L=1 by default.    
    - frw_bcwd : if true, the forward backward symmetrization is performed. 
    - impr     : if true, the correlator is improved.  
    - impr_set : either "1" or "2", to select set1 or set2 of improvement coefficients accordingly.  
    - cons     : if true, the function returns both the local-local and the local-conserved G88 correlators, else only the local-local is returned.  
    - std      : if true, standard symmetric derivatives are used for the vector-tensor correlator. If false, improved derivatives are used. The latter are thought to improve the short distance cutoff effects.

Examples:
```@example
gC_ll, gC_lc = corrC_conn(path, ensinfo, path_rw=path_rw, frw_bcwd=true, impr=true, cons=true, std=false) 
```
"""
function corrC_conn(path_data::String, ens::EnsInfo; path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, impr_set::String="1", cons::Bool=true, std::Bool=false, plus::Bool=false)

    cstr = plus ? "charm_plus" : "charm"
    
    Gamma_l = ["V1V1", "V2V2", "V3V3", "V1T10", "V2T20", "V3T30"]
    Gamma_c = ["V1V1c", "V2V2c", "V3V3c", "V1cT10", "V2cT20", "V3cT30"]

    v1v1 = get_corr(path_data, ens, cstr, Gamma_l[1], path_rw=path_rw, frw_bcwd=false, L=L)
    v2v2 = get_corr(path_data, ens, cstr, Gamma_l[2], path_rw=path_rw, frw_bcwd=false, L=L)
    v3v3 = get_corr(path_data, ens, cstr, Gamma_l[3], path_rw=path_rw, frw_bcwd=false, L=L)
    
    vvobs = -(v1v1.obs .+ v2v2.obs .+ v3v3.obs)/3.

    if impr
        v1t10 = get_corr(path_data, ens, cstr, Gamma_l[4], path_rw=path_rw, frw_bcwd=false, L=L)
        v2t20 = get_corr(path_data, ens, cstr, Gamma_l[5], path_rw=path_rw, frw_bcwd=false, L=L)
        v3t30 = get_corr(path_data, ens, cstr, Gamma_l[6], path_rw=path_rw, frw_bcwd=false, L=L)

        vtobs = -(v1t10.obs .+ v2t20.obs .+ v3t30.obs)/3.

        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta)
        elseif impr_set == "1old"
            cv_l = cv_loc_old(beta)
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
        end
        improve_corr_vkvk!(vvobs, vtobs, 2*cv_l, std=std)
    end
    gC = Corr(vvobs, v1v1.id, "GCll")
    if frw_bcwd
        frwd_bckwrd_symm!(gC)
    end

    if cons
        v1v1_c = get_corr(path_data, ens, cstr, Gamma_c[1], path_rw=path_rw, frw_bcwd=false, L=L)
        v2v2_c = get_corr(path_data, ens, cstr, Gamma_c[2], path_rw=path_rw, frw_bcwd=false, L=L)
        v3v3_c = get_corr(path_data, ens, cstr, Gamma_c[3], path_rw=path_rw, frw_bcwd=false, L=L)
    
        vvobs_c = -(v1v1_c.obs .+ v2v2_c.obs .+ v3v3_c.obs)/3.

        if impr
            v1t10_c = get_corr(path_data, ens, cstr, Gamma_c[4], path_rw=path_rw, frw_bcwd=false, L=L)
            v2t20_c = get_corr(path_data, ens, cstr, Gamma_c[5], path_rw=path_rw, frw_bcwd=false, L=L)    
            v3t30_c = get_corr(path_data, ens, cstr, Gamma_c[6], path_rw=path_rw, frw_bcwd=false, L=L)

            vtobs_c = -(v1t10_c.obs .+ v2t20_c.obs .+ v3t30_c.obs)/3.

            beta = ens.beta
            if impr_set == "1"
                cv_l = cv_loc(beta) 
                cv_c = cv_cons(beta)
            elseif impr_set == "1old"
                cv_l = cv_loc_old(beta)
                cv_c = cv_cons_old(beta)
            elseif impr_set == "2"
                cv_l = cv_loc_set2(beta)
                cv_c = cv_cons_set2(beta)
            end
            improve_corr_vkvk_cons!(vvobs_c, vtobs, vtobs_c, cv_l, cv_c, std=std)
        end
        GC_c = Corr(vvobs_c, v1v1_c.id, "GClc")
        if frw_bcwd
            frwd_bckwrd_symm!(GC_c)
        end
    end

    !cons ? (return gC) : (return gC, GC_c)
end

