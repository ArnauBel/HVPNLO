#============= Integrand structure =============#

struct Integrand
    obs::Union{Vector{uwreal},Matrix{uwreal}}
    id::String
    gamma::String

    function Integrand(a::Union{Vector{uwreal},Matrix{uwreal}}, cd::CData)
       
        return new(a, cd.id, cd.gamma)
    end
    Integrand(a::Union{Vector{uwreal},Matrix{uwreal}}, id::String, gamma::String) = new(a, id, gamma) 
end
function Base.show(io::IO, integrand::Integrand)
    println(io, "Integrand")
    println(io, " - Ensemble ID: ", integrand.id)
    println(io, " - Gamma:       ", integrand.gamma)
end
export Integrand

## Label function

function label_func(string::String)

    mystring = collect(string)
    body = mystring[2]*mystring[3]
    tail = mystring[end-1]*mystring[end]
    nnum = mystring[4]

    length(mystring)==4 ? (return latexstring("G^{(R)}_{$tail}(t)")) : nothing

    length(mystring)==5 ? (return latexstring("G^{($body)}_{$tail}(t)")) : nothing

    length(mystring)==6 ? (return latexstring("G^{($body)}_{$tail}(t)\\ \\vec{n}^2\\leq$nnum")) : nothing
end

#============= HVP NLO computation =============#

function amuHVPNLO(diagram::String, ens::EnsInfo; path_HVP::String=joinpath(@__DIR__, "..", "LatticeData", "HVP_data"), path_rw::String=joinpath(@__DIR__, "..", "LatticeData", "rwf_deflated"), path_ms::String=joinpath(@__DIR__, "..", "LatticeData", "ms_t0_dat"), path_coef::String=joinpath(@__DIR__, "..", "Coefficients"), pl::Bool=false, errmult::Union{Int64,Float64}=1)::Tuple{uwreal, uwreal}

    gll, glc = corrR(path_HVP, ens, path_rw = path_rw, frw_bcwd = true)

    sym_points = Int64(length(gll.obs)/2+1)

    # t0_ph = uwreal([0.1439, 0.0006], "sqrtt0 [fm]")
    t0ens = get_t0(path_ms, ens, path_rw = path_rw, pl = false)
    aens = t0_ph ./ sqrt.( t0ens )
    
    # alpha = 1/137.035999084
    # massmu = 0.10565837
    # hbarc = 0.1973269804        # GeV * fm

    factor = hbarc * sqrt(t0ens)/t0_ph     # convertion from GeV to 1/a units
    tfm =  aens .* collect(0:sym_points-1)

    functions_dict =  Dict("a" => Tildef4a, "b" => Tildef4b, "c" => Tildef4c)

    if haskey(functions_dict, diagram)
        diagram == "c" ? (Tildef = factor^4 .* functions_dict[diagram]((massmu/factor) .* collect(0:sym_points-1), path_coef)) : (Tildef = factor^2 .* functions_dict[diagram]((massmu/factor) .* collect(0:sym_points-1), path_coef))     # kernel in lattice unitis
    else
        error("Invalid diagram: '$diagram'. \n Please choose between diagrams 'a', 'b' or 'c'. ")
    end

    if diagram in ["a", "b"]
        intll = - gll.obs[1:sym_points] .* Tildef
        intlc = - glc.obs[1:sym_points] .* Tildef

        amuNLOll = (alpha/pi)^3 * sum(intll) * 1e10
        amuNLOlc = (alpha/pi)^3 * sum(intlc) * 1e10

    else
        intll = (gll.obs[1:sym_points] .* hcat(gll.obs[1:sym_points]...)) .* Tildef
        intlc = (glc.obs[1:sym_points] .* hcat(glc.obs[1:sym_points]...)) .* Tildef

        diagll, diaglc = [diag(intll), diag(intlc)]
        triagll, triaglc = [tril(intll, -1), tril(intlc, -1)]

        amuNLOll = (alpha/pi)^3 * (sum(diagll) + 2*sum(triagll)) * 1e10
        amuNLOlc = (alpha/pi)^3 * (sum(diaglc) + 2*sum(triaglc)) * 1e10

    end

    uwerr(amuNLOll)
    uwerr(amuNLOlc)

    if pl
        uwerr.(intll)
        uwerr.(intlc)
        if diagram in ["a", "b"]
            errorbar(collect(0:sym_points-1).-0.2, value.(intlc), ADerrors.err.(intlc).*errmult, fmt="s", label="Local-conserved integrand", color = "red", capsize=2)
            errorbar(collect(0:sym_points-1).+0.2, value.(intll), ADerrors.err.(intll).*errmult, fmt="s", label="Local-local integrand", color = "green", capsize=2)
            axis("tight")
            errmult == 1 ? (PyPlot.title("Integrand for diagram (" * diagram * "); Ensamble = " * ens.id)) : (PyPlot.title("Integrand for diagram (" * diagram * "); Ensamble = " * ens.id * " (Errors x $errmult)"))
            xlabel(L"$t$ [a]")
            ylabel(L"$G(t)$ $\tilde{f}_4(m_\mu t)$")
            legend(["Local-local integrand","Local-conserved integrand"], loc  = "best")
            #grid("on")
            display(gcf())      #display the figure
            close("all")
        else
            figc, axc = subplots(figsize=(8, 8))
            myheatmap = axc.pcolormesh(value.(tfm), value.(tfm), value.(intlc), cmap="Reds")
            colorbar(myheatmap, label=L"$\tilde{f}_4^{(c)}(t,\tau)$")
            xlabel(L"$t$ [fm]")
            ylabel(L"$\tau$ [fm]")
            title("Integrand for diagram (c) & local-conserved currents; Ensamble = " * ens.id)
            display(gcf())      #display the figure
            close()
        end
    end

    return amuNLOll, amuNLOlc
end
function amuHVPNLO(diagram::String, corr::Union{Corr,Vector{Corr}}, t0ens::uwreal; path_coef::String=joinpath(@__DIR__, "..", "Coefficients"), pl::Bool=false, int::Bool=false)

    typeof(corr) == Corr ? (corr=[corr]) : nothing
    
    # t0_ph = uwreal([0.1439, 0.0006], "sqrtt0 [fm]")
    aens = t0_ph ./ sqrt.(t0ens)

    lencheck = length(corr[1].obs)
    for subcorr in corr
        length(subcorr.obs) != lencheck ? (error("All introduced correlators sohuld be the same length")) : nothing
    end

    sym_points = Int64(lencheck/2+1)
    
    # alpha = 1/137.035999084
    # massmu = 0.10565837
    # hbarc = 0.1973269804        # GeV * fm

    factor = hbarc * sqrt(t0ens)/t0_ph     # convertion from GeV to 1/a units
    tfm =  aens .* collect(0:sym_points-1)

    functions_dict =  Dict("a" => Tildef4a, "b" => Tildef4b, "c" => Tildef4c)

    if haskey(functions_dict, diagram)
        diagram == "c" ? (Tildef = factor^4 .* functions_dict[diagram]((massmu/factor) .* collect(0:sym_points-1), path_coef)) : (Tildef = factor^2 .* functions_dict[diagram]((massmu/factor) .* collect(0:sym_points-1), path_coef)) 
    else
        error("Invalid diagram: '$diagram'. \n Please choose between diagrams 'a', 'b' or 'c'. ")
    end

    integrand = Integrand[]
    amuNLO = uwreal[]

    if diagram in ["a", "b"]
        for subcorr in corr
            subintegrand = - subcorr.obs[1:sym_points] .* Tildef
            subamuNLO = (alpha/pi)^3 * sum(subintegrand) * 1e10

            push!(integrand, Integrand(subintegrand,subcorr.id,subcorr.gamma))
            push!(amuNLO, subamuNLO)
        end
    else
        for subcorr in corr
            subintegrand = (subcorr.obs[1:sym_points] .* hcat(subcorr.obs[1:sym_points]...)) .* Tildef
            subamuNLO = (alpha/pi)^3 * (sum(diag(subintegrand)) + 2*sum(tril(subintegrand, -1))) * 1e10

            push!(integrand, Integrand(subintegrand,subcorr.id,subcorr.gamma))
            push!(amuNLO, subamuNLO)
        end
    end
    
    uwerr.(amuNLO)

    if pl
        if diagram in ["a", "b"]
            label=[]
            for subintegrand in integrand
                uwerr.(subintegrand.obs)
                fill_between(value.(tfm), value.(subintegrand.obs)-ADerrors.err.(subintegrand.obs), value.(subintegrand.obs)+ADerrors.err.(subintegrand.obs), alpha=0.8, label=label_func(subintegrand.gamma))
                push!(label,label_func(subintegrand.gamma))
            end
            axis("tight")
            PyPlot.title(ens.id*": Integrands for diagram (" * diagram * ")")
            xlabel(L"$t$ [fm]")
            diagram == "a" ? (ylabel(L"$G(t)$ $\tilde{f}^{(a)}_4(m_\mu t)$")) : (ylabel(L"$G(t)$ $\tilde{f}^{(b)}_4(m_\mu t)$"))
            legend(label, loc  = "best")
            #grid("on")
            display(gcf())      #display the figure
            close("all")
        else
            for subintegrand in integrand
                figc, axc = subplots(figsize=(8, 8))
                myheatmap = axc.pcolormesh(value.(tfm), value.(tfm), value.(subintegrand.obs), cmap="Reds")
                colorbar(myheatmap, label=L"$\tilde{f}_4^{(c)}(t,\tau)$")
                xlabel(L"$t$ [fm]")
                ylabel(L"$\tau$ [fm]")
                title(ens.id*": Integrand for diagram (c); " * label_func(subintegrand.gamma))
                display(gcf())      #display the figure
                close()
            end
            
        end
    end

    int ? (return amuNLO, integrand) : (return amuNLO)
end
function amuHVPNLO(diagram::String, corr::Union{Corr,Vector{Corr}}, Tildef::Union{Vector{uwreal},Matrix{uwreal}}; t0ens::Union{uwreal,Nothing}=nothing, pl::Bool=false, int::Bool=false)

    typeof(corr) == Corr ? (corr=[corr]) : nothing

    lencheck = length(corr[1].obs)
    for subcorr in corr
        length(subcorr.obs) != lencheck ? (error("All introduced correlators sohuld be the same length")) : nothing
    end

    pl == true && isnothing(t0ens) ? (error("t0ens is required to perform the plot")) : nothing

    sym_points = Int64(lencheck/2+1)
    
    # alpha = 1/137.035999084

    integrand = Integrand[]
    amuNLO = uwreal[]

    if diagram in ["a", "b"]
        for subcorr in corr
            subintegrand = - subcorr.obs[1:sym_points] .* Tildef
            subamuNLO = (alpha/pi)^3 * sum(subintegrand) * 1e10

            push!(integrand, Integrand(subintegrand,subcorr.id,subcorr.gamma))
            push!(amuNLO, subamuNLO)
        end
    else
        for subcorr in corr
            subintegrand = (subcorr.obs[1:sym_points] .* hcat(subcorr.obs[1:sym_points]...)) .* Tildef
            subamuNLO = (alpha/pi)^3 * (sum(diag(subintegrand)) + 2*sum(tril(subintegrand, -1))) * 1e10

            push!(integrand, Integrand(subintegrand,subcorr.id,subcorr.gamma))
            push!(amuNLO, subamuNLO)
        end
    end
    
    uwerr.(amuNLO)

    if pl
        # t0_ph = uwreal([0.1439, 0.0006], "sqrtt0 [fm]")
        aens = t0_ph ./ sqrt.( t0ens )
        tfm =  aens .* collect(0:sym_points-1)
        if diagram in ["a", "b"]
            label=[]
            for subintegrand in integrand
                uwerr.(subintegrand.obs)
                if subintegrand.gamma in ["G08ll","G08lc"]
                    if !all(x -> x == 0.0, value.(subintegrand.obs))
                        fill_between(value.(tfm), -value.(subintegrand.obs)+ADerrors.err.(subintegrand.obs), -value.(subintegrand.obs)-ADerrors.err.(subintegrand.obs), alpha=0.8, label="-"*label_func(subintegrand.gamma))
                        push!(label,"-"*label_func(subintegrand.gamma))
                    else
                        nothing
                    end
                else
                    fill_between(value.(tfm), value.(subintegrand.obs)-ADerrors.err.(subintegrand.obs), value.(subintegrand.obs)+ADerrors.err.(subintegrand.obs), alpha=0.8, label=label_func(subintegrand.gamma))
                    push!(label,label_func(subintegrand.gamma))
                end
            end
            axis("tight")
            PyPlot.title(ens.id*": Integrands for diagram (" * diagram * ")")
            xlabel(L"$t$ [fm]")
            diagram == "a" ? (ylabel(L"$G(t)$ $\tilde{f}^{(a)}_4(m_\mu t)$")) : (ylabel(L"$G(t)$ $\tilde{f}^{(b)}_4(m_\mu t)$"))
            legend(label, loc  = "best")
            #grid("on")
            display(gcf())      #display the figure
            close("all")
        else
            for subintegrand in integrand
                figc, axc = subplots(figsize=(8, 8))
                myheatmap = axc.pcolormesh(value.(tfm), value.(tfm), value.(subintegrand.obs), cmap="Reds")
                colorbar(myheatmap, label=L"$\tilde{f}_4^{(c)}(t,\tau)$")
                xlabel(L"$t$ [fm]")
                ylabel(L"$\tau$ [fm]")
                title(ens.id*": Integrand for diagram (c); " * label_func(subintegrand.gamma))
                display(gcf())      #display the figure
                close()
            end
            
        end
    end

    int ? (return amuNLO, integrand) : (return amuNLO)
end

#============= Finite volume corrections =============#

function amu∆G(diagram::String, ∆corr::Union{Corr,Vector{Corr}}, Tildef::Union{Vector{uwreal},Matrix{uwreal}};  corr::Union{Corr,Nothing}=nothing, t0ens::Union{uwreal,Nothing}=nothing, pl::Bool=false, int::Bool=false)

    lencheck = length(∆corr[1].obs)
    idcheck = ∆corr[1].id
    for sub∆corr in ∆corr
        length(sub∆corr.obs) != lencheck ? (error("All introduced fv corrections sohuld be the same length")) : nothing
        typeof(corr) == Corr ? (length(corr.obs)/2 != lencheck ? (error("The corrrelator and the fv corrections sohuld be the same length")) : nothing) : nothing

        sub∆corr.id != idcheck ? (error("All fv corrections should correspond to the same ensamble")) : nothing
        typeof(corr) == Corr ? (corr.id != idcheck ? (error("The corrrelator and the fv corrections sohuld come from the same ensamble")) : nothing) : nothing
    end

    pl == true && isnothing(t0ens) ? (error("t0ens is required to perform the plot")) : nothing
    diagram == "c" && isnothing(corr) ? (error("For the FV corrections of diagram 'c' the correlator must be given")) : nothing

    ens = EnsInfo(∆corr[1].id)
    ens.kappa_l == ens.kappa_s ? (mult=1.5) : (mult=1.)
    
    typeof(∆corr) == Corr ? (∆corr=[∆corr]) : nothing

    sym_points = Int64(lencheck+1)
    # alpha = 1/137.035999084

    integrand = Integrand[]
    amuNLO = uwreal[]

    if diagram in ["a", "b"]
        for subcorr in ∆corr
            subintegrand = - mult .* vcat(Tildef[1], subcorr.obs .* Tildef[2:end])
            subamuNLO = (alpha/pi)^3 * sum(subintegrand) * 1e10

            push!(integrand, Integrand(subintegrand,subcorr.id,subcorr.gamma))
            push!(amuNLO, subamuNLO)
        end
        if typeof(corr) == Corr
            background = Vector{Float64}
            corrint = - corr.obs[1:sym_points] .* Tildef
            uwerr.(corrint)
            
            diagram == "a" ? (background=-ADerrors.err.(corrint)) : (background=ADerrors.err.(corrint))
        end
        
    else
        for subcorr in ∆corr
            fullsubcorr = vcat(Tildef[1,1],subcorr.obs)
            subintegrand = (mult^2 .* (fullsubcorr .* hcat(fullsubcorr...)) .+ mult .* (corr.obs[1:sym_points] .* hcat(fullsubcorr...) .+ fullsubcorr .* hcat(corr.obs[1:sym_points]...))) .* Tildef
            subamuNLO = (alpha/pi)^3 * (sum(diag(subintegrand)) + 2*sum(tril(subintegrand, -1))) * 1e10

            push!(integrand, Integrand(subintegrand,subcorr.id,subcorr.gamma))
            push!(amuNLO, subamuNLO)
        end
        if typeof(corr) == Corr
            background = Matrix{Float64}
            corrint = (corr.obs[1:sym_points] .* hcat(corr.obs[1:sym_points]...)) .* Tildef
            uwerr.(corrint)

            background = ADerrors.err.(corrint)
        end
    end
    
    uwerr.(amuNLO)

    if pl
        # t0_ph = uwreal([0.1439, 0.0006], "sqrtt0 [fm]")
        aens = t0_ph ./ sqrt.( t0ens )
        tfm =  aens .* collect(0:sym_points-1)
        if diagram in ["a", "b"]
            label=[]
            if typeof(corr) == Corr
                fill_between(value.(tfm), 0, background, alpha=0.2, color="gray", label=L"$\delta G(t)$")
                push!(label,L"$\delta G(t)$")
            end
            for subintegrand in integrand
                uwerr.(subintegrand.obs)
                fill_between(value.(tfm), value.(subintegrand.obs)-ADerrors.err.(subintegrand.obs), value.(subintegrand.obs)+ADerrors.err.(subintegrand.obs), alpha=0.4, label=label_func(subintegrand.gamma))
                push!(label,label_func(subintegrand.gamma))
            end
            axis("tight")
            PyPlot.title(ens.id*": FV corrections for diagram (" * diagram * ")")
            xlabel(L"$t$ [fm]")
            diagram == "a" ? (ylabel(L"$\Delta G(t)$ $\tilde{f}^{(a)}_4(m_\mu t)$")) : (ylabel(L"$\Delta G(t)$ $\tilde{f}^{(b)}_4(m_\mu t)$"))
            legend(label, loc  = "best")
            #grid("on")
            display(gcf())      #display the figure
            close("all")
        else
            figc, axc = subplots(figsize=(8, 8))
            myheatmap = axc.pcolormesh(value.(tfm), value.(tfm), value.(integrand[end].obs), cmap="Reds")
            colorbar(myheatmap, label=L"$\Delta \tilde{f}_4^{(c)}(t,\tau)$")
            xlabel(L"$t$ [fm]")
            ylabel(L"$\tau$ [fm]")
            title(ens.id*": FV corrections for diagram (c); "*label_func(integrand[end].gamma))
            display(gcf())      #display the figure
            close()            
        end
    end

    int ? (return amuNLO, integrand) : (return amuNLO)
end


#============= Bounding Method  ==============#

function Eeff(tstar::Int64, obs::Vector{uwreal})

    E_eff = 0.5 * log((obs[tstar+1] / obs[tstar+2]) ^2)
    return E_eff
end
function Eeff(tstar::Int64, obs::Matrix{uwreal})

    E_eff = 0.5 * log((obs[tstar+1] / obs[tstar+2]) ^2)
    return E_eff
end

function Gb(t::Vector{Int64}, tcut::Int64, corr::Corr, Eeff::uwreal)
    T = 2*t[end]

    Gcut = corr.obs[tcut+1]
    GPBC(x0) = exp(-Eeff*x0) + exp(-Eeff*(T-x0))
    G(x0) = exp(-Eeff*x0)

    split(corr.id,"")[end-1]=="5" ? (UBarray = (Gcut/GPBC(tcut)) .*  GPBC.(t)) : (UBarray = (Gcut/G(tcut)) .*  G.(t))
    return UBarray[tcut+1:end]
end
function Gb(t::Vector{Int64}, tcut::Int64, obs::Vector{uwreal}, ens::EnsInfo, Eeff::uwreal)
    T = 2*t[end]

    Gcut = obs[tcut+1]
    GPBC(x0) = exp(-Eeff*x0) + exp(-Eeff*(T-x0))
    G(x0) = exp(-Eeff*x0)

    split(ens.id,"")[end-1]=="5" ? (UBarray = (Gcut/GPBC(tcut)) .*  GPBC.(t)) : (UBarray = (Gcut/G(tcut)) .*  G.(t))
    return UBarray[tcut+1:end]
end

function amu_BM(diagram::String, corr::Vector{Corr}, Tildef::Union{Vector{uwreal},Matrix{uwreal}}; bound_impr::Bool=false, t_step::Int64=1, t0ens::Union{uwreal,Nothing}=nothing, pl::Bool=false)

    pl && isnothing(t0ens) ? (error("t0ens is required to perform the plot")) : nothing

    lencheck = length(corr[1].obs)
    any(subcorr -> length(subcorr.obs) != lencheck, corr[2:end]) ? error("All introduced correlators should be the same length") : nothing

    tcut_initial, sym_points = Int64(length(corr[1].obs)/8+1), Int64(length(corr[1].obs)/2+1)
    t = collect(0:sym_points-1)
    ens = EnsInfo(corr[1].id)

    # alpha = 1/137.035999084

    mycorr33 = corr[1].obs
    mycorr88 = corr[2].obs
    mycorr08 = corr[3].obs
    mycorrR  = corr[end].obs

    mpi = m_ens[ens.id]["m_pi"] 
    mrho = m_ens[ens.id]["m_rho"]
    L = ens.L
    E2pi = 2*sqrt(mpi^2 + (2*pi/L)^2)

    ub = uwreal[]
    lb = uwreal[]
    bound_impr ? (lb_impr = uwreal[]) : nothing
    
    amuNLO = uwreal[]

    if diagram in ["a", "b"]
        int = - mycorrR[1:sym_points] .* Tildef 
        for tcut in tcut_initial:t_step:(sym_points-t_step)

            mrho < E2pi ? (UB33 = Gb(t, tcut, mycorr33, ens, mrho)) : (UB33 = Gb(t, tcut, mycorr33, ens, E2pi))
            UB88 = Gb(t, tcut, mycorr88, ens, mrho)
            value(mycorr08[1]) != 0.0 ? (UB08 = Gb(t, tcut, mycorr08, ens, mrho)) : (UB08 = fill(0.,length(UB33)))

            UBInt = (UB33 .+ UB88 .+ UB08) .* Tildef[tcut+1:end]

            ub_amuNLO = (alpha/pi)^3 * (sum(int[1:tcut])+sum(UBInt)) * 1e10
            lb_amuNLO = (alpha/pi)^3 * sum(int[1:tcut]) * 1e10

            push!(lb, lb_amuNLO)
            push!(ub, ub_amuNLO)

            if bound_impr
                E_eff33=Eeff(tcut, mycorr33)
                E_eff88=Eeff(tcut, mycorr88)
                E_eff08=Eeff(tcut, mycorr08)

                LB33 = Gb(t, tcut, mycorr33, ens, E_eff33)
                LB88 = Gb(t, tcut, mycorr88, ens, E_eff88)
                value(mycorr08[1]) != 0.0 ? (LB08 = Gb(t, tcut, mycorr08, ens, E_eff08)) : (LB08 = fill(0.,length(UB33)))

                LBint = (LB33 .+ LB88 .+ LB08) .* Tildef[tcut+1:end]
                lb_impr_amuNLO = (alpha/pi)^3 * (sum(int[1:tcut])+sum(LBint)) * 1e10

                push!(lb_impr, lb_impr_amuNLO)
            end
        end
        amuNLO = (alpha/pi)^3 * sum(int) * 1e10
    end
    if diagram == "c"
        int = (mycorrR[1:sym_points] .* hcat(mycorrR[1:sym_points]...)) .* Tildef
        for tcut in tcut_initial:t_step:(sym_points-t_step)
            
            mrho < E2pi ? (UB33 = Gb(t, tcut, mycorr33, ens, mrho)) : (UB33 = Gb(t, tcut, mycorr33, ens, E2pi))
            UB88 = Gb(t, tcut, mycorr88, ens, mrho)
            value(mycorr08[1]) != 0.0 ? (UB08 = Gb(t, tcut, mycorr08, ens, mrho)) : (UB08 = fill(0.,length(UB33)))

            UBInt = ((UB33 .+ UB88 .+ UB08) .* hcat(UB33 .+ UB88 .+ UB08)) .* Tildef[tcut+1:end,tcut+1:end]

            ub_amuNLO = (alpha/pi)^3 * (sum(int[1:tcut,1:tcut]) + 2*sum(int[1:tcut,tcut:end]) + sum(UBInt)) * 1e10
            lb_amuNLO = (alpha/pi)^3 * (sum(int[1:tcut,1:tcut]) + 2*sum(int[1:tcut,tcut:end])) * 1e10

            #(alpha/pi)^3 * (sum(diag(int[1:tcut,1:tcut])) + 2*sum(tril(int[1:tcut,1:tcut], -1)) + 2*sum(int[1:tcut,tcut:end]) + sum(diag(UBInt)) + 2*sum(tril(UBInt, -1))) * 1e10
            #(alpha/pi)^3 * (sum(diag(int[1:tcut,1:tcut])) + 2*sum(tril(int[1:tcut,1:tcut], -1)) + 2*sum(int[1:tcut,tcut:end])) * 1e10

            push!(lb, lb_amuNLO)
            push!(ub, ub_amuNLO)

            if bound_impr
                E_eff33=Eeff(tcut, mycorr33)
                E_eff88=Eeff(tcut, mycorr88)
                E_eff08=Eeff(tcut, mycorr08)

                LB33 = Gb(t, tcut, mycorr33, ens, E_eff33)
                LB88 = Gb(t, tcut, mycorr88, ens, E_eff88)
                value(mycorr08[1]) != 0.0 ? (LB08 = Gb(t, tcut, mycorr08, ens, E_eff08)) : (LB08 = fill(0.,length(UB33)))

                LBint = ((LB33 .+ LB88 .+ LB08) .* hcat(LB33 .+ LB88 .+ LB08)) .* Tildef[tcut+1:end,tcut+1:end]
                lb_impr_amuNLO = (alpha/pi)^3 * (sum(int[1:tcut,1:tcut]) + 2*sum(int[1:tcut,tcut:end]) + sum(LBint)) * 1e10

                #(alpha/pi)^3 * (sum(diag(int[1:tcut,1:tcut])) + 2*sum(tril(int[1:tcut,1:tcut], -1)) + 2*sum(int[1:tcut,tcut:end]) + sum(diag(LBint)) + 2*sum(tril(LBint, -1))) * 1e10

                push!(lb_impr, lb_impr_amuNLO)
            end
        end
        amuNLO = (alpha/pi)^3 * (sum(diag(int)) + 2*sum(tril(int, -1))) * 1e10
    end

    if pl
        uwerr.(ub)
        uwerr.(lb)

        uwerr(amuNLO)

        t0_ph = uwreal([0.1439, 0.0006], "sqrtt0 [fm]")
        aens = t0_ph / sqrt(t0ens)
        tcut_fm =  aens .* collect(tcut_initial:t_step:sym_points-t_step)
        label = ["No bounding",L"Upper bound ($E_0$)","Lower bound (zero)"]

        fill_between(value.(tcut_fm), value(amuNLO)-ADerrors.err(amuNLO), value(amuNLO)+ADerrors.err(amuNLO), alpha=0.2, color="gray", label=label[1])
        errorbar(value.(tcut_fm), value.(ub), ADerrors.err.(ub), marker="o" ,fmt="s", mfc="none", label=label[2], color="red", capsize=2)
        errorbar(value.(tcut_fm), value.(lb), ADerrors.err.(lb), marker="d", fmt="s", mfc="none", label=label[3], color="green", capsize=2)
        if bound_impr
            uwerr.(lb_impr)
            averb = (ub+lb_impr)/2
            uwerr.(averb)

            errorbar(value.(tcut_fm), value.(lb_impr), ADerrors.err.(lb_impr), marker="d" ,fmt="s", mfc="none", label=L"Lower bound ($E_{eff}$)", color="limegreen", capsize=2)
            push!(label, L"Lower bound ($E_{eff}$)")
            errorbar(value.(tcut_fm), value.(averb), ADerrors.err.(averb), marker="s" ,fmt="s", mfc="none", label="Average", color="purple", capsize=2)
            push!(label, "Average")
        end
        
        axis("tight")
        PyPlot.title(ens.id*": Bounding Method (" * diagram * ")")
        xlabel(L"$t_{cut}$ [fm]")
        ylabel(L"$a_\mu^{HVP}[NLO](T/2)$")
        ylim((value(amuNLO)-6*ADerrors.err(amuNLO), value(amuNLO)+4*ADerrors.err(amuNLO)))
        legend(label, loc="best")
        display(gcf())      #display the figure
        close("all")
    end
     
    bound_impr ? (return amuNLO, lb, lb_impr, ub) : (return amuNLO, lb, ub)
end
