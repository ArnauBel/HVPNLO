
using IterTools
using ADerrors
import ADerrors: err

include("HVPtools/Const.jl")

# Permutated terms:

a2(x)       = x[:,1]                                                              # (a^2/8t0)
a2loga(x)   = x[:,1] .* log.(x[:,1])                                              # (a^2/8t0)*log(a^2/8t0)
a3(x)       = x[:,1].^(3/2)                                                       # (a^2/8t0)^{3/2}
a4(x)       = x[:,1].^(2)                                                         # (a^2/8t0)^{2}

a2phi2(x)   = x[:,1] .* (x[:,2] .- value(phi2_ph))                                # (a^2/8t0)*(ϕ2 - ϕ2^{ph})     
a2phi4(x)   = x[:,1] .* (x[:,3] .- value.(phi4_ph))                               # (a^2/8t0)*(ϕ4 - ϕ4^{ph})

a3phi2(x)   = x[:,1].^(3/2) .* (x[:,2] .- value.(phi2_ph))                        # (a^2/8t0)^{3/2}*(ϕ2 - ϕ2^{ph})   

phi2(x)     = (x[:,2] .- value.(phi2_ph))
phi2sqr(x)  = (x[:,2].^2 .- value.(phi2_ph).^2)                                   # (ϕ2^2 - (ϕ2^{ph})^2) 
phi2inv(x)  = (x[:,2].^(-1) .- value.(phi2_ph).^(-1))                             # (1/ϕ2 - 1/(ϕ2^{ph}))
phi2log(x)  = x[:,2] .* log.(x[:,2]) .- value.(phi2_ph) .* log.(value.(phi2_ph))  # (ϕ2log(ϕ2) - ϕ2^{ph}log(ϕ2^{ph}))
logphi2(x)  = log.(x[:,2]) .- log.(value.(phi2_ph))                               # (log(ϕ2) - log(ϕ2^{ph}))

phi4(x)     = (x[:,3] .- value.(phi4_ph))
phi4sqr(x)  = (x[:,3].^2 .- value.(phi4_ph).^2)                                   # (ϕ4^2 - (ϕ4^{ph})^2) 
phi4inv(x)  = (x[:,3].^(-1) .- value.(phi4_ph).^(-1))                             # (1/ϕ4 - 1/(ϕ4^{ph}))
phi4log(x)  = x[:,3] .* log.(x[:,3]) .- value.(phi4_ph) .* log.(value.(phi4_ph))  # (ϕ4log(ϕ4) - ϕ4^{ph}log(ϕ4^{ph}))
logphi4(x)  = log.(x[:,3]) .- log.(value.(phi4_ph))                               # (log(ϕ4) - log(ϕ4^{ph}))


# Multiplier terms

deltaphi(x) = x[:,3] .- 3/2 .* x[:,2]                                             # (ϕ4 - (3/2)*ϕ2)



# Model generator:

function generate_model_map(dimensions::Int)
    map = BitVector[]
    for t in IterTools.product([0:1 for _ in 1:dimensions]...)
        push!(map,BitVector(t))
    end
    return map
end

# Call all model permutations

function call_models(model_var_list::Vector{Function},
    EnsList::Vector{EnsInfo},
    mdof::Int64;
    SimpleBase::Bool=false,
    MultFunc::Union{Nothing,Function}=nothing,
    RESCAL::Bool=true,
    nphi2max::Int64=1,
    nphi4max::Int64=1,
    namax::Int64=3
    )

    nens = length(EnsList)

    isov_basemodel(x,p) = SimpleBase ? (p[1] .+ p[2] .* a2(x) .+ p[3] .* phi2(x)) : (p[1] .+ p[2] .* a2(x) .+ p[3] .* phi2(x) .+ p[4] .* phi4(x))
    if RESCAL
        beta_list = getfield.(EnsList,:beta)
        isov_basemodel_resc(x,p) = SimpleBase ? (p[1] .+ p[2] .* a2_rescaling.(beta_list) .* a2(x) .+ p[3] .* phi2(x)) : (p[1] .+ p[2] .* a2_rescaling.(beta_list) .* a2(x) .+ p[3] .* phi2(x) .+ p[4] .* phi4(x))
    end
    n_par_tot_isov_base = SimpleBase ? 3 : 4

    model_var_label = string.(model_var_list)
    model_map = generate_model_map(length(model_var_list))

    # Method:

    n_par_var = length(model_var_list) # number of extra parameters
    n_par_tot_isov =  [n_par_tot_isov_base]
    label_tot_isov = Vector{Vector{String}}(undef, 0)
    baseStr = SimpleBase ? "baseSimp" : "base"
    push!(label_tot_isov, [baseStr])
    if RESCAL
        push!(n_par_tot_isov,n_par_tot_isov_base)
        push!(label_tot_isov, [baseStr*"Resc"])
    end

    f_tot_len  = minimum([nens-(n_par_tot_isov[1]+mdof)+1,n_par_var+1])
    f_tot_isov = Vector{Vector{Function}}(undef, f_tot_len)
    f_tot_isov[1] = !RESCAL ? Vector{Function}(undef, 1) : Vector{Function}(undef, 2)
    f_tot_isov[1][1] = (x,p) -> isov_basemodel(x,p)
    if RESCAL
        f_tot_isov[1][2] = (x,p) -> isov_basemodel_resc(x,p)
    end

    for n = 2:minimum([nens-(n_par_tot_isov[1]+mdof)+1,n_par_var+1])
        aux = filter(x->sum(x)==n-1, model_map)
        f_tot_isov[n] = !RESCAL ? Vector{Function}(undef, length(aux)) : Vector{Function}(undef, 2*length(aux))

        for (k, a) in enumerate(aux)
            nphi2_term = sum([elem in ["phi2inv","phi2sqr","phi2log","logphi2"] for elem in model_var_label[a]])
            nphi4_term = sum([elem in ["phi4inv","phi4sqr","phi4log","logphi4"] for elem in model_var_label[a]])
            na_term    = sum([elem in ["a2loga","a3","a4"] for elem in model_var_label[a]])  # a2phi2  a2phi4   a3phi2
            if (("a2loga" ∈ model_var_label[a] || "a4" ∈ model_var_label[a] || "a3phi2" ∈ model_var_label[a]) && "a3" ∉ model_var_label[a]) || (nphi2_term > nphi2max) || (nphi4_term > nphi4max) || (na_term > namax)
                continue
            end
            push!(n_par_tot_isov, n_par_tot_isov[1]+n-1)
            push!(label_tot_isov, union([baseStr],model_var_label[a]))
            if !RESCAL
                f_tot_isov[n][k] = (x,p) -> isov_basemodel(x,p) .+ sum([p[i+n_par_tot_isov[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a]))
            else
                push!(n_par_tot_isov, n_par_tot_isov[1]+n-1)
                push!(label_tot_isov, union([baseStr*"Resc"],model_var_label[a]))
                f_tot_isov[n][2*k-1] = (x,p) -> isov_basemodel(x,p)      .+ sum([p[i+n_par_tot_isov[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a]))
                f_tot_isov[n][2*k]   = (x,p) -> isov_basemodel_resc(x,p) .+ sum([p[i+n_par_tot_isov[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a]))
            end
        end
        f_tot_isov[n] = [f_tot_isov[n][l] for l in eachindex(f_tot_isov[n]) if isassigned(f_tot_isov[n], l)]
    end

    f_tot_isov = vcat(f_tot_isov...)

    if !isnothing(MultFunc)
        f_tot_isov = [(x,p) -> MultFunc(x) .* f(x,p) for f in f_tot_isov]
    end

    return [f_tot_isov, n_par_tot_isov, label_tot_isov]
end
