
using IterTools
using ADerrors
import ADerrors: err

# include("uwConst.jl")

# Permutated terms:

a2(x)       = x[:,1]                                                              # (a^2/8t0)
a2loga(x)   = x[:,1] .* log.(x[:,1])                                              # (a^2/8t0)*log(a^2/8t0)
a3(x)       = x[:,1].^(3/2)                                                       # (a^2/8t0)^{3/2}
a4(x)       = x[:,1].^(2)                                                         # (a^2/8t0)^{2}


phi2(x)     = (x[:,2] .- value.(phi2_ph))
phi2sqr(x)  = (x[:,2].^2 .- value.(phi2_ph).^2)                                   # (ϕ2^2 - (ϕ2^{ph})^2) 
phi2inv(x)  = (x[:,2].^(-1) .- value.(phi2_ph).^(-1))                             # (1/ϕ2 - 1/(ϕ2^{ph}))
phi2log(x)  = x[:,2] .* log.(x[:,2]) .- value.(phi2_ph) .* log.(value.(phi2_ph))  # (ϕ2log(ϕ2) - ϕ2^{ph}log(ϕ2^{ph}))
logphi2(x)  = log.(x[:,2]) .- log.(value.(phi2_ph))                               # (log(ϕ2) - log(ϕ2^{ph}))

y(x)        = (x[:,2] .- value.(y_ph))
ysqr(x)     = (x[:,2].^2 .- value.(y_ph).^2)                                      # (y^2 - (y^{ph})^2) 
yinv(x)     = (x[:,2].^(-1) .- value.(y_ph).^(-1))                                # (1/y - 1/(y^{ph}))
ylog(x)     = x[:,2] .* log.(x[:,2]) .- value.(y_ph) .* log.(value.(y_ph))        # (ylog(y) - y^{ph}log(y^{ph}))
logy(x)     = log.(x[:,2]) .- log.(value.(y_ph))                                  # (log(y) - log(y^{ph}))


phi4(x)     = (x[:,3] .- value.(phi4_ph))
phi4sqr(x)  = (x[:,3].^2 .- value.(phi4_ph).^2)                                   # (ϕ4^2 - (ϕ4^{ph})^2) 
phi4inv(x)  = (x[:,3].^(-1) .- value.(phi4_ph).^(-1))                             # (1/ϕ4 - 1/(ϕ4^{ph}))
phi4log(x)  = x[:,3] .* log.(x[:,3]) .- value.(phi4_ph) .* log.(value.(phi4_ph))  # (ϕ4log(ϕ4) - ϕ4^{ph}log(ϕ4^{ph}))
logphi4(x)  = log.(x[:,3]) .- log.(value.(phi4_ph))                               # (log(ϕ4) - log(ϕ4^{ph}))

z(x)        = (x[:,3] .- value.(z_ph))
zsqr(x)     = (x[:,3].^2 .- value.(z_ph).^2)                                      # (z^2 - (z^{ph})^2) 
zinv(x)     = (x[:,3].^(-1) .- value.(z_ph).^(-1))                                # (1/z - 1/(z^{ph}))
zlog(x)     = x[:,3] .* log.(x[:,3]) .- value.(z_ph) .* log.(value.(z_ph))        # (zlog(z) - z^{ph}log(z^{ph}))
logz(x)     = log.(x[:,3]) .- log.(value.(z_ph))                                  # (log(z) - log(z^{ph}))


a2phi2(x)   = x[:,1] .* (x[:,2] .- value.(phi2_ph))                               # (a^2/8t0)*(ϕ2 - ϕ2^{ph})     
a2phi4(x)   = x[:,1] .* (x[:,3] .- value.(phi4_ph))                               # (a^2/8t0)*(ϕ4 - ϕ4^{ph})

a2y(x)      = x[:,1] .* (x[:,2] .- value.(y_ph))                                  # (a^2/8t0)*(y - y^{ph})     
a2z(x)      = x[:,1] .* (x[:,3] .- value.(z_ph))                                  # (a^2/8t0)*(z - z^{ph})


a3phi2(x)   = x[:,1].^(3/2) .* (x[:,2] .- value.(phi2_ph))                        # (a^2/8t0)^{3/2}*(ϕ2 - ϕ2^{ph})   

a3y(x)      = x[:,1].^(3/2) .* (x[:,2] .- value.(y_ph))                           # (a^2/8t0)^{3/2}*(ϕ2 - ϕ2^{ph})   

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
    fPiresc::Bool=false,
    a2resc::Bool=true,
    nmPi_max::Union{Int64,Vector{Int64}}=1,
    nmK_max::Int64=1,
    na_max::Int64=3
    )

    nens = length(EnsList)

    model_var_label = string.(model_var_list)
    model_map = generate_model_map(length(model_var_list))

    if !fPiresc
        if any(x -> x in ["a2y","a2z","a3y","y","ysqr","ylog","yinv","logy","z","zsqr","zlog","zinv","logz"], model_var_label)
            error("Found term in model list not allowed without rescaling")
        end
    else
        if any(x -> x in ["a2phi2","a2phi4","a3phi2","phi2","phi2sqr","phi2log","phi2inv","logphi2","phi4","phi4sqr","phi4log","phi4inv","logphi4"], model_var_label)
            error("Found term in model list not allowed when rescaling")
        end
    end

    if typeof(nmPi_max) == Int64
        nmPi_max = [nmPi_max,Inf]
    end
    
    isov_basemodel_t0(x,p) = SimpleBase ? (p[1] .+ p[2] .* a2(x) .+ p[3] .* phi2(x)) : (p[1] .+ p[2] .* a2(x) .+ p[3] .* phi2(x) .+ p[4] .* phi4(x))
    isov_basemodel_fπ(x,p) = SimpleBase ? (p[1] .+ p[2] .* a2(x) .+ p[3] .* y(x))    : (p[1] .+ p[2] .* a2(x) .+ p[3] .* y(x) .+ p[4] .* z(x))

    isov_basemodel = !fPiresc ? isov_basemodel_t0 : isov_basemodel_fπ

    if a2resc
        beta_list = getfield.(EnsList,:beta)
        isov_basemodel_a2resc_t0(x,p) = SimpleBase ? (p[1] .+ p[2] .* a2_rescaling.(beta_list) .* a2(x) .+ p[3] .* phi2(x)) : (p[1] .+ p[2] .* a2_rescaling.(beta_list) .* a2(x) .+ p[3] .* phi2(x) .+ p[4] .* phi4(x))
        isov_basemodel_a2resc_fπ(x,p) = SimpleBase ? (p[1] .+ p[2] .* a2_rescaling.(beta_list) .* a2(x) .+ p[3] .* y(x))    : (p[1] .+ p[2] .* a2_rescaling.(beta_list) .* a2(x) .+ p[3] .* y(x) .+ p[4] .* z(x))

        isov_basemodel_a2resc = !fPiresc ? isov_basemodel_a2resc_t0 : isov_basemodel_a2resc_fπ
    end
    n_par_tot_isov_base = SimpleBase ? 3 : 4

    # Method:

    n_par_var = length(model_var_list) # number of extra parameters
    n_par_tot_isov =  [n_par_tot_isov_base]
    label_tot_isov = Vector{Vector{String}}(undef, 0)
    baseStr = SimpleBase ? "baseSimp" : "base"
    push!(label_tot_isov, [baseStr])
    if a2resc
        push!(n_par_tot_isov,n_par_tot_isov_base)
        push!(label_tot_isov, [baseStr*"Resc"])
    end

    f_tot_len  = minimum([nens-(n_par_tot_isov[1]+mdof)+1,n_par_var+1])
    f_tot_isov = Vector{Vector{Function}}(undef, f_tot_len)
    f_tot_isov[1] = !a2resc ? Vector{Function}(undef, 1) : Vector{Function}(undef, 2)
    f_tot_isov[1][1] = (x,p) -> isov_basemodel(x,p)
    if a2resc
        f_tot_isov[1][2] = (x,p) -> isov_basemodel_a2resc(x,p)
    end

    for n = 2:minimum([nens-(n_par_tot_isov[1]+mdof)+1,n_par_var+1])
        aux = filter(x->sum(x)==n-1, model_map)
        f_tot_isov[n] = !a2resc ? Vector{Function}(undef, length(aux)) : Vector{Function}(undef, 2*length(aux))

        for (k, a) in enumerate(aux)
            if !fPiresc
                nphi2_term     = sum([elem in ["phi2inv","phi2sqr","phi2log","logphi2","a2phi2","a3phi2"] for elem in model_var_label[a]])
                nphi4_term     = sum([elem in ["phi4inv","phi4sqr","phi4log","logphi4","a2phi4"] for elem in model_var_label[a]])
                na_term        = sum([elem in ["a2loga","a3","a4","a2phi2","a2phi4","a3phi2"] for elem in model_var_label[a]])
                nphi2pure_term = sum([elem in ["phi2inv","phi2sqr","phi2log","logphi2"] for elem in model_var_label[a]])

                if (any([term ∈ model_var_label[a] for term in ["a4","a3phi2"]]) && "a3" ∉ model_var_label[a]) || ("a3phi2" ∈ model_var_label[a] && "a2phi2" ∉ model_var_label[a]) || (any([term ∈ model_var_label[a] for term in ["a2phi4","phi4inv","phi4sqr","phi4log","logphi4","a2phi4"]]) && SimpleBase && "phi4" ∉ model_var_label[a]) || (nphi2_term > nmPi_max[1]) || (nphi2pure_term > nmPi_max[2]) || (nphi4_term > nmK_max) || (na_term > na_max)
                    continue
                end
            else
                ny_term     = sum([elem in ["yinv","ysqr","ylog","logy","a2y","a3y"] for elem in model_var_label[a]])
                nz_term     = sum([elem in ["zinv","zsqr","zlog","logz","a2z"] for elem in model_var_label[a]])
                na_term     = sum([elem in ["a2loga","a3","a4","a2y","a2z","a3y"] for elem in model_var_label[a]])
                nypure_term = sum([elem in ["yinv","ysqr","ylog","logy"] for elem in model_var_label[a]])
                
                if (any([term ∈ model_var_label[a] for term in ["a4","a3y"]]) && "a3" ∉ model_var_label[a]) || ("a3y" ∈ model_var_label[a] && "a2y" ∉ model_var_label[a]) || (any([term ∈ model_var_label[a] for term in ["a2z","zinv","zsqr","zlog","logz","a2z"]]) && SimpleBase && "z" ∉ model_var_label[a]) || (ny_term > nmPi_max[1]) || (nypure_term > nmPi_max[2]) || (nz_term > nmK_max) || (na_term > na_max)
                    continue
                end
            end
            push!(n_par_tot_isov, n_par_tot_isov[1]+n-1)
            push!(label_tot_isov, union([baseStr],model_var_label[a]))
            if !a2resc
                f_tot_isov[n][k] = (x,p) -> isov_basemodel(x,p) .+ sum([p[i+n_par_tot_isov[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a]))
            else
                push!(n_par_tot_isov, n_par_tot_isov[1]+n-1)
                push!(label_tot_isov, union([baseStr*"Resc"],model_var_label[a]))
                f_tot_isov[n][2*k-1] = (x,p) -> isov_basemodel(x,p)        .+ sum([p[i+n_par_tot_isov[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a]))
                f_tot_isov[n][2*k]   = (x,p) -> isov_basemodel_a2resc(x,p) .+ sum([p[i+n_par_tot_isov[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a]))
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
