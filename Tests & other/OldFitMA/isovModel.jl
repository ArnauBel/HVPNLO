
include("HVPtools/Const.jl")

# <<<<===============================================>>>> #
# <<<<==========>>>> ISOVECTOR CHANNEL <<<<==========>>>> #
# <<<<===============================================>>>> # 

# Permutated terms:

phi2sqr(x)  = (x[:,2].^2 .- value.(phi2_ph).^2)                                   # (ϕ2^2 - (ϕ2^{ph})^2) 
phi2inv(x)  = (x[:,2].^(-1) .- value.(phi2_ph).^(-1))                             # (1/ϕ2 - 1/(ϕ2^{ph}))                # not included in 'all'
logphi2(x)  = log.(x[:,2]) .- log.(value.(phi2_ph))                               # (log(ϕ2) - log(ϕ2^{ph}))            # not included in 'all'
phi2log(x)  = x[:,2] .* log.(x[:,2]) .- value.(phi2_ph) .* log.(value.(phi2_ph))  # (ϕ2log(ϕ2) - ϕ2^{ph}log(ϕ2^{ph}))

a2phi2(x)   = x[:,1] .* (x[:,2] .- value(phi2_ph))                                # (a^2/8t0)*(ϕ2 - ϕ2^{ph})     
a2phi4(x)   = x[:,1] .* (x[:,3] .- value.(phi4_ph))                               # (a^2/8t0)*(ϕ4 - ϕ4^{ph})

a3phi2(x)   = x[:,1].^(3/2) .* (x[:,2] .- value.(phi2_ph))                        # (a^2/8t0)^{3/2}*(ϕ2 - ϕ2^{ph})   

a2loga(x)   = x[:,1] .* log.(x[:,1])                                              # (a^2/8t0)*log(a^2/8t0)
a3cutoff(x) = x[:,1].^(3/2)                                                       # (a^2/8t0)^{3/2}
a4cutoff(x) = x[:,1].^(2)                                                         # (a^2/8t0)^{2}

# Call all model permutations

function call_models(type_basemodel::String,type_DA::String,mdof::Int64,nens::Int64)

    # Base model and model combinations:

    if type_basemodel == "phi4"
        isov_basemodel(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+ p[4] .* (x[:,3] .- value.(phi4_ph))
        n_par_tot_isov_base = 4
    elseif type_basemodel == "simple"
        isov_basemodel(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph))
        n_par_tot_isov_base = 3
    end


    if type_DA == "Only(phi2inv,logphi2)"
        model_var_list  = [phi2inv, logphi2]
        model_var_label = ["phi2inv", "logphi2"]
        model_map = [Bool.([i,j]) for i=0:1 for j=0:1]

    elseif type_DA == "Only(a3cutoff,phi2inv,logphi2)"
        model_var_list  = [a3cutoff, phi2inv, logphi2]
        model_var_label = ["a3", "phi2inv", "logphi2"]
        model_map = [Bool.([i,j,k]) for i=0:1 for j=0:1 for k=0:1] 

    elseif type_DA == "All-(a4,a2phi4)"
        model_var_list  = [a3cutoff, a2phi2, a3phi2, phi2sqr, phi2log]
        model_var_label = ["a3", "a2phi2", "a3phi2", "phi2sqr", "phi2log"]
        model_map = [Bool.([i,j,k,l,m]) for i=0:1 for j=0:1 for k=0:1 for l=0:1 for m=0:1] 

    elseif type_DA == "All-(a4,a2loga)"
        model_var_list  = [a3cutoff, a2phi2, a3phi2, a2phi4, phi2sqr, phi2log]
        model_var_label = ["a3", "a2phi2", "a3phi2", "a2phi4", "phi2sqr", "phi2log"]
        model_map = [Bool.([i,j,k,l,m,n]) for i=0:1 for j=0:1 for k=0:1 for l=0:1 for m=0:1 for n=0:1] 

    elseif type_DA == "All-(a4)"
        model_var_list  = [a2loga, a3cutoff, a2phi2, a3phi2, a2phi4, phi2sqr, phi2log]
        model_var_label = ["a2loga", "a3", "a2phi2", "a3phi2", "a2phi4", "phi2sqr", "phi2log"]
        model_map = [Bool.([i,j,k,l,m,n,0]) for i=0:1 for j=0:1 for k=0:1 for l=0:1 for m=0:1 for n=0:1 for o=0:1] 

    elseif type_DA == "All"
        model_var_list  = [a2loga, a3cutoff, a4cutoff, a2phi2, a3phi2, a2phi4, phi2sqr, phi2log]
        model_var_label = ["a2loga", "a3", "a4", "a2phi2", "a3phi2", "a2phi4", "phi2sqr", "phi2log"]
        model_map = [Bool.([i,j,k,l,m,n,o,p]) for i=0:1 for j=0:1 for k=0:1 for l=0:1 for m=0:1 for n=0:1 for o=0:1 for p=0:1]
    else
        error("DA type: $typeDA not considered in the 'isovModel.jl' file")
    end

    # Method:

    n_par_var = length(model_var_list) # number of extra parameters
    n_par_tot_isov =  [n_par_tot_isov_base]
    label_tot_isov = Vector{Vector{String}}(undef, 0)
    push!(label_tot_isov, ["base"])

    # f_tot_isov = Vector{Vector{Function}}(undef, n_par_var+1)
    f_tot_isov = Vector{Vector{Function}}(undef, minimum([nens-(n_par_tot_isov[1]+mdof)+1,n_par_var+1]))
    f_tot_isov[1] = Vector{Function}(undef, 1)
    f_tot_isov[1][1] = (x,p) -> isov_basemodel(x,p)

    # for n = 2:n_par_var+1
    for n = 2:minimum([nens-(n_par_tot_isov[1]+mdof)+1,n_par_var+1])
        aux = filter(x->sum(x)==n-1, model_map)
        f_tot_isov[n] = Vector{Function}(undef, length(aux))
        f_aux = []
        for (k, a) in enumerate(aux)
            if "a4" ∈ model_var_label[a] && "a3" ∉ model_var_label[a]
                continue
            end
            push!(n_par_tot_isov, n_par_tot_isov[1]+n-1)
            push!(label_tot_isov, model_var_label[a])
            f_tot_isov[n][k] = (x,p) -> isov_basemodel(x,p) .+ sum([p[i+n_par_tot_isov[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a]))
            # push!(f_aux,  (x,p) -> isov_basemodel(x,p) .+ sum([p[i+n_par_tot_isov[1]] for i=1:(n-1)] .* (fill(x, n-1) .|> model_var_list[a])))
        end
        f_tot_isov[n] = [f_tot_isov[n][l] for l in eachindex(f_tot_isov[n]) if isassigned(f_tot_isov[n], l)]
    end

    f_tot_isov = vcat(f_tot_isov...)

    return [f_tot_isov, n_par_tot_isov, label_tot_isov]
end
