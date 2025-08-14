


##-----------------------------------------------------------------------

function apply_rw(data::Array{Float64}, W::Matrix{Float64}, cdidm::Union{Nothing, Vector{Int64}}=nothing)
    nc =  isnothing(cdidm) ? collect(1:size(data, 1)) : cdidm
    W1 = W[1, nc]
    W2 = W[2, nc]

    data_r = data .* W1 .* W2
    return (data_r, W1 .* W2)
end

function apply_rw(data::Array{Float64}, W::Vector{Matrix{Float64}}, cdidm::Vector{Int64}, rep_len::Vector{Int64}, mask::Vector{Bool})

    chunk(arr, n::Vector{Int64}) = [arr[1+ sum(n[1:i-1]):sum(n[1:i])] for i in eachindex(n)]
    idm = chunk(cdidm, rep_len)
    rw1 = []
    rw2 = []

    W = [W[i] for (i, mask_val) in enumerate(mask) if mask_val]

    rw1 = [W[k][1, idm[k]] for k in eachindex(idm)]
    rw2 = [W[k][2, idm[k]] for k in eachindex(idm)]

    rw = vcat([rw1[k] .* rw2[k] for k in eachindex(idm)]...)
    data_r = data .* rw 
    return (data_r, rw)
end


##-----------------------------------------------------------------------


ens = "H101"
corr_light = get_corr(path_HVP, EnsInfo(ens), "light", "V1V1",frw_bcwd=false, path_rw=path_rw)
corr_charm = get_corr(path_HVP, EnsInfo(ens), "charm", "V1V1",frw_bcwd=false, path_rw=path_rw)
corr_disc = get_corr_disc(path_HVP, EnsInfo("N101"), "88", frw_bcwd=false, path_rw=path_rw)


corr_disc(path_HVP, EnsInfo("N101"), path_rw = path_rw, frw_bcwd = true, impr = IMPR, impr_set = IMPR_SET, discr = ["ll","lc","cc"], std = STD_DERIV)

##

ens = "N101"
# cd = get_data_disc(path_HVP, ens, "cc")["VV"]
cd = get_data(path_HVP, "N101", "charm", "V1V1")
KEYS = collect(keys(cdata))
rw = isnothing(path_rw) ? nothing : get_rw(path_rw, ens)

#---------------------------------------------------------

reporder  = [parse(Int64,collect(keys(cd.rep_len))[i][end]) for i in 1:length(cd.rep_len)]
myreplicatot = OrderedDict{String, Int64}()
myrw = Vector{Matrix{Float64}}()
if !all([reporder[i+1]>reporder[i] ? true : false for i in 1:(length(reporder)-1)])
   for (i,rep) in enumerate(collect(keys(cd.rep_len)))
       myreplicatot[rep] = cd.replicatot[rep]
       push!(myrw,rw[findall(x -> x == rep, collect(keys(cd.replicatot)))][1])
   end
else
   myreplicatot = cd.replicatot
   myrw = rw
end

cd.replicatot
myreplicatot
cd.rep_len

rw
myrw

data = cd.re_data
tvals   = size(data, 2)
replen  = collect(values(cd.rep_len))
reptot  = collect(values(myreplicatot))
vcfg    = [cd.idm[1+sum(replen[1:k-1]):sum(replen[1:k])] for k in eachindex(replen)]
replica = Int64.(maximum.(vcfg))
nms     = sum(replica)

idm = cd.idm[:] 
mask = [elem in keys(cd.rep_len) for elem in keys(myreplicatot)]
maskedrep = [haskey(cd.rep_len, key) ? cd.rep_len[key] : 0 for key in keys(myreplicatot)]
idm_sum = [fill(sum(reptot[1:k-1]), (maskedrep)[k]) for k in eachindex(reptot)]
idm .+= vcat(idm_sum...)

if isnothing(rw)
    obs = [uwreal(data[:,t], cd.id, collect(values(cd.replicatot)), idm, cd.nms) for t in 1:tvals]
else
    if length(replen)  == length(reptot) == 1 
        data_r, W = apply_rw(data, myrw, cd.idm)
    else
        data_r, W = apply_rw(data, myrw, cd.idm, replen, mask)
    end
    
    ow = [uwreal(data_r[:,t], cd.id, collect(values(cd.replicatot)), idm, cd.nms) for t in 1:tvals]
    W_obs = uwreal(W, cd.id, collect(values(cd.replicatot)), idm, cd.nms)
    obs = [ow[t] / W_obs for t in 1:tvals]
end

Corr(obs, cd)

##-----------------------------------------------------------------------


function apply_myrw(data::Array{Float64}, W::Vector{Matrix{Float64}}, cdidm::Vector{Int64}, rep_len::Vector{Int64})

    chunk(arr, n::Vector{Int64}) = [arr[1+ sum(n[1:i-1]):sum(n[1:i])] for i in eachindex(n)]
    idm = chunk(cdidm, rep_len)
    rw1 = []
    rw2 = []

    rw1 = [W[k][1, idm[k]] for k in eachindex(idm)]
    rw2 = [W[k][2, idm[k]] for k in eachindex(idm)]

    rw = vcat([rw1[k] .* rw2[k] for k in eachindex(idm)]...)
    data_r = data .* rw 
    return (data_r, rw)
end




##-----------------------------------------------------------------------

ens = "N101"
mycorr_light = get_corr(path_HVP, EnsInfo(ens), "light", "V1V1",frw_bcwd=false, path_rw=path_rw)
mycorr_charm = get_corr(path_HVP, EnsInfo(ens), "charm", "V1V1",frw_bcwd=false, path_rw=path_rw)
mycorr_disc = get_corr_disc(path_HVP, EnsInfo("N101"), "88", frw_bcwd=false, path_rw=path_rw)


cd = get_data(path_HVP, ens, "charm", "V1V1")
cd = get_data_disc(path_HVP, ens, "88")["VV"]

rw = isnothing(path_rw) ? nothing : get_rw(path_rw, ens)

cd.replicatot
cd.rep_len

#---------------------------------------------------------

# The idm values for uwreal() are creating by "chunking" and then sorting the cd.idm values
replen   = collect(values(cd.rep_len))
reptot   = collect(values(cd.replicatot))

reporder  = [parse(Int64,collect(keys(cd.rep_len))[i][end]) for i in 1:length(cd.rep_len)]
if !all([reporder[i+1]>reporder[i] ? true : false for i in 1:(length(reporder)-1)])            # check if order is preserved when reading the data (my be a problem for disconnected)
    chunk(arr, n::Vector{Int64}) = [arr[1+ sum(n[1:i-1]):sum(n[1:i])] for i in eachindex(n)]
    vec_idm = chunk(cd.idm, replen)
    vec_idm_sorted = Vector{Vector{Int64}}()
    for rep in collect(keys(cd.replicatot))
        push!(vec_idm_sorted,vec_idm[findall(x -> x == rep, collect(keys(cd.rep_len)))][1])
    end
    idm = vcat(vec_idm_sorted...)
else
    idm = cd.idm[:] 
end

maskedrep = [haskey(cd.rep_len, key) ? cd.rep_len[key] : 0 for key in keys(cd.replicatot)]
idm_sum = [fill(sum(reptot[1:k-1]), (maskedrep)[k]) for k in eachindex(reptot)]
idm .+= vcat(idm_sum...)

# The reweighting is obtained by defining a new reweighting Vector{matrix} which follows the order and existance of the data
real ? data = cd.re_data ./ L^3 : data = cd.im_data ./ L^3
tvals    = size(data, 2)

#myreplicatot = OrderedDict{String, Int64}()
myrw = Vector{Matrix{Float64}}()
for rep in collect(keys(cd.rep_len))
    #myreplicatot[rep] = cd.replicatot[rep]
    push!(myrw,rw[findall(x -> x == rep, collect(keys(cd.replicatot)))][1])
end

#myreptot = collect(values(myreplicatot))

if isnothing(rw)
    obs = [uwreal(data[:,t], cd.id, collect(values(cd.replicatot)), idm, cd.nms) for t in 1:tvals]
else
    data_r, W = apply_myrw(data, myrw, cd.idm, replen)
    
    ow = [uwreal(data_r[:,t], cd.id, collect(values(cd.replicatot)), idm, cd.nms) for t in 1:tvals]
    W_obs = uwreal(W, cd.id, collect(values(cd.replicatot)), idm, cd.nms)
    obs = [ow[t] / W_obs for t in 1:tvals]
end

return Corr(obs, cd)


##

uwerr.(mycorr.obs)
uwerr.(mycorr_light.obs)
uwerr.(mycorr_charm.obs)
uwerr.(mycorr_disc["VV"].obs)

mycorr.obs[1:4]
mycorr_light.obs[1:4]
mycorr_charm.obs[1:4]
mycorr_disc["VV"].obs[1:4]

mchist(mycorr.obs[10], "N101")#[275:284]
mchist(mycorr_disc["VV"].obs[10], "N101")#[275:284]
findall(x -> x != 0.0, mchist(mycorr.obs[10], "N101"))
findall(x -> x != 0.0, mchist(mycorr_disc["VV"].obs[10], "N101"))