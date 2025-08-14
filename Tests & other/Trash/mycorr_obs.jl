using ADerrors
using HVPobs

include("data_management.jl")

julia_script_directory = @__DIR__
 
path_HVP  = joinpath(julia_script_directory, "..", "LatticeData", "HVP_data")
path_rw   = joinpath(julia_script_directory, "..", "LatticeData", "rwf_deflated")

##

function myapply_rw(data::Array{Float64}, W::Matrix{Float64}, idm::Union{Nothing, Vector{Int64}}=nothing)
    nc =  isnothing(idm) ? collect(1:size(data, 1)) : idm
    W1 = W[1, nc]
    W2 = W[2, nc]

    data_r = data .* W1 .* W2
    return (data_r, W1 .* W2)
end

function myapply_rw(data::Array{Float64}, W::Vector{Matrix{Float64}}, cdidm::Vector{Int64}, rep_len::Vector{Int64}; mask::Vector{Bool}=ones(length(rep_len)))

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

function mycorr_obs(cd::CData; real::Bool=true, rw::Union{Array{Float64,2}, Vector{Array{Float64,2}}, Nothing}=nothing, L::Int64=1, nms::Union{Int64, Nothing}=nothing)
    
    real ? data = cd.re_data ./ L^3 : data = cd.im_data ./ L^3
    tvals   = size(data, 2)
    replen  = collect(values(cd.rep_len))
    reptot  = collect(values(cd.replicatot))
    vcfg    = [cd.idm[1+sum(replen[1:k-1]):sum(replen[1:k])] for k in eachindex(replen)]
    replica = Int64.(maximum.(vcfg))
    nms     = isnothing(nms) ?  sum(replica) : nms

    idm = cd.idm[:] 

    mask = [elem in keys(cd.rep_len) for elem in keys(cd.replicatot)]
    idm_sum = [fill(sum(reptot[1:k-1]), (replen.*mask)[k]) for k in eachindex(reptot)]
    idm .+= vcat(idm_sum...)

    if isnothing(rw)
        obs = [uwreal(data[:,t], cd.id, collect(values(cd.replicatot)), idm, cd.nms) for t in 1:tvals]
    else
        if length(replen)  == length(reptot) == 1 
            data_r, W = myapply_rw(data, rw, cd.idm)
        else
            data_r, W = myapply_rw(data, rw, cd.idm, replen,  mask)
        end
        
        ow = [uwreal(data_r[:,t], cd.id, collect(values(cd.replicatot)), idm, cd.nms) for t in 1:tvals]
        W_obs = uwreal(W, cd.id, collect(values(cd.replicatot)), idm, cd.nms)
        obs = [ow[t] / W_obs for t in 1:tvals]
    end

    return Corr(obs, cd)
end

function myget_corr(path::String, ens::EnsInfo, fl::String, g::String; path_rw::Union{String, Nothing}=nothing, L::Int64=1, frw_bcwd::Bool=false)

    cdata = get_data(path, ens.id, fl, g)
    rw = isnothing(path_rw) ? nothing : get_rw(path_rw, ens.id)
    corr = mycorr_obs(cdata, real=true, rw=rw, L=L)
    if frw_bcwd
        frwd_bckwrd_symm!(corr)
    end

    return corr 
end

function mycorr33(path_data::String, ens::EnsInfo; sector::String="light", path_rw::Union{Nothing,String}=nothing, L::Int64=1, frw_bcwd::Bool=true, impr::Bool=true, impr_set::String="1", cons::Bool=true, std::Bool=false)
    
    Gamma_l = ["V1V1", "V2V2", "V3V3", "V1T10", "V2T20", "V3T30"]
    Gamma_c = ["V1V1c", "V2V2c", "V3V3c", "V1cT10", "V2cT20", "V3cT30"]

    v1v1 = myget_corr(path_data, ens, sector, Gamma_l[1], path_rw=path_rw, frw_bcwd=false, L=L)
    v2v2 = myget_corr(path_data, ens, sector, Gamma_l[2], path_rw=path_rw, frw_bcwd=false, L=L)
    v3v3 = myget_corr(path_data, ens, sector, Gamma_l[3], path_rw=path_rw, frw_bcwd=false, L=L)

    if impr
        v1t10 = myget_corr(path_data, ens, sector, Gamma_l[4], path_rw=path_rw, frw_bcwd=false, L=L)
        v2t20 = myget_corr(path_data, ens, sector, Gamma_l[5], path_rw=path_rw, frw_bcwd=false, L=L)
        v3t30 = myget_corr(path_data, ens, sector, Gamma_l[6], path_rw=path_rw, frw_bcwd=false, L=L)

        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta)
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
        end
        improve_corr_vkvk!(v1v1, v1t10, 2*cv_l, std=std)
        improve_corr_vkvk!(v2v2, v2t20, 2*cv_l, std=std)
        improve_corr_vkvk!(v3v3, v3t30, 2*cv_l, std=std)
    end
    g33 = Corr(0.5 .* (v1v1.obs .+ v2v2.obs .+ v3v3.obs)./3, v1v1.id, "G33ll")
    if frw_bcwd
        frwd_bckwrd_symm!(g33)
    end

    if cons
        v1v1_c = myget_corr(path_data, ens, sector, Gamma_c[1], path_rw=path_rw, frw_bcwd=false, L=L)
        v2v2_c = myget_corr(path_data, ens, sector, Gamma_c[2], path_rw=path_rw, frw_bcwd=false, L=L)
        v3v3_c = myget_corr(path_data, ens, sector, Gamma_c[3], path_rw=path_rw, frw_bcwd=false, L=L)

        if impr
            v1t10_c = myget_corr(path_data, ens, sector, Gamma_c[4], path_rw=path_rw, frw_bcwd=false, L=L)
            v2t20_c = myget_corr(path_data, ens, sector, Gamma_c[5], path_rw=path_rw, frw_bcwd=false, L=L)    
            v3t30_c = myget_corr(path_data, ens, sector, Gamma_c[6], path_rw=path_rw, frw_bcwd=false, L=L)
            
            beta = ens.beta
            if impr_set == "1"
                cv_l = cv_loc(beta) 
                cv_c = cv_cons(beta)
            elseif impr_set == "2"
                cv_l = cv_loc_set2(beta)
                cv_c = cv_cons_set2(beta)
            end
            improve_corr_vkvk_cons!(v1v1_c, v1t10, v1t10_c, cv_l, cv_c, std=std)
            improve_corr_vkvk_cons!(v2v2_c, v2t20, v2t20_c, cv_l, cv_c, std=std)
            improve_corr_vkvk_cons!(v3v3_c, v3t30, v3t30_c, cv_l, cv_c, std=std)
        end
        g33_c = Corr(0.5 .* (v1v1_c.obs .+ v2v2_c.obs .+ v3v3_c.obs)./3, v1v1_c.id, "G33lc")
        if frw_bcwd
            frwd_bckwrd_symm!(g33_c)
        end
    end

    !cons ? (return g33) : (return g33, g33_c)
end

##

cd = get_data(path_HVP, "N101", "charm", "V1V1")
real = true
rw = get_rw(path_rw, "N101")
L = 1
nms = nothing

cd.rep_len
cd.replicatot

##

real ? data = cd.re_data ./ L^3 : data = cd.im_data ./ L^3
tvals   = size(data, 2)
replen  = collect(values(cd.rep_len))
reptot  = collect(values(cd.replicatot))
vcfg    = [cd.idm[1+sum(replen[1:k-1]):sum(replen[1:k])] for k in eachindex(replen)]
replica = Int64.(maximum.(vcfg))
nms     = isnothing(nms) ?  sum(replica) : nms

idm = cd.idm[:] 

mask = [elem in keys(cd.rep_len) for elem in keys(cd.replicatot)]
idm_sum = [fill(sum(reptot[1:k-1]), (replen.*mask)[k]) for k in eachindex(reptot)]
idm .+= vcat(idm_sum...)

if isnothing(rw)
    obs = [uwreal(data[:,t], cd.id, collect(values(cd.replicatot)), idm, cd.nms) for t in 1:tvals]
else
    if length(replen)  == length(reptot) == 1 
        data_r, W = myapply_rw(data, rw, cd.idm)
    else
        data_r, W = myapply_rw(data, rw, cd.idm, replen,  mask)
    end
    
    ow = [uwreal(data_r[:,t], cd.id, collect(values(cd.replicatot)), idm, cd.nms) for t in 1:tvals]
    W_obs = uwreal(W, cd.id, collect(values(cd.replicatot)), idm, cd.nms)
    obs = [ow[t] / W_obs for t in 1:tvals]
end

resultV1V1 = Corr(obs, cd)
-0.5*resultV1V1.obs
-0.5*(resultV1V1.obs[1:6]+resultV2V2.obs[1:6]+resultV3V3.obs[1:6])/3

##

IMPR      = true
IMPR_SET  = "1" 
RENORM    = true
STD_DERIV = false

Z3 = get_Z3(EnsInfo("H101"), impr_set=IMPR_SET) 

G33ll, G33lc = corr33(path_HVP, EnsInfo("H101"), sector="light", path_rw=path_rw, L=1, frw_bcwd=true, impr=true, impr_set=IMPR_SET, cons=true, std=true)

uwerr.(G33ll.obs); uwerr.(G33lc.obs)
G33ll.obs[2:6]
G33lc.obs[2:6]

if RENORM
    renormalize!(G33ll, Z3^2)
    renormalize!(G33lc, Z3)
end


G33ll.obs[1:6]



corr33(path_HVP, EnsInfo("H101"), sector="light", path_rw=path_rw, L=1, frw_bcwd=true, impr=true, impr_set=IMPR_SET, cons=true, std=true)

Ztest = Zvc_l["J500"]

uwerr(Ztest)
details(Ztest)


##

D452corr33 = corr33(path_HVP, EnsInfo("D452"), sector="light", path_rw=path_rw, L=1, frw_bcwd=true, impr=true, impr_set=IMPR_SET, cons=true, std=true)
uwerr(D452corr33[1].obs[10])

details(D452corr33[1].obs[10])