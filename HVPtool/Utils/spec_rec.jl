
function findfirst_uninterrupted(bools)
    i = 0
    stop = false
    while !stop
        arg = findfirst(bools[i+1:end])
        i = !isnothing(arg) ? (i+arg) : nothing
        if isnothing(i) || isnothing(findfirst(.!bools[i+1:end]))  # stops if there are unterrupted trues or if it reached the end
            stop = true
        end
    end
    return i
end


function get_spectr_data(path::String,ens::EnsInfo)
    ph5 = joinpath(path,"$(ens.id)-pSq0-T1up.h5")
    file = h5open(ph5, "r")

    # file["summary"]

    # Read the energies
    EnKeys = HDF5.keys(file["pSq0-T1up/spectrum"])
    E = uwreal[]
    for n in EnKeys
        En = read(file["pSq0-T1up/spectrum/$n"])
        push!(E,uwreal(jackknife_err(En),n))
    end

    # Read the overlaps
    ZKeys = HDF5.keys(file["pSq0-T1up/overlaps"])
    Z = uwreal[]
    Z_impr = uwreal[]
    for p in ZKeys
        Zp = read(file["pSq0-T1up/overlaps/$p"])
        Zp_impr = read(file["pSq0-T1up/overlaps_imp/$p"])
        push!(Z,uwreal(jackknife_err(Zp),p))
        push!(Z_impr,uwreal(jackknife_err(Zp_impr),p))
    end

    return E, Z, Z_impr
end

# function corr_n(n::Int64,Union{Vector{Int64},Vector{Float64}},E::Vector{uwreal},Z::Vector{uwreal},L::Int64)
#     return Z[n] / (L^3) .* exp.(-E[n] .* t)
# end
function corr_n(n::Int64,t::Union{Vector{Int64},Vector{Float64}},E::Vector{uwreal},Z::Vector{uwreal},L::Int64)
    return (1/2) * Z[n] / (L^3) .* exp.(-E[n] .* t) # the 1/2 factor convers light-light corr to the iso-vector corr
end

function reconstr_corr(
        ens::EnsInfo,
        E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};
        nmax::Union{Nothing,Int64}=nothing,
        impr_set::String="1",
        IMPR::Bool=true,
        RENORM::Bool=true,
        FreeEN::Bool=false,
        total::Bool=false,
    )
    length(Z) != length(Zimpr) ? error("Length of Z and Zimpr are not the same") : nothing
    # nmax chooses the tower of states to be added
    # if nmax is not defined, then the max number of available states are used
    # FreeEN ensures that the last energy level is not used when saturating the corr; usefull for impr BM
    if isnothing(nmax)
        FreeEN ? nmax = minimum([length(E)-1,length(Z)]) : nmax = minimum([length(E),length(Z)])
    else
        if FreeEN
            nmax > length(E)-1 ? error("Input E data not large enough to tower $nmax (-1) states") : nothing
        else
            nmax > length(E) ? error("Input E data not large enough to tower $nmax states") : nothing
        end
        nmax > length(Z) ? error("Input Z data not large enough to tower $nmax states") : nothing
    end
    
    t = collect(1:Int64(HVPobs.Data.get_T(ens.id))/2+1)
    corrVec_PiPi = [corr_n(1,t.-1,E,Z.^2 ,ens.L)]
    corrVec_Impr = [corr_n(1,t.-1,E,Zimpr.*Z,ens.L)]
    for n=2:nmax
        push!(corrVec_PiPi,corrVec_PiPi[n-1] .+ corr_n(n,t.-1,E,Z.^2 ,ens.L))
        push!(corrVec_Impr,corrVec_Impr[n-1] .+ corr_n(n,t.-1,E,Zimpr.*Z,ens.L))
    end

    # corrVec = deepcopy(corrVec_PiPi)
    if IMPR
        if impr_set == "1"
            cv_l = cv_loc(ens.beta)
        elseif impr_set == "1old"
            cv_l = cv_loc_old(ens.beta)
        elseif impr_set == "2"
            cv_l = cv_loc_set2(ens.beta)
        else
            error("Impr Set $impr_set not recoginsed, please choose between '1', '1old' and '2'.")
        end
        # [1:Int64(HVPobs.Data.get_T(ens.id)/2+1)]
        # [improve_corr_vkvk!(corrVec[n], corrVec_JPi[n], 2*cv_l, std=std, treelevel=true) for n in collect(1:nmax)]
        corrVec = [corrVec_PiPi[n] .+ (2*cv_l).*corrVec_Impr[n] for n=1:nmax]
    else
        corrVec = corrVec_PiPi
    end

    if RENORM
        Z3 = get_Z3(ens, impr_set=impr_set)
        [renormalize!(corr, Z3^2) for corr in corrVec]
        if !total
            [renormalize!(corr, Z3^2) for corr in corrVec_PiPi]
            [renormalize!(corr, Z3^2) for corr in corrVec_Impr]
        end
    end

    !total ? (return corrVec) : (return corrVec, corrVec_PiPi, corrVec_Impr)
end
reconstr_corr(ensid::String,E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};nmax::Union{Nothing,Int64}=nothing,impr_set::String="1",IMPR::Bool=true,RENORM::Bool=true,FreeEN::Bool=false,total::Bool=false,) = reconstr_corr(EnsInfo(ensid),E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};nmax=nmax,impr_set=impr_set,IMPR=IMPR,RENORM=RENORM,FreeEN=FreeEN,total=total,)
