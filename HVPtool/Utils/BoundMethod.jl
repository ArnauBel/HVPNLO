
function Eeff(tstar::Int64, obs::Vector{uwreal})

    E_eff = 0.5 * log((obs[tstar+1] / obs[tstar+2]) ^2)
    return E_eff
end

# function corr_bound(t::Vector{Int64}, tcut::Int64, obs::Vector{uwreal}, ens::EnsInfo, Eeff::Union{uwreal,Float64}; allowPBC::Bool=false)
#     T = 2*t[end]

#     Gcut = obs[tcut]
#     GPBC(x0) = exp(-Eeff*x0) + exp(-Eeff*(T-x0))
#     G(x0) = exp(-Eeff*x0)

#     if split(ens.id,"")[end-1]=="5" && allowPBC
#         UBarray = (Gcut/GPBC(tcut)) .*  GPBC.(t[tcut:end])
#     else
#         UBarray = Gcut/G(tcut) .*  G.(t[tcut:end])
#     end
#     return UBarray[2:end]
# end
# corr_bound(t::Vector{Int64}, tcut::Int64, corr::Corr, ens::EnsInfo, Eeff::Union{uwreal,Float64}; allowPBC::Bool=false) = corr_bound(t, tcut, corr.obs, ens, Eeff; allowPBC=allowPBC)

function corr_bound(t::Vector{Int64}, tcut::Int64, obs::Vector{uwreal}, Eeff::Union{uwreal,Float64})
    Gcut = obs[tcut]
    G(x0) = exp(-Eeff*x0)
    UBarray = Gcut/G(tcut) .*  G.(t[tcut:end])
    return UBarray[2:end]
end
corr_bound(t::Vector{Int64}, tcut::Int64, corr::Corr, Eeff::Union{uwreal,Float64}) = corr_bound(t, tcut, corr.obs, Eeff)


# function buonding_method(ub::Vector{uwreal},lb::Vector{uwreal};PLAT::Bool=false,AVER::Bool=false)
#     averb = (ub.+lb)./2; uwerr.(averb)
#     x0 = findfirst(abs.(value.(ub).-value.(lb)) .< 0.5.*err.(averb))
#     ∆x = findfirst(abs.(averb[x0:end].-averb[x0]) .> err(averb[x0]))
#     if isnothing(∆x) || ∆x > 5
#         hvp  = averb[x0]
#         syst = 0.0
#         xend = isnothing(∆x) ? length(averb) : (x0-1) + ∆x
#     else
#         xend=x0+3
#         x0 > 3 ? (x0-=3) : (x0=1)
#         hvp  = sum(averb[x0:xend])/length(averb[x0:xend])
#         aux1 = sum(averb[x0:xend].^2)/length(averb[x0:xend])
#         aux2 = hvp^2
#         syst = sqrt(abs(value(aux1 - aux2)))
#     end
#     # plateau_fm = value(aEns).*(collect(x0:xend).+tcut0.-2)

#     returnVec = [hvp,syst]
#     PLAT ? push!(returnVec,[x0,xend]) : nothing
#     AVER ? push!(returnVec,averb) : nothing
#     return returnVec
# end

function bounding_method(ub::Vector{uwreal},lb::Vector{uwreal},aEns::uwreal;PLAT::Bool=false,AVER::Bool=false,tcut0::Union{Int64,Nothing}=nothing,tstep::Union{Int64,Nothing}=nothing)
    if isnothing(tcut0)
        tcut0 = 1
    end
    if isnothing(tstep)
        tstep = 1
    end
    averb = (ub.+lb)./2; uwerr.(averb)
    x0    = findfirst(abs.(value.(ub).-value.(lb)) .< err.(averb))
    # println("x0   = $(x0)")
    xend  = (4*tstep*aEns.mean < 0.25) ? Int64(x0 + round(0.25/(tstep*aEns.mean),RoundUp)) : Int64(x0 + round(4/tstep,RoundUp))
    # println("xend = $(xend)")
    if xend > length(averb)
        xend = length(averb)
        @warn(" - BM converging too slow! t=tmax (T/2) was reached. Average will be produced with shorter plateau.")
    end
    hvp  = mean(averb[x0:xend])
    aux1 = mean(averb[x0:xend].^2)
    aux2 = hvp^2
    syst = sqrt(abs(value(aux1 - aux2)))

    if PLAT
        !isnothing(tcut0) ? plateau_fm = aEns.mean.*(tstep.*collect(x0:xend) .+ (tcut0-1-tstep)) : error("tcut0 required to output plateau")
    end

    returnVec = [hvp,syst]
    PLAT ? push!(returnVec,plateau_fm) : nothing
    AVER ? push!(returnVec,averb) : nothing
    return returnVec
end
