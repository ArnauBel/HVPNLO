# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using HDF5

using Statistics

using TimerOutputs
using Suppressor

using Plots
using PyPlot
using Colors

# include

include("../tools/const.jl")

# Path definition

julia_script_directory = @__DIR__

path_bdio = joinpath(julia_script_directory, "..",  "..", "ObsBDIO")
path_spectro = joinpath(julia_script_directory, "..", "..", "LMEData", "spectroscopy")

# Ensamble choice

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ensList = ["D200", "E250", "J303"]
ensInfo = EnsInfo.(ensList)

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

@info("Ready")


##-------------------------------

function get_spectr_data(path::String,ens::EnsInfo)
    ph5 = joinpath(path,"$(ens.id)-pSq0-T1up.h5")
    file = h5open(ph5, "r")

    file["summary"]

    # Read the energies
    EnKeys = keys(file["pSq0-T1up/spectrum"])
    E = uwreal[]
    for n in EnKeys
        En = read(file["pSq0-T1up/spectrum/$n"])
        push!(E,uwreal([mean(En),std(En)],n))
    end

    # Read the overlaps
    ZKeys = keys(file["pSq0-T1up/overlaps"])
    Z = uwreal[]
    Z_impr = uwreal[]
    for p in ZKeys
        Zp = read(file["pSq0-T1up/overlaps/$p"])
        Zp_impr = read(file["pSq0-T1up/overlaps_imp/$p"])
        push!(Z,uwreal([mean(Zp),std(Zp)],p))
        push!(Z_impr,uwreal([mean(Zp_impr),std(Zp_impr)],p))
    end

    return E, Z, Z_impr
end

function corr_n(n::Int64,t::Vector{Int64},E::Vector{uwreal},Z::Vector{uwreal},L::Int64)
    return Z[n] / (L^3) .* exp.(-E[n] .* t)
    # return Z[n] / (L^3) .* exp.(-E[n] .* t)
end

function reconstr_corr(
        ens::EnsInfo,
        E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};
        nmax::Union{Nothing,Int64}=nothing,
        FreeEN::Bool=false,total::Bool=false,
        std::Bool=false,
        impr_set::String="1"
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
    
    t = collect(1:Int64(HVPobs.Data.get_T(ens.id)))
    corrVec_PiPi = [corr_n(1,t,E,Z.^2,ens.L)]
    corrVec_JPi  = [corr_n(1,t,E,Z.*Zimpr,ens.L)]
    for n=2:nmax
        push!(corrVec_PiPi,corrVec_PiPi[n-1].+corr_n(n,t,E,Z,ens.L))
        push!(corrVec_JPi,corrVec_JPi[n-1].+corr_n(n,t,E,Zimpr,ens.L))
    end
    if !total
        return corrVec_PiPi, corrVec_JPi
    else
        if impr_set == "1"
            cv_l = cv_loc(ens.beta)
        elseif impr_set == "2"
            cv_l = cv_loc_set2(ens.beta)
        else
            error("Impr Set $impr_set not recoginsed, please choose between '1' and '2'.")
        end
        corrVec = Vector{Vector{uwreal}}(); [push!(corrVec,improve_corr_vkvk!(corrVec_PiPi[n], corrVec_JPi[n], 2*cv_l, std=std)[1:Int64(HVPobs.Data.get_T(ens.id)/2+1)]) for n in collect(1:nmax)]
        return  corrVec_PiPi, corrVec_JPi, corrVec
    end
end
reconstr_corr(ensid::String,E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};nmax::Union{Nothing,Int64}=nothing,FreeEN::Bool=false,total::Bool=false,std::Bool=false,impr_set::String="1") = reconstr_corr(EnsInfo(ensid),E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};nmax=nmax,FreeEN=FreeEN,total=total,std=std,impr_set=impr_set)

## test :

diag = "NLOa"
ens = EnsInfo("E250")
impr_set = "1"

E, Z, Z_impr = get_spectr_data(path_spectro,ens)

corrVec_PiPi, corrVec_JPi, corrVec = reconstr_corr(ens,E,Z,Z_impr,FreeEN=true,total=true,impr_set="1")


##---------------------------



fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_TMR_NLO"),"r") 
TMRDict = Dict{String, Any}()
partial_res = Vector{Dict}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    push!(partial_res,ALPHAdobs_read_next(fb, size=sz, keys=ks))
end
BDIO_close!(fb)
for dict in partial_res
    merge!(TMRDict, dict)
end
if diag == "NLOa&b"
    TMR = TMRDict["TMRa"] .+ TMRDict["TMRb"]
else
    TMR = TMRDict["TMR$(diag[end])"]
end


fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_corr_set"*impr_set),"r")
corr = Dict{String, Array{uwreal}}()
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    sz = tuple(d["size"]...)
    ks = collect(d["keys"])
    corr = ALPHAdobs_read_next(fb, size=sz, keys=ks)
end
BDIO_close!(fb)

corr_ll = corr["g33_ll"][1:Int64(HVPobs.Data.get_T(ens.id)/2 + 1)]
corr_lc = corr["g33_lc"][1:Int64(HVPobs.Data.get_T(ens.id)/2 + 1)]



##

uwerr.(corr_ll); uwerr.(corr_lc); [uwerr.(corr) for corr in corrVec]

t = collect(1:Int64(HVPobs.Data.get_T(ens.id)/2 + 1))

colors = ["blue","lightblue","green","limegreen","red","orange","yellow"]
legends = ["Data"]

title("Ens: $(ens.id)")
errorbar(t, value.(corr_ll), err.(corr_ll), color = "black", fmt="d", label=legends[1], capsize=2)
# errorbar(t, value.(corr_lc), err.(corr_lc), color = "", fmt="s", mfc="none", label="G lc", capsize=2)
for n=1:min(length(E),length(Z))
    errorbar(t, value.(corrVec[n]), err.(corrVec[n]), color=colors[n], fmt="o", mfc="none", label="Recons. (n=$n)", capsize=2)
    push!(legends,"Recons. (n=$n)")
end
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
yscale("log")
ylim([0.5*minimum(minimum.(value.(hcat([corr_ll,corrVec[1]]...)))),2*maximum(minimum.(value.(hcat([corr_ll,corrVec[min(length(E),length(Z))]]...))))])
xlabel("t/a")
ylabel(latexstring("G(t)"))
legend(legends, loc  = "best")
display(gcf())      #display the figure
close()

##

∆corr_ll = corrVec[end] - corr_ll; uwerr.(∆corr_ll)
∆corr_ll

##--------------------



for tcut in tcut0:t[end-1]  # compute the upper and lower corr bounds 

    if aens*tcut < tcut_fix  # we fix the eff energy at some point
        Eeff_=Eeff(tcut, obs)
    end
    if key[2:3] == "33"
        UB = mrho < E2pi ? corrBound(t, tcut, obs, ens, mrho) : corrBound(t, tcut, obs, ens, E2pi)
    elseif key[2:3] == "88"
        UB = corrBound(t, tcut, obs, ens, mrho)
    end
    LB = corrBound(t, tcut, obs, ens, Eeff_)

    UBInt = UB .* TMRw[tcut+1:end]
    LBint = LB .* TMRw[tcut+1:end]

    ub_amuNLO = (alpha/pi)^exp * (sum(int[1:tcut])+sum(UBInt)) * 1e10
    lb_amuNLO = (alpha/pi)^exp * (sum(int[1:tcut])+sum(LBint)) * 1e10

    push!(ub, ub_amuNLO)
    push!(lb, lb_amuNLO)
end

averb = (ub.+lb)./2; uwerr.(averb)
x0    = findfirst(abs.(value.(ub).-value.(lb)) .< 0.75.*err.(averb))
xend_x0 = findfirst(abs.(averb[x0:end].-averb[x0]) .> 0.5*err(averb[x0]))
xend_ = isnothing(xend_x0) ? x0+15  : x0-1 + xend_x0
if xend_-x0 > 5
    xend = xend_>length(averb) ? length(averb) : xend_
else
    xend=x0+3
    x0 > 3 ? x0-=3 : x0=1
end
plateau_fm = value(aens).*(collect(x0:xend).+tcut0.-2)
amu = sum(averb[x0:xend])/length(averb[x0:xend])

aux1 = sum(averb[x0:xend].^2)/length(averb[x0:xend])
aux2 = amu^2
syst = sqrt(abs(value(aux1 - aux2)))

HVP[key] = amu
HVPsyst[key] = syst
plateau[key] = [plateau_fm[1],plateau_fm[end]]
