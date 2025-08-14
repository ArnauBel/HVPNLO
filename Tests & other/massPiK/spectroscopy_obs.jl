######################################################################################
# This file was created by Alessandro Conigli 
# Here we compute the finite-a pion and kaon masses required in the chiral-continuum extrapolations
# The observables are then stored in a BDIO file, following the order 
########################################################################################

using Revise
using HVPobs
using PyPlot, LaTeXStrings
using OrderedCollections
using TimerOutputs
using BDIO, ADerrors
import ADerrors: err

#================ SET UP VARIABLES ===============#

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] =13
rcParams["axes.labelsize"] =22
rcParams["axes.titlesize"] = 18

plt.rc("text", usetex=true) # set to true if a LaTeX installation is present

#======== INCLUDE, PATHS AND CONSTANTS ==========#

include("./utils/const.jl")
include("./utils/types.jl")
include("./utils/tools.jl")
# include("data_management.jl")
include("./utils/plot_utils.jl")

path_data = "/Users/alessandroconigli/Lattice/data/HVP/mesons"
path_rw   = "/Users/alessandroconigli/Lattice/data/HVP/rwf_deflated"
path_ms   = "/Users/alessandroconigli/Lattice/data/HVP/ms_t0_dat"
path_fvc  = "/Users/alessandroconigli/Lattice/data/HVP/FVC"

path_plot = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/plots/ensembles/"
path_bdio = "/Users/alessandroconigli/MyDrive/postdoc-mainz/projects/deltalpha/data"

IMPR      = true
RENORM    = true
STD_DERIV = false

enslist = ["J306"]
ensinfo = EnsInfo.(enslist)


#============ OBSERVABLE ALLOCATIONS ============#

g_pi = Vector{Corr}(undef, length(ensinfo))
g_k  = Vector{Corr}(undef, length(ensinfo))

g_pi_a0p = Vector{Corr}(undef, length(ensinfo))
g_k_a0p = Vector{Corr}(undef, length(ensinfo))

obs = Vector(undef, length(ensinfo))

##


#================ READING DATA ====================#
wpmm = Dict{String, Vector{Float64}}()
wpmm["F300"]  = [5.0, -2.0, -1.0, -1.0]

@info("Reading data...")
@time begin
    
    for (k, ens) in enumerate(ensinfo)
        println("\n - Ensemble: $(ens.id)")

        println("   - Pion correlator")
        g_pi[k] = get_mesons_corr(path_data, ens, "ll", "PP", path_rw=path_rw, frw_bcwd=false, L=1 )
        g_pi_a0p[k] = get_mesons_corr(path_data, ens, "ll", "A0P", path_rw=path_rw, frw_bcwd=false, L=1 )
        println("   - Kaon correlator")
        g_k[k] = get_mesons_corr(path_data, ens, "ls", "PP", path_rw=path_rw, frw_bcwd=false, L=1 )
        g_k_a0p[k] = get_mesons_corr(path_data, ens, "ls", "A0P", path_rw=path_rw, frw_bcwd=false, L=1 )

        println("   - Gradient flow t0")
        obs[k] = OrderedDict()
        obs[k]["t0"] = get_t0(path_ms, ens, path_rw=path_rw, pl=true, wpm=wpmm)
    end
end

##
yy = g_pi_a0p[1].obs; uwerr.(yy); yy
errorbar(collect(1:length(yy)), value.(yy), err.(yy), fmt="s")
display(gcf())
close("all")
##
aa = meff(g_pi[1],[10,40], pl=true)
meff(g_k[1],[25,45], pl=true)
##
# =========== COMPUTING Meff ============#
wpmm = Dict{String, Vector{Float64}}()
wpmm["H101"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["H102r002"] = [5.0, -2.0, -1.0, -1.0]
wpmm["H400"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["N202"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N200"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N203"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N300"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["J303"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J304"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["F300"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J306"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J307"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["E300"]     = [5.0, -1.0, -1.0, -1.0]

for k in eachindex(ensinfo)
    if ensinfo[k].kappa_l == ensinfo[k].kappa_s
        obs[k]["mpi"] = get_meff_BMA(g_pi[k:k], ensinfo[k], path_plt=path_plot, ll=L"$m_{\pi}$", wpm=wpmm)[1] 
        obs[k]["mk"] = obs[k]["mpi"]
    else
        obs[k]["mpi"] = get_meff_BMA(g_pi[k:k], ensinfo[k], path_plt=path_plot, ll=L"$m_{\pi}$", wpm=wpmm)[1] 
        obs[k]["mk"]  = get_meff_BMA(g_k[k:k], ensinfo[k], path_plt=path_plot, ll=L"$m_{K}$", wpm=wpmm)[1] 
    end
end

## test with mpcac

mps = obs[1]["mpi"]
m12 = -  mpcac(g_pi_a0p[1].obs, g_pi[1].obs, [20,170], ca=-ca(ensinfo[1].beta))


fpi = dec_const_mpcac(m12, mps, ensinfo[1], renorm=true); uwerr(fpi); fpi
## test with decay constant
fpi = get_dec_const_BMA(g_pi[1], g_pi_a0p[1], obs[1]["mpi"], ensinfo[1], pl=true, ll = L"$f_{\pi}$") 
fpi_r = fpi * Za_l_sub(ensinfo[1].beta); uwerr(fpi_r); fpi_r
##

#======== SAVE TO BDIO =========#
for (k, o) in enumerate(obs)

    ens = ensinfo[k].id
    p = joinpath(path_bdio, ens,  string(ens, "_spectrum.bdio"))
    fb = BDIO_open(p, "w", ens)
    uinfo = 0 

    for (key, value) in o
        if isa(value, uwreal)
            write_uwreal(value, fb, uinfo)
        elseif isa(value, Vector{uwreal})
            [write_uwreal(v, fb, uinfo) for v in value]
        end
        uinfo +=1
    end

    BDIO_close!(fb)
end

## 
function dec_const_mpcac(m_pcac::uwreal, mps::uwreal, ens::EnsInfo; renorm::Bool=true)

    rave = ( ZP(ens.beta))/ mps^2 * m_pcac
    # rave = (2 )/ mps^2 * m_pcac

    if renorm
        ba = 1 + 0.0472 * (6/ens.beta)
        rave = Za_l_sub(ens.beta) * (1 + ba * m_pcac  ) * rave
    end

    return rave
end


## playing with masses with PBC
Thalf = Int64(length(g_pi[1].obs)/2) +1
T = length(g_pi[1].obs)
ydata = g_pi_a0p[1].obs[Thalf-20:Thalf+20] ; uwerr.(ydata)
xdata = collect(0:length(ydata)-1) .+ Thalf .-20 .-1
idx = findfirst(x->x==0.0, value.(ydata))

#ydata[idx].err = 1e-10  # uwreal([0.0, 1e-10], "test"); uwerr.(ydata)
ydata[idx] =  uwreal([0.0, 1e-10], "test"); uwerr.(ydata)

errorbar(xdata, value.(ydata), err.(ydata), fmt="s")
display(gcf())
close("all")

# @. model(x,p) =  p[2] * (-exp(-p[1]*(x)) + exp(-p[1]*(T-x)) )#  + p[4] * (exp(-(p[3]+p[1])*(x)) + exp(-(p[3]+p[1])*(Thalf*2-x)) )
@. model(x,p) =  p[2] * (exp(-p[1]*(x)) - exp(-p[1]*(T-x)) )#  + p[4] * (exp(-(p[3]+p[1])*(x)) + exp(-(p[3]+p[1])*(Thalf*2-x)) )

fit = fit_routine(model, xdata, ydata, 2, pval=true )
uwerr.(fit.param)
fit.param
#fit.pval
##




