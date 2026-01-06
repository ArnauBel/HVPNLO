# Import packages

include("HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using SpecialFunctions

using BDIO

using Plots
using PyPlot
using Colors

# include uwreal constants

# include("HVPtool/uwConst.jl")

# Path definition

julia_script_directory = @__DIR__

path_hvp_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "HVPData"),
    "SSD"   => joinpath(julia_script_directory, "..", "ObsExternal","PortableSSD","HVPData"),
    "clust" => joinpath(julia_script_directory, "..", "ObsExternal", "mogon_mount", "HVPdata")
)

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_meson = joinpath(julia_script_directory, "..", "HVPData", "meson_data")

path_rw    = joinpath(path_hvp_dict["local"], "reweight")
path_ms    = joinpath(path_hvp_dict["local"], "ms_t0_dat")

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

# path to ObsBDIO

path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")

# Ensamble choice

# Ens with problems: H105

ensList = ["A653","A654","B450","C101","C102","D150","D200","D201","D251","D450","D451","D452","E250","E300","F300","H101","H102","H200","H650","J303","J304","J306","J307","J500","J501","N101","N200","N202","N203","N300","N302","N451","N452","S400"]

# bc = obc
# ensList = ["C101","C102","D200","D201","E300","F300","H101","H102","H200","J303","J304","J306","J307","J500","J501","N101","N200","N202","N203","N300","N302","S400"]

# bc = pbc
# ensList = ["A653","A654","B450","D150","D251","D450","D451","D452","E250","H650","N451","N452"]

ensInfo = EnsInfo.(ensList)

wpmm = Dict{String, Vector{Float64}}()
wpmm["D450"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["H101"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["H102r002"] = [5.0, -2.0, -1.0, -1.0]
wpmm["H400"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["N202"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N200"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N203"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N300"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["N302"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["J303"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J304"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["F300"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J306"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J307"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J500"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["E300"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["A654"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["E300"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["H650"]     = [5.0, -2.0, -1.0, -1.0]

# Finite volume corrrections

function Z00(k::Integer)
    k < 0 && return 0  # No solutions for negative k
    bound = isqrt(k)    # Largest integer ≤ √k
    count = 0
    
    for m1 in -bound:bound
        for m2 in -bound:bound
            for m3 in -bound:bound
                if m1^2 + m2^2 + m3^2 == k
                    count += 1
                end
            end
        end
    end
    
    return -count  # Z00(0, k) = -r3(k)
end

m(k::Integer) = -Z00(k)

function g1(x::Float64;N2::Int64=40)
    return sum([4*m(n2)/(sqrt(n2)*x) * besselk(1,sqrt(n2)*x) for n2=1:N2])
end

# tree level for cross-check
function ∆Mπ_tl(ens::EnsInfo)
    return - 3/8π^2 * meson_ens[ens.id]["mPi"].mean^2 / (meson_ens[ens.id]["fPi"].mean^2*ens.L) * besselk(1,meson_ens[ens.id]["mPi"].mean*ens.L)
end
∆Mπ_tl(ensid::String) = ∆Mπ_tl(EnsInfo(ensid))

function ∆Fπ_tl(ens::EnsInfo)
    return 3/2π^2 * meson_ens[ens.id]["mPi"].mean / (meson_ens[ens.id]["fPi"].mean*ens.L) * besselk(1,meson_ens[ens.id]["mPi"].mean*ens.L)
end
∆Fπ_tl(ensid::String) = ∆Fπ_tl(EnsInfo(ensid))

# NLO
function ∆Mπ_NLO(ens::EnsInfo;N2::Int64=40)
    mπ = meson_ens[ens.id]["mPi"].mean
    mη = m_ens[ens.id]["mRho"].mean
    fπ = meson_ens[ens.id]["fPi"].mean
    L  = ens.L
    Nf = 2

    Xiπ = mπ^2 / (4π*fπ)^2
    Xiη = mη^2 / (4π*fπ)^2
    lambdaπ = mπ*L
    lambdaη = mη*L

    return mπ * (- 1/(2*Nf) * Xiπ * g1(lambdaπ,N2=N2)) + 1/12 * Xiη * g1(lambdaη,N2=N2)
end
∆Mπ_NLO(ensid::String) = ∆Mπ_NLO(EnsInfo(ensid))

function ∆MK_NLO(ens::EnsInfo;N2::Int64=40)
    mη = meson_ens[ens.id]["mPi"].mean
    mK = meson_ens[ens.id]["mK"].mean
    fπ = meson_ens[ens.id]["fPi"].mean
    L  = ens.L

    Xiη = mη^2 / (4π*fπ)^2
    lambdaη = mη*L

    return mK * (- 1/6 * Xiη * g1(lambdaη,N2=N2))
end
∆MK_NLO(ensid::String) = ∆MK_NLO(EnsInfo(ensid))

function ∆Fπ_NLO(ens::EnsInfo;N2::Int64=40)
    mπ = meson_ens[ens.id]["mPi"].mean
    mK = meson_ens[ens.id]["mK"].mean
    fπ = meson_ens[ens.id]["fPi"].mean
    L  = ens.L
    Nf = 2

    Xiπ = mπ^2 / (4π*fπ)^2
    XiK = mK^2 / (4π*fπ)^2
    lambdaπ = mπ*L
    lambdaK = mK*L

    return fπ * (Nf/2 * Xiπ * g1(lambdaπ,N2=N2) + 1/2 * XiK * g1(lambdaK,N2=N2))
end
∆Fπ_NLO(ensid::String) = ∆Fπ_NLO(EnsInfo(ensid))

function ∆Fπ_NNLO(ens::EnsInfo;N2::Int64=40)
    factor = hbarc * sqrt(t0sym(ens.beta))/sqrtt0_ph  
    mπ_ph = mpi_ph * (1e-3 / factor)  # from MeV to 1/a

    mπ = meson_ens[ens.id]["mPi"].mean
    fπ = meson_ens[ens.id]["fPi"].mean

    L  = ens.L

    Xiπ = mπ^2 / (4π*fπ)^2
    lambdaπ = mπ*L

    B0(n2) = 2besselk(1,sqrt(n2)*lambdaπ)
    B2(n2) = 2besselk(2,sqrt(n2)*lambdaπ)/(sqrt(n2)*lambdaπ)

    lbar_ph = [uwreal([-0.4,0.6],"l1bar ph"), uwreal([4.3,0.1],"l2bar ph"), uwreal([2.9,2.4],"l3bar ph"), uwreal([4.4,0.2],"l4bar ph")]
    lbar    = lbar_ph .+ 2log(mπ_ph/mπ)

    g = [2-π/2, π/4-1/2, 1/2-π/8, 3π/16-1/2]

    SFπ(n2) = 1/6*(8g[1] - 13g[2])*B0(n2) - 1/3*(40g[1]-12g[2]-8g[3]-13g[4])*B2(n2)

    IFπ(n2) = B0(n2)*(-7/9 + 2*lbar[1] + 4/3*lbar[2] - 3*lbar[4]) + B2(n2)*(112/9 - 8/3*lbar[1] - 32/3*lbar[2])  + SFπ(n2)

    RFπNLO = 1/lambdaπ * Xiπ^2 * sum([m(n2)/sqrt(n2) * IFπ(n2) for n2=1:N2])

    return ∆Fπ_NLO(ens,N2=N2) - fπ * RFπNLO
end
∆Fπ_NNLO(ensid::String) = ∆Fπ_NNLO(EnsInfo(ensid))

function RFK_NLO(ens::EnsInfo;N2::Int64=40)
    mπ = meson_ens[ens.id]["mPi"].mean
    mK = meson_ens[ens.id]["mK"].mean
    mη = m_ens[ens.id]["mRho"].mean
    fπ = meson_ens[ens.id]["fPi"].mean
    L  = ens.L
    # Nf = 2

    Xiπ = mπ^2 / (4π*fπ)^2
    XiK = mK^2 / (4π*fπ)^2
    Xiη = mη^2 / (4π*fπ)^2
    lambdaπ = mπ*L
    lambdaK = mK*L
    # lambdaη = mη*L

    # there is a -1 factor with respect to the usual definition
    return 3/8 * Xiπ * g1(lambdaπ,N2=N2) + 3/4 * XiK * g1(lambdaK,N2=N2) # + 3/8 * Xiη * g1(lambdaη,N2=N2)
end
RFK_NLO(ensid::String) = RFK_NLO(EnsInfo(ensid))

@info("Ready")

## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##

##==========================> Data reading and testing <==========================##

ens = ""

ens = EnsInfo(ens)

# use LMA and deflated when possible: 
p_mes = !isempty(filter(x->occursin(ens.id, x), readdir(joinpath(path_meson,"LMA"), join=true))) ? joinpath(path_meson,"LMA") : path_meson
p_rw  = !isempty(filter(x->occursin(ens.id, x), readdir(joinpath(path_rw,"reweight_deflated"), join=true))) ? joinpath(path_rw,"reweight_deflated") : path_rw

@info("Reading data")

println("   - Pion correlator")
corr_pi_pp  = get_mesons_corr(p_mes, ens, "ll", "PP", path_rw=p_rw, frw_bcwd=false, L=1 )
corr_pi_a0p = get_mesons_corr(p_mes, ens, "ll", "A0P", path_rw=p_rw, frw_bcwd=false, L=1 )

if ens.kappa_l != ens.kappa_s
    println("   - Kaon correlator")
    corr_k_pp  = get_mesons_corr(p_mes, ens, "ls", "PP", path_rw=p_rw, frw_bcwd=false, L=1 )
    corr_k_a0p = get_mesons_corr(p_mes, ens, "ls", "A0P", path_rw=p_rw, frw_bcwd=false, L=1 )
end

@info("Ready")

## pbc : 

t = collect(1:Int(length(corr_pi_pp.obs)/2 + 1))

obs_pp = -corr_pi_pp.obs[t]; uwerr.(obs_pp)
obs_a0p = corr_pi_a0p.obs[t]; uwerr.(obs_a0p)

fig = figure(figsize=(8,6))
errorbar(t, value.(obs_pp), err.(obs_pp), capsize=2, fmt="o", mfc="none", color="black", label = "Pi-PP 2pt. function")
errorbar(t, value.(obs_a0p), err.(obs_a0p), capsize=2, fmt="o", mfc="none", color="blue", label = "Pi-A0P 2pt. function")
# ylim(-0.02,0.05)
xlabel("t/a")
legend()
tight_layout()
display(gcf())
close()

## obc :

t = collect(1:Int(length(corr_pi_pp.obs)))

obs_pp  = -corr_pi_pp.obs[t]
obs_a0p = corr_pi_a0p.obs[t]

m_obs = meff(corr_pi_pp.obs); uwerr.(m_obs)

fig = figure(figsize=(8,6))
errorbar(t[2:end-2].+1/2, value.(m_obs), err.(m_obs), capsize=2, fmt="o", mfc="none", color="gray", label = "Pi eff. mass")

yscale("log")
# ylim(-0.02,0.05)
xlabel("t/a")
legend()
tight_layout()
display(gcf())
close()

## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##


##==========================> mPi & mK fit <==========================##

MESON = ["Pion","Kaon"]  #  ["Pion"]  ["Kaon"]  ["Pion","Kaon"]

ensid = ""

mPi_fitinfo = Dict{String,Dict{String,Any}}(
    "A653" => Dict{String,Any}("2state" => true , "plat" => [0.40,0.50], "mdof" => 4),
    # "A654" => Dict{String,Any}("2state" => true , "plat" => [0.30,0.40], "mdof" => 4),
    "A654" => Dict{String,Any}("2state" => false, "plat" => [0.50,0.60], "mdof" => 4),
    "H650" => Dict{String,Any}("2state" => false, "plat" => [0.40,0.80], "mdof" => 4),
    "D150" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.40], "mdof" => 4),
    "B450" => Dict{String,Any}("2state" => true , "plat" => [0.30,0.40], "mdof" => 4),
    "N452" => Dict{String,Any}("2state" => false, "plat" => [0.40,0.70], "mdof" => 4),
    "N451" => Dict{String,Any}("2state" => true , "plat" => [0.15,0.20], "mdof" => 4),
    "D450" => Dict{String,Any}("2state" => false, "plat" => [0.60,0.70], "mdof" => 4),
    "D451" => Dict{String,Any}("2state" => false, "plat" => [0.60,0.70], "mdof" => 4),
    "D452" => Dict{String,Any}("2state" => false, "plat" => [0.30,0.35], "mdof" => 4),
    "D251" => Dict{String,Any}("2state" => false, "plat" => [0.40,0.50], "mdof" => 4),
    "E250" => Dict{String,Any}("2state" => false, "plat" => [0.30,0.50], "mdof" => 4),

    "H101" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 25), # "mdof" => 20 for closer to LD paper mass
    "H102" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.75], "mdof" => 10),
    "N101" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.65], "mdof" => 5),
    "C101" => Dict{String,Any}("2state" => false, "plat" => [0.32,0.50], "mdof" => 5),
    "C102" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.75], "mdof" => 5),
    "S400" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.75], "mdof" => 30),
    "H200" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "N202" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 65),
    "N203" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 40),
    "N200" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 20),
    # "N200" => Dict{String,Any}("2state" => false, "plat" => [0.35,0.65], "mdof" => 20), # closer to LD paper
    "D200" => Dict{String,Any}("2state" => false, "plat" => [0.30,0.50], "mdof" => 10),
    "D201" => Dict{String,Any}("2state" => false, "plat" => [0.30,0.60], "mdof" => 10), # "mdof" too small ?
    "N300" => Dict{String,Any}("2state" => false, "plat" => [0.30,0.70], "mdof" => 10),
    "J307" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "N302" => Dict{String,Any}("2state" => false, "plat" => [0.30,0.65], "mdof" => 10),
    "J306" => Dict{String,Any}("2state" => false, "plat" => [0.22,0.75], "mdof" => 95),
    "J303" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "J304" => Dict{String,Any}("2state" => false, "plat" => [0.15,0.80], "mdof" => 90),
    # "E300" => Dict{String,Any}("2state" => false, "plat" => [0.10,0.70], "mdof" => 40), # wtf
    "E300" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.70], "mdof" => 40),
    "F300" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 90),
    "J500" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "J501" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 30),
)

mK_fitinfo = Dict{String,Dict{String,Any}}(
    "A654" => Dict{String,Any}("2state" => false, "plat" => [0.60,0.80], "mdof" => 4),
    "H650" => Dict{String,Any}("2state" => false , "plat" => [0.30,0.70], "mdof" => 4),
    "D150" => Dict{String,Any}("2state" => false, "plat" => [0.50,0.60], "mdof" => 4),
    "N452" => Dict{String,Any}("2state" => false, "plat" => [0.50,0.60], "mdof" => 4),
    "N451" => Dict{String,Any}("2state" => false, "plat" => [0.40,0.50], "mdof" => 4),
    "D450" => Dict{String,Any}("2state" => false, "plat" => [0.50,0.55], "mdof" => 4),
    "D451" => Dict{String,Any}("2state" => false, "plat" => [0.35,0.45], "mdof" => 4),
    "D452" => Dict{String,Any}("2state" => false, "plat" => [0.40,0.50], "mdof" => 4),
    "D251" => Dict{String,Any}("2state" => false, "plat" => [0.50,0.60], "mdof" => 4),
    "E250" => Dict{String,Any}("2state" => false, "plat" => [0.40,0.50], "mdof" => 4),

    "H102" => Dict{String,Any}("2state" => false,"plat" => [0.25,0.75], "mdof" => 10),
    "N101" => Dict{String,Any}("2state" => false,"plat" => [0.25,0.75], "mdof" => 5),
    "C101" => Dict{String,Any}("2state" => false,"plat" => [0.25,0.75], "mdof" => 5),
    "C102" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.75], "mdof" => 5),
    "S400" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.75], "mdof" => 5),
    "N203" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 40),
    "N200" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 20),
    "D200" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.65], "mdof" => 10),
    "D201" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.70], "mdof" => 30),
    "N302" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    # "J306" => Dict{String,Any}("2state" => false, "plat" => [0.10,0.75], "mdof" => 80),
    "J306" => Dict{String,Any}("2state" => false, "plat" => [0.22,0.75], "mdof" => 95),
    "J303" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 10),
    "J304" => Dict{String,Any}("2state" => false, "plat" => [0.15,0.80], "mdof" => 90),
    "E300" => Dict{String,Any}("2state" => false, "plat" => [0.25,0.70], "mdof" => 40),
    "F300" => Dict{String,Any}("2state" => false, "plat" => [0.10,0.70], "mdof" => 70),
    "J501" => Dict{String,Any}("2state" => false, "plat" => [0.20,0.80], "mdof" => 30),
)

AIC       = true  # always

PVAL      = false

PLOT      = true
WRITE     = false
OVERWRITE = false

path_bdio_w = path_bdio_dict["local"]

@info("Starting mP (and cP) calculation")

mDict = Dict()

@. cosh_model(x,p) = 0.0  # the function needs to be initialized

for ens in  [EnsInfo(ensid)]  # ensInfo  [EnsInfo(ensid)]

    println("- ens: $(ens.id)")
    println("   - Reading data...")

    # use LMA and deflated when possible: 
    p_mes = !isempty(filter(x->occursin(ens.id, x), readdir(joinpath(path_meson,"LMA"), join=true))) ? joinpath(path_meson,"LMA") : path_meson
    p_rw  = !isempty(filter(x->occursin(ens.id, x), readdir(joinpath(path_rw,"reweight_deflated"), join=true))) ? joinpath(path_rw,"reweight_deflated") : path_rw

    println("      - Pion correlator")
    corr_pi_pp = get_mesons_corr(p_mes, ens, "ll", "PP", path_rw=p_rw, frw_bcwd=false, L=1 )

    bc = ens.bc

    if bc == "pbc"
        frwd_bckwrd_symm!(corr_pi_pp)
    end

    if ens.kappa_l != ens.kappa_s
        println("      - Kaon correlator")
        corr_k_pp = get_mesons_corr(p_mes, ens, "ls", "PP", path_rw=p_rw, frw_bcwd=false, L=1 )

        if bc == "pbc"
            frwd_bckwrd_symm!(corr_k_pp)
        end
    end

    mDict[ens.id] = Dict{String,uwreal}()
    for meson in (ens.kappa_l != ens.kappa_s ? MESON : ["Pion"])
        if meson == "Pion"
            println("   - Starting Pion PP...")

            pl    = mPi_fitinfo[ens.id]["plat"]
            mdof  = mPi_fitinfo[ens.id]["mdof"]
            STATE = mPi_fitinfo[ens.id]["2state"]

            corr_pp = deepcopy(corr_pi_pp)
        elseif meson == "Kaon"
            println("   - Starting Kaon PP...")
            
                pl    = mK_fitinfo[ens.id]["plat"]
                mdof  = mK_fitinfo[ens.id]["mdof"]
                STATE = mK_fitinfo[ens.id]["2state"]

            corr_pp = deepcopy(corr_k_pp)
        else
            error("Meson type $meson not recognized.")
        end

        obs_pp = corr_pp.obs
        T = HVPobs.Data.get_T(corr_pp.id)

        if bc == "obc"
            println("      - Fitting for bc = obc...")

            m_obs = meff(obs_pp) # first and last data points are lost (-2) and derivative is taken (-1); length(meff_) = length(obs) - 3
            len = length(m_obs)

            @. const_model(x,p) = p[1] + 0*x

            plconst_vec_left  = collect(floor(Int64,pl[1]*len):2:ceil(Int64,len*(pl[1]+pl[2])/2)-ceil(Int64,mdof/2))
            plconst_vec_right = collect(ceil(Int64,len*(pl[1]+pl[2])/2)+floor(Int64,mdof/2):2:ceil(Int64,pl[2]*len))

            if isempty(plconst_vec_left) || isempty(plconst_vec_right)
                error("Not enough data to fit with d.o.f. ≥ $mdof. Decrease $mdof or increase the plateau search range. ")
            end

            wpm = corr_pp.id in keys(wpmm) ? wpmm : nothing

            fitconst_vec = Vector{FitRes}()
            pl_vec = Vector{Vector{Float64}}()
            for p0 in plconst_vec_left
                for pf in plconst_vec_right
                    # plateau = [p0,pl_f]
                    m_data = m_obs[p0:pf] 
                    fit = fit_routine(const_model,collect(p0:pf), m_data, 1, pval=PVAL, wpm=wpm, info=false, lineprint=false)
                    push!(fitconst_vec,fit)
                    push!(pl_vec,[p0,pf])
                end
            end

            fit_vec = fitconst_vec

        elseif bc == "pbc"
            println("      - Fitting for bc = pbc...")
            
            len = length(obs_pp)/2+1
            obs_pp = -obs_pp[1:Int(len)]

            if STATE
                @. cosh_model(x,p) = p[2] * (exp(-p[1]*x) + exp(-p[1]*(T-x))) +  p[4] * (exp(-p[3]*x) + exp(-p[3]*(T-x)))
                np = 4
            else
                @. cosh_model(x,p) = p[2] * (exp(-p[1]*x) + exp(-p[1]*(T-x)))
                np = 2
            end

            plcosh_vec = collect(floor(Int64,pl[1]*len):1:ceil(Int64,min(len-mdof,pl[2]*len)))
            pl_f = ceil(Int64,len)

            if isempty(plcosh_vec)
                error("Not enough data to fit with d.o.f. ≥ $mdof. Decrease $mdof or increase the plateau search range. ")
            end

            wpm = corr_pp.id in keys(wpmm) ? wpmm : nothing

            fitcosh_vec_pp = Vector{FitRes}()
            for p0 in plcosh_vec
                data_pp = obs_pp[p0:end] 

                fit_pp = fit_routine(cosh_model,collect(p0:pl_f).-1, data_pp, np, pval=PVAL, wpm=wpm, info=false, lineprint=false)
                push!(fitcosh_vec_pp,fit_pp)
            end

            p0_vec   = plcosh_vec
            fit_vec = fitcosh_vec_pp
        end

        println("      - Computing average...")

        w = get_w_from_fitres(vcat(fit_vec...), AIC=AIC)

        m_res_vec = (STATE && bc == "pbc") ? [[par[1],par[3]][argmin(value.([par[1],par[3]]))] for par in getfield.(vcat(fit_vec...),:param)] : [par[1] for par in getfield.(vcat(fit_vec...),:param)]

        m_res, m_sys = model_average(m_res_vec, w)
        m_res = m_res[1] + uwreal([0.0,m_sys],"mPP MA syst")

        if meson == "Pion"
            m_res += uwreal([1.0,0.5].*∆Mπ_NLO(ens),"mPi fvc syst")
            mkeys = ["mPi","cPi_PP"]
        elseif meson == "Kaon"
            m_res += uwreal([1.0,0.5].*∆MK_NLO(ens),"mK fvc syst")
            mkeys = ["mK","cK_PP"]
        end
        uwerr(m_res) 

        mDict[ens.id][mkeys[1]] = m_res

        println("         ⟹ mPP = $(print_uwreal(m_res))")

        if bc == "pbc"
            cpp_res_vec = STATE ? [[par[2],par[4]][argmin(value.([par[1],par[3]]))] for par in getfield.(vcat(fit_vec...),:param)] : [par[2] for par in getfield.(vcat(fit_vec...),:param)]

            cpp_res, cpp_sys = model_average(cpp_res_vec, w)
            cpp_res = cpp_res[1] + uwreal([0.0,cpp_sys],"cpp MA syst"); uwerr(cpp_res)

            mDict[ens.id][mkeys[2]] = cpp_res
            
            println("         ⟹ cPP = $(print_uwreal(cpp_res))")
        end

        bestW = 20

        if PLOT
            println("      - Plotting...")

            res = m_res
            if meson == "Pion"
                res_vec = m_res_vec .+ uwreal([1.0,0.5].*∆Mπ_NLO(ens),"mPi fvc syst")
            elseif meson == "Kaon"
                res_vec = m_res_vec .+ uwreal([1.0,0.5].*∆MK_NLO(ens),"mK fvc syst")
            end

            fig = figure(figsize=(16,12))

            if bc == "obc"

                gs = fig.add_gridspec(4, 1, height_ratios=[4, 1, 1, 1])  # Adjust the height_ratios as needed

                ax1 = fig.add_subplot(gs[1, 1])
                title("PP $(corr_pp.id) ($meson)")

                x0 = max(Int64(plconst_vec_left[1]-5),1):min(Int64(plconst_vec_right[end]+5),len)
                m_vec = m_obs[x0]; uwerr.(m_vec)

                errorbar(collect(x0).+1.5, value.(m_vec), err.(m_vec), fmt="o", capsize=2, color="black")
                fill_between([plconst_vec_left[1]-0.5,plconst_vec_right[end]+0.5].+1.5, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")
                
                WARG = sortperm(w,rev=true)[1:(length(w) >= bestW ? bestW : length(w))]
                uwerr.(res_vec[WARG])
                for (i,warg) in enumerate(WARG)
                    y_penal = i/length(WARG)
                    fill_between(pl_vec[warg].+1.5, value.(res_vec[warg])+y_penal*err.(res_vec[warg]), value.(res_vec[warg])-y_penal*err.(res_vec[warg]), alpha=w[warg], color="blue")
                end

                errorbar([plconst_vec_left[1],plconst_vec_left[end]].+1.5  , 0.5*(res.mean + maximum(value.(m_vec))).*[1,1], 0.01*(-res.mean + maximum(value.(m_vec))).*[1,1], fmt="", color="limegreen")
                errorbar([plconst_vec_right[1],plconst_vec_right[end]].+1.5, 0.5*(res.mean + maximum(value.(m_vec))).*[1,1], 0.01*(-res.mean + maximum(value.(m_vec))).*[1,1], fmt="", color="firebrick")

                axvline(x=plconst_vec_left[1]+1.0, color="cyan", linestyle="--")
                axvline(x=plconst_vec_right[end]+2.0, color="cyan", linestyle="--")

                axis("tight")
                xlabel(L"t/a")
                if meson == "Pion"
                    ylabel(L"$m^{\rm{eff}}_\pi(t)$")
                elseif meson == "Kaon"
                    ylabel(L"$m^{\rm{eff}}_K(t)$")
                end
                # ylim(res.mean-5*res.err,maximum(value.(m_vec) .+ 1.5 .* err.(m_vec)))

                ax2 = fig.add_subplot(gs[2, 1])
                # pl_vec[WARG]
                errorbar(collect(range(x0[1],x0[end],length(WARG))), value.(res_vec[WARG]), err.(res_vec[WARG]), fmt="d", mfc="none", color="blue")
                fill_between([x0[1],x0[end]], value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")

                if meson == "Pion"
                    ylabel(L"$m_\pi$")
                elseif meson == "Kaon"
                    ylabel(L"$m_K$")
                end
                setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels

                ax3 = fig.add_subplot(gs[3, 1])

                PyPlot.plot(collect(range(x0[1],x0[end],length(WARG))), w[WARG], linestyle="none", marker="o", mfc="none", color="blue")
                if AIC
                    ylabel(latexstring("\\rm{w}\\ [\\rm{AIC}]"))
                else
                    ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))
                end

                if PVAL
                    setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels

                    ax4 = fig.add_subplot(gs[4, 1])

                    pval_vec = getfield.(fit_vec,:pval)

                    PyPlot.plot(collect(range(x0[1],x0[end],length(WARG))), pval_vec[WARG], linestyle="none", marker="o", mfc="none", color="blue")

                    ylabel(L"$\rm{p-values}$")
                    ax4.set_xticks(collect(range(x0[1],x0[end],length(WARG))))
                    ax4.set_xticklabels(string.([Int.(pl) for pl in pl_vec[WARG]]))
                else
                    ax3.set_xticks(collect(range(x0[1],x0[end],length(WARG))))
                    ax3.set_xticklabels(string.([Int.(pl) for pl in pl_vec[WARG]]))
                end
                xlabel("Best fits")


                tight_layout()
                display(fig)
                close("all")

            elseif bc == "pbc"

                gs = fig.add_gridspec(5, 1, height_ratios=[4, 1, 1, 1, 1])  # Adjust the height_ratios as needed

                obs = obs_pp

                ax1 = fig.add_subplot(gs[1, 1])
                title("PP $(corr_pp.id) ($meson)")

                x0 = collect(max(floor(Int64,len*(pl[1]))-2,1):len) .- 1
                vec = obs[Int(x0[1])+1:end]; uwerr.(vec)

                errorbar(x0, value.(vec), err.(vec), fmt="o", capsize=2, color="black")
                maxw_arg = argmax(w[1:length(fit_vec)])
                par = getfield.(fit_vec,:param)[maxw_arg]

                x_ = collect(max(floor(Int64,len*(pl[1])),1)-0.5:0.1:len).-1
                ycosh_fit = cosh_model(x_,par); uwerr.(ycosh_fit)

                fill_between(x_, value.(ycosh_fit)+err.(ycosh_fit), value.(ycosh_fit)-err.(ycosh_fit), alpha=0.3, color="orange")

                axvline(x=p0_vec[1]-1.5, color="red", linestyle="--")
                axvline(x=p0_vec[end]-0.5, color="red", linestyle="--")

                axis("tight")
                if meson == "Pion"
                    ylabel(L"$G^{\rm{PP}}_\pi(t)$")
                elseif meson == "Kaon"
                    ylabel(L"$G^{\rm{PP}}_K(t)$")
                end

                p0_vec = Float64.(p0_vec) .- (1.5 + 1.0)
                
                ax21 = fig.add_subplot(gs[2, 1])
                uwerr.(res_vec)

                errorbar(p0_vec .+ 1.5, value.(res_vec), err.(res_vec), fmt="d", mfc="none", color="blue")
                fill_between(x0, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")

                if meson == "Pion"
                    ylabel(L"$m^{\rm{PP}}_\pi$")
                elseif meson == "Kaon"
                    ylabel(L"$m^{\rm{PP}}_K$")
                end
                setp(ax21.get_xticklabels(),visible=false) # Disable x tick labels

                ax22 = fig.add_subplot(gs[3, 1])
                uwerr.(cpp_res_vec)

                errorbar(p0_vec .+ 1.5, value.(cpp_res_vec), err.(cpp_res_vec), fmt="d", mfc="none", color="blue")
                fill_between(x0, value.(cpp_res)+err.(cpp_res), value.(cpp_res)-err.(cpp_res), alpha=0.4, color="green")

                if meson == "Pion"
                    ylabel(L"$C^{\rm{PP}}_\pi$")
                elseif meson == "Kaon"
                    ylabel(L"$C^{\rm{PP}}_K$")
                end
                setp(ax22.get_xticklabels(),visible=false) # Disable x tick labels


                ax3 = fig.add_subplot(gs[4, 1])

                fill_between(x0, maximum(w)/2, maximum(w)/2, alpha=0.0, color="white")
                PyPlot.plot(p0_vec .+ 1.5, w, linestyle="none", marker="o", mfc="none", color="blue")

                ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))

                if PVAL
                    setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels

                    ax4 = fig.add_subplot(gs[5, 1])

                    pval_vec = getfield.(fit_vec,:pval)
                    fill_between(x0, maximum(vcat(pval_vec...))/2, maximum(vcat(pval_vec...))/2, alpha=0.0, color="white")
                    PyPlot.plot(p0_vec .+ 1.5, pval_vec, linestyle="none", marker="o", mfc="none", color="blue")

                    ylabel(L"$\rm{p-values}$")
                end
                xlabel(L"$t/a$")


                tight_layout()
                display(fig)
                close()
            end
        end
    end # end meson loop
    if WRITE
        println("   - Printing BDIO...")

        p = create_path(path_bdio_w,["mass&dec",ens.id,"$(ens.id)_mPP"],OVERWRITE=OVERWRITE)

        io = IOBuffer()
        write(io, "$(ens.id) mPP")
        fb = ALPHAdobs_create(p, io)

        extra = Dict{String, Any}("ens" => ens.id)
        ALPHAdobs_write(fb, mDict[ens.id], extra=extra)
        
        ALPHAdobs_close(fb)
    end
end # end ens loop

## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##

##==========================> fPi & fK fit <==========================##


MESON = ["Pion","Kaon"]  #  ["Pion"]  ["Kaon"]  ["Pion","Kaon"]

ensid = "F300"

fPi_fitinfo = Dict{String,Dict{String,Any}}(
    "A653" => Dict{String,Any}("MPCAC_plat" => [4,42]  ),
    "A654" => Dict{String,Any}("MPCAC_plat" => [4,42]  ),
    "H650" => Dict{String,Any}("MPCAC_plat" => [4,90]  ),
    "D150" => Dict{String,Any}("MPCAC_plat" => [4,122] ),
    "B450" => Dict{String,Any}("MPCAC_plat" => [9,53]  ),
    "N452" => Dict{String,Any}("MPCAC_plat" => [10,115]),
    "N451" => Dict{String,Any}("MPCAC_plat" => [8,118] ),
    "D450" => Dict{String,Any}("MPCAC_plat" => [5,121] ),
    "D451" => Dict{String,Any}("MPCAC_plat" => [8,117] ),
    "D452" => Dict{String,Any}("MPCAC_plat" => [7,119] ),
    "D251" => Dict{String,Any}("MPCAC_plat" => [58,68] ),
    "E250" => Dict{String,Any}("MPCAC_plat" => [40,150]),
    
    # "H101" => Dict{String,Any}("MPCAC_plat" => [16,78] , "2state" => false, "plf" => 1.00, "plat" => [0.75,1.00], "mdof" => 4),
    "H101" => Dict{String,Any}("MPCAC_plat" => [16,78] , "2state" => false, "plf" => 1.00, "plat" => [0.50,0.60], "mdof" => 4),
    "H102" => Dict{String,Any}("MPCAC_plat" => [20,78] , "2state" => false, "plf" => 1.00, "plat" => [0.90,1.00], "mdof" => 4),
    "N101" => Dict{String,Any}("MPCAC_plat" => [25,95] , "2state" => false, "plf" => 1.00, "plat" => [0.85,1.00], "mdof" => 4),
    # "N101" => Dict{String,Any}("MPCAC_plat" => [25,95] , "2state" => false, "plf" => 1.00, "plat" => [0.50,1.00], "mdof" => 4),
    "C101" => Dict{String,Any}("MPCAC_plat" => [20,75] , "2state" => false, "plf" => 0.85, "plat" => [0.77,0.85], "mdof" => 4),
    "C102" => Dict{String,Any}("MPCAC_plat" => [17,76] , "2state" => false, "plf" => 1.00, "plat" => [0.77,1.00], "mdof" => 4),
    "S400" => Dict{String,Any}("MPCAC_plat" => [20,100], "2state" => false, "plf" => 1.00, "plat" => [0.75,1.00], "mdof" => 4),
    "H200" => Dict{String,Any}("MPCAC_plat" => [15,80] , "2state" => false, "plf" => 1.00, "plat" => [0.60,1.00], "mdof" => 4),
    # "N202" => Dict{String,Any}("MPCAC_plat" => [15,105], "2state" => false, "plf" => 1.00, "plat" => [0.45,1.00], "mdof" => 4),
    "N202" => Dict{String,Any}("MPCAC_plat" => [15,105], "2state" => true , "plf" => 1.00, "plat" => [0.20,0.25], "mdof" => 4),
    "N203" => Dict{String,Any}("MPCAC_plat" => [15,105], "2state" => false, "plf" => 0.85, "plat" => [0.70,0.85], "mdof" => 4),
    "N200" => Dict{String,Any}("MPCAC_plat" => [15,105], "2state" => false, "plf" => 0.95, "plat" => [0.75,0.95], "mdof" => 4),
    "D200" => Dict{String,Any}("MPCAC_plat" => [14,106], "2state" => false, "plf" => 1.00, "plat" => [0.85,1.00], "mdof" => 4),
    "D201" => Dict{String,Any}("MPCAC_plat" => [14,106], "2state" => false, "plf" => 1.00, "plat" => [0.75,1.00], "mdof" => 4),
    # "N300" => Dict{String,Any}("MPCAC_plat" => [15,105], "2state" => false, "plf" => 0.70, "plat" => [0.60,0.70], "mdof" => 4),
    "N300" => Dict{String,Any}("MPCAC_plat" => [15,105], "2state" => false, "plf" => 1.00, "plat" => [0.50,0.80], "mdof" => 4),
    "J307" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.55,1.00], "mdof" => 4),    
    "N302" => Dict{String,Any}("MPCAC_plat" => [15,105], "2state" => false, "plf" => 1.00, "plat" => [0.80,1.00], "mdof" => 4),
    "J306" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.40,1.00], "mdof" => 4),    
    "J303" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.50,1.00], "mdof" => 4),    
    "J304" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.80,1.00], "mdof" => 4),    
    "E300" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.85,1.00], "mdof" => 4),    
    "F300" => Dict{String,Any}("MPCAC_plat" => [20,230], "2state" => false, "plf" => 1.00, "plat" => [0.67,1.00], "mdof" => 4),    
    "J500" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.50,1.00], "mdof" => 4),    
    "J501" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.50,1.00], "mdof" => 4),
    )

fK_fitinfo = Dict{String,Dict{String,Any}}(
    "A654" => Dict{String,Any}("MPCAC_plat" => [4,42]  ),
    "H650" => Dict{String,Any}("MPCAC_plat" => [4,90]  ),
    "D150" => Dict{String,Any}("MPCAC_plat" => [10,116]),
    "N452" => Dict{String,Any}("MPCAC_plat" => [10,115]),
    "N451" => Dict{String,Any}("MPCAC_plat" => [8,118] ),
    "D450" => Dict{String,Any}("MPCAC_plat" => [23,103]),
    "D451" => Dict{String,Any}("MPCAC_plat" => [8,117] ),
    "D452" => Dict{String,Any}("MPCAC_plat" => [7,119] ),
    "D251" => Dict{String,Any}("MPCAC_plat" => [58,68] ),
    "E250" => Dict{String,Any}("MPCAC_plat" => [40,150]),
    
    "H102" => Dict{String,Any}("MPCAC_plat" => [20,75] , "2state" => false, "plf" => 1.00, "plat" => [0.90,1.00], "mdof" => 4),
    "N101" => Dict{String,Any}("MPCAC_plat" => [25,95] , "2state" => false, "plf" => 1.00, "plat" => [0.85,1.00], "mdof" => 4),
    "C101" => Dict{String,Any}("MPCAC_plat" => [20,75] , "2state" => false, "plf" => 0.85, "plat" => [0.77,0.85], "mdof" => 4),
    "C102" => Dict{String,Any}("MPCAC_plat" => [17,76] , "2state" => false, "plf" => 1.00, "plat" => [0.75,1.00], "mdof" => 4),
    "S400" => Dict{String,Any}("MPCAC_plat" => [20,100], "2state" => false, "plf" => 1.00, "plat" => [0.75,1.00], "mdof" => 4),
    "N203" => Dict{String,Any}("MPCAC_plat" => [15,105], "2state" => false, "plf" => 0.85, "plat" => [0.70,0.85], "mdof" => 4),
    "N200" => Dict{String,Any}("MPCAC_plat" => [15,105], "2state" => false, "plf" => 0.95, "plat" => [0.75,0.95], "mdof" => 4),
    "D200" => Dict{String,Any}("MPCAC_plat" => [14,106], "2state" => false, "plf" => 1.00, "plat" => [0.85,1.00], "mdof" => 4),
    "D201" => Dict{String,Any}("MPCAC_plat" => [14,106], "2state" => false, "plf" => 1.00, "plat" => [0.75,1.00], "mdof" => 4),
    "N302" => Dict{String,Any}("MPCAC_plat" => [15,105], "2state" => false, "plf" => 1.00, "plat" => [0.80,1.00], "mdof" => 4),
    "J306" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.40,1.00], "mdof" => 4),    
    "J303" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.50,1.00], "mdof" => 4),    
    "J304" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.80,1.00], "mdof" => 4),    
    "E300" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.85,1.00], "mdof" => 4),    
    "F300" => Dict{String,Any}("MPCAC_plat" => [20,230], "2state" => false, "plf" => 1.00, "plat" => [0.70,1.00], "mdof" => 4),    
    "J501" => Dict{String,Any}("MPCAC_plat" => [20,170], "2state" => false, "plf" => 1.00, "plat" => [0.50,1.00], "mdof" => 4),  
)

AIC       = true  # always

PVAL      = false

PLOT      = [true,true]
WRITE     = false
OVERWRITE = false

path_bdio_w = path_bdio_dict["local"]

@info("Starting fP calculation")

fDict = Dict()

for ens in [EnsInfo(ensid)]  # ensInfo  [EnsInfo(ensid)]

    println("- ens: $(ens.id)")
    println("   - Reading data...")

    # use LMA and deflated when possible: 
    p_mes = !isempty(filter(x->occursin(ens.id, x), readdir(joinpath(path_meson,"LMA"), join=true))) ? joinpath(path_meson,"LMA") : path_meson
    p_rw  = !isempty(filter(x->occursin(ens.id, x), readdir(joinpath(path_rw,"reweight_deflated"), join=true))) ? joinpath(path_rw,"reweight_deflated") : path_rw
    
    println("      - Pion correlator")
    corr_pi_pp  = get_mesons_corr(p_mes, ens, "ll", "PP" , path_rw=p_rw, frw_bcwd=false, L=1 )
    corr_pi_a0p = get_mesons_corr(p_mes, ens, "ll", "A0P", path_rw=p_rw, frw_bcwd=false, L=1 )

    bc = ens.bc

    if bc == "pbc"
        frwd_bckwrd_symm!(corr_pi_pp)
        frwd_bckwrd_antisymm!(corr_pi_a0p)
    end

    if ens.kappa_l != ens.kappa_s
        println("      - Kaon correlator")
        corr_k_pp  = get_mesons_corr(p_mes, ens, "ls", "PP" , path_rw=p_rw, frw_bcwd=false, L=1 )
        corr_k_a0p = get_mesons_corr(p_mes, ens, "ls", "A0P", path_rw=p_rw, frw_bcwd=false, L=1 )

        if bc == "pbc"
            frwd_bckwrd_symm!(corr_k_pp)
            frwd_bckwrd_antisymm!(corr_k_a0p)
        end
    end
    fDict[ens.id] = Dict{String,uwreal}()
    for meson in (ens.kappa_l != ens.kappa_s ? MESON : ["Pion"])
        if meson == "Pion"
            println("   - Starting Pion PS...")

            if bc != "pbc"
                pl    = fPi_fitinfo[ens.id]["plat"]
                plf   = fPi_fitinfo[ens.id]["plf"]
                mdof  = fPi_fitinfo[ens.id]["mdof"]
                STATE = fPi_fitinfo[ens.id]["2state"]
            end
            plrs = fPi_fitinfo[ens.id]["MPCAC_plat"]

            corr_pp = deepcopy(corr_pi_pp)
            corr_a0p = deepcopy(corr_pi_a0p)
        elseif meson == "Kaon"
            println("   - Starting Kaon PS...")
            
            if bc != "pbc"
                pl    = fK_fitinfo[ens.id]["plat"]
                plf   = fK_fitinfo[ens.id]["plf"]
                mdof  = fK_fitinfo[ens.id]["mdof"]
                STATE = fK_fitinfo[ens.id]["2state"]
            end
            plrs = fK_fitinfo[ens.id]["MPCAC_plat"]

            corr_pp = deepcopy(corr_k_pp)
            corr_a0p = deepcopy(corr_k_a0p)
        else
            error("Meson type $meson not recognized.")
        end

        obs_pp = corr_pp.obs
        obs_a0p = corr_a0p.obs
        T = HVPobs.Data.get_T(corr_pp.id)

        mrs = mpcac(-obs_a0p, obs_pp, plrs, ca=ca(ens.beta), pl=PLOT[1])
        
        mDict = BDIOread_mPP(path_bdio_w,ens.id)
        if meson == "Pion"
            mP = mDict["mPi"] - uwreal([1.0,0.5].*∆Mπ_NLO(ens),"mPi fvc syst")
            if bc == "pbc"
                cP = mDict["cPi_PP"]
            end
        elseif meson == "Kaon"
            mP = mDict["mK"]  - uwreal([1.0,0.5].*∆MK_NLO(ens),"mK fvc syst")
            if bc == "pbc"
                cP = mDict["cK_PP"]
            end
        end

        improve_corr_vkvk!(obs_a0p, obs_pp, ca(ens.beta))

        if bc == "obc"
            println("      - Fitting for bc = obc...")

            _, y0_ens = findmax(value.(-obs_pp)) # get the source
            y0 = y0_ens-1

            R = sqrt(2/mP) .* (((obs_a0p .* reverse(obs_a0p)) ./ (obs_pp[T-y0]) ).^2).^0.25
            R = R[2:Int(T/2)-1]

            len = length(R)

            @. const_model(x,p) = !STATE ? p[1] + 0*x : p[1] + p[2]*exp(p[3]*x)
            np = !STATE ? 1 : 3

            pl_f = minimum([ceil(Int64,plf*len),len])
            p0_vec = collect(floor(Int64,pl[1]*len):1:ceil(Int64,min(pl_f-mdof,pl[2]*len)))

            if isempty(p0_vec)
                error("Not enough data to fit with d.o.f. ≥ $mdof. Decrease $mdof or increase the plateau search range. ")
            end

            wpm = corr_a0p.id in keys(wpmm) ? wpmm : nothing

            fit_vec = Vector{FitRes}()

            for p0 in p0_vec
                R_data = R[p0:pl_f]

                fit_ = fit_routine(const_model,collect(p0:pl_f).-1, R_data, np, pval=PVAL, wpm=wpm, info=false, lineprint=false)
                push!(fit_vec,fit_)
            end

            println("      - Computing average...")

            w = get_w_from_fitres(vcat(fit_vec...), AIC=AIC)

            Rres_vec = [par[1] for par in getfield.(vcat(fit_vec...),:param)]

            R_res, R_sys = model_average(Rres_vec, w)

            R_res = R_res[1] + uwreal([0.0,R_sys],"Ratio MA syst"); uwerr(R_res)
        elseif bc == "pbc"
            println("      - Fitting for bc = pbc...")

            # len_pp = Int64(length(obs_pp)/2+1)

            # obs_pp  = obs_pp[1:len_pp]
            # obs_a0p = obs_a0p[1:len_pp]
            
            # x_ = collect(1:len_pp) .- 1.
            # R = sqrt(2/mP) .* (obs_a0p ./ (obs_pp.^2).^(1/4)) .* (sqrt.(exp.(-mP.*x_) .+ exp.(-mP.*(T.-x_)))./(exp.(-mP.*x_) .- exp.(-mP.*(T.-x_))))
            # R = R[1:end-1] # -1 to avoid the NaN

            ZP = sqrt(2mP*cP)

            R_res = 2ZP/mP^2 * mrs
        end

        ba = 1 + 0.0472 * (6/ens.beta)
        f_res = Za_l_sub(ens.beta) * (1 + ba * mrs) * R_res

        if meson == "Pion"
            f_res += uwreal([1.0,0.5].*∆Fπ_NNLO(ens).mean,"fPi fvc syst")
            fkey = "fPi"
        elseif meson == "Kaon"
            f_res += uwreal([1.0,0.5].*(f_res.mean*RFK_NLO(ens)),"fK fvc syst")
            fkey = "fK"
        end
        uwerr(f_res)

        fDict[ens.id][fkey] = f_res

        println("         ⟹ fPP = $(print_uwreal(f_res))")

        if PLOT[2] && bc == "obc"
            println("      - Plotting...")

            bestW = 20

            fig = figure(figsize=(16,12))

            obs = R
            res_vec = Za_l_sub(ens.beta) * (1 + ba * mrs) .* Rres_vec
            if meson == "Pion"
                res_vec .+= ∆Fπ_NLO(ens)
            end
            res = f_res

            gs = fig.add_gridspec(4, 1, height_ratios=[4, 1, 1, 1])  # Adjust the height_ratios as needed

            ax1 = fig.add_subplot(gs[1, 1])
            title("$(corr_a0p.id) ($meson)")

            x0 = collect(max(floor(Int64,len*(pl[1]))-5,1):len) .- 1
            vec = obs[Int(x0[1])+1:end]; uwerr.(vec)

            errorbar(x0, value.(vec), err.(vec), fmt="o", capsize=2, color="black")
            maxw_arg = argmax(w[1:length(fit_vec)])
            par = getfield.(fit_vec,:param)[maxw_arg]

            if plf == 1.0
                x_ = collect(max(floor(Int64,len*(pl[1])),1)-0.4:0.1:pl_f).-1
            else
                x_ = collect(max(floor(Int64,len*(pl[1])),1)-0.4:0.1:pl_f+0.6).-1
            end
            yconst_fit = const_model(x_,par); uwerr.(yconst_fit)


            fill_between(x_, value.(yconst_fit)+err.(yconst_fit), value.(yconst_fit)-err.(yconst_fit), alpha=0.3, color="orange")

            axvline(x=p0_vec[1]-1.0-0.4, color="red", linestyle="--")
            axvline(x=p0_vec[end]-1.0+0.4, color="red", linestyle="--")
            if plf != 1.0
                axvline(x=pl_f-1.0+0.6, color="orange", linestyle="dotted")
            end

            axis("tight")
            if meson == "Pion"
                ylabel(L"$R_\pi(t)$")
            elseif meson == "Kaon"
                ylabel(L"$R_K(t)$")
            end
            
            ax2 = fig.add_subplot(gs[2, 1])
            uwerr.(res_vec)

            errorbar(p0_vec .- 1.0, value.(res_vec), err.(res_vec), fmt="d", mfc="none", color="blue")
            fill_between(x0, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")

            if meson == "Pion"
                ylabel(L"$f_\pi$")
            elseif meson == "Kaon"
                ylabel(L"$f_K$")
            end
            setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels


            ax3 = fig.add_subplot(gs[3, 1])

            fill_between(x0, maximum(w)/2, maximum(w)/2, alpha=0.0, color="white")
            PyPlot.plot(p0_vec .- 1., w, linestyle="none", marker="o", mfc="none", color="blue")

            ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))
            
            if PVAL
                setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels
            
                ax4 = fig.add_subplot(gs[4, 1])

                pval_vec = getfield.(fit_vec,:pval)
                fill_between(x0, maximum(vcat(pval_vec...))/2, maximum(vcat(pval_vec...))/2, alpha=0.0, color="white")
                PyPlot.plot(p0_vec .- 1.0, pval_vec, linestyle="none", marker="o", mfc="none", color="blue")

                ylabel(L"$\rm{p-values}$")
            end
            xlabel(L"$t/a$")


            tight_layout()
            display(fig)
            close("all")
        end
    end # end mesoon loop
    if WRITE
        println("   - Printing BDIO...")

        p = create_path(path_bdio_w,["mass&dec",ens.id,"$(ens.id)_fPS"],OVERWRITE=OVERWRITE)

        io = IOBuffer()
        write(io, "$(ens.id) fPS")
        fb = ALPHAdobs_create(p, io)

        extra = Dict{String, Any}("ens" => ens.id)
        ALPHAdobs_write(fb, fDict[ens.id], extra=extra)
        
        ALPHAdobs_close(fb)
    end
end # end ens loop

## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##


##==========================> Reading test <==========================##

ensid = ""

path_bdio_r = path_bdio_dict["local"]

mDict = BDIOread_mPP(path_bdio_r,ensid)

fDict = BDIOread_fPS(path_bdio_r,ensid)

@info("Data ready")

## <<----------------------------------------------------------------------------------------------------------------------->> ##
## <<----------------------------------------------------------------------------------------------------------------------->> ##

## rough interpolation (usefull to approximate FVC for new ensembles)

path_bdio_r = path_bdio_dict["local"]

for ens in ensInfo
    mDict = BDIOread_mPP(path_bdio_r,ens.id); uwerr(mDict["mPi"])
    fDict = BDIOread_fPS(path_bdio_r,ens.id); uwerr(fDict["fPi"])

    println("ens $(ens.id) => mPi fvc/err = $(round(∆Mπ_NLO(ens)/mDict["mPi"].err,digits=2)); fPi fvc/err = $(round(∆Fπ_NLO(ens)/fDict["fPi"].err,digits=2))")
end

## 

y1 = m_ens["A653"]["mK"]
y2 = m_ens["A654"]["mK"]

x1 = 0.1365716
x2 = 0.136750
x  = 0.136850
                   
par, _ = lin_fit([x1,x2],[y1,y2],wpm=wpmm,lineprint=false)
y = y_lin_fit(par,x)



