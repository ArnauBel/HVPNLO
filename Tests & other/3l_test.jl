# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

using QuadGK

using ProgressBars
using Suppressor
using TimerOutputs

using PyPlot

# Include necessary functions

include("../HVPtools/Utils.jl")

include("../HVPtools/Reader.jl")

include("../HVPtools/Writer.jl")

include("../HVPtools/TMRKernel.jl")


# BDIO path definition (set 'STD_DERIV = true' to use the standard sym. derivative in the impr.)

julia_script_directory = @__DIR__

path_bdio = joinpath(julia_script_directory, "..", "..", "ObsBDIO")
path_tl   = joinpath(julia_script_directory, "..", "..", "LMEData", "tree_level_improv")
path_coef = joinpath(julia_script_directory, "..", "..", "KernelCoeff")

# Blind analysis (Simon K.) safe ensembles: 
# SU(3) sym.  H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

# All considered ensembles are:
ensList = ["A653","A654","B450","C101","C102","D150","D200","D201","D450","D451","D452","E250","E300","F300","H101","H102","H200","J303","J304","J306","J307","J500","J501","N101","N200","N202","N203","N300","N302","N451","N452","S400"]

ensInfo = EnsInfo.(ensList)

# We do not have charm or disconnected data for some of the ensembles
ensNOcharm = ["J501","N451","D150","D451","J304","C102","D251","D201","J306","J307","F300","H200","N452"]
ensNOdisc  = ["D251","J306","J307","F300","D450"]

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##


##==========================> Tree-level computation <==========================##


diag = "LO"  # LO  NLOa  NLOb  NLOa&b
IMPR_SET = ["1","1old","2"]

Qlist = [3.5, 4.0, 5.0, 6.0, 7.0, 8.0]

STD_DERIV = false

MASSLESS  = false

C4 = Dict("LO" => π^2/9, "NLOa" => NaN, "NLOb" => NaN, "NLOa&b" => NaN) # to be added in constants

exp_diag = diag == "LO" ? 2 : 3

t0 = Dict()
HVPQ_ll = Dict(); HVPQ_lc = Dict()
HVPQ_ll_P = Dict(); HVPQ_lc_P = Dict()
HVPQ_ll_NP = Dict(); HVPQ_lc_NP = Dict()
for ens in ensInfo
    println("ens: $(ens.id)")

    corr33tl_ll, corr33tl_lc = read_tree_level_v33(path_tl, cons=true)
    corr33tl_v3s03_ll, corr33tl_v3s03_lc = read_tree_level_v3sig03(path_tl, cons=true, massless=MASSLESS)

    corr33tl_ll = uwreal.(corr33tl_ll)
    corr33tl_lc = uwreal.(corr33tl_lc)


    corr33tl_ll_P = deepcopy(corr33tl_ll)
    corr33tl_lc_P = deepcopy(corr33tl_lc)

    cv_l = 0.0
    cv_c = 0.5

    improve_corr_vkvk!(corr33tl_ll_P, corr33tl_v3s03_ll, 2*cv_l, std=STD_DERIV, treelevel=true)
    improve_corr_vkvk_cons!(corr33tl_lc_P, corr33tl_v3s03_ll, corr33tl_v3s03_lc, cv_l, cv_c, std=STD_DERIV, treelevel=true)

    corr33tl_ll_NP = Dict(); corr33tl_lc_NP = Dict()
    for impr_set in IMPR_SET
        corr33tl_ll_NP[impr_set] = deepcopy(corr33tl_ll)
        corr33tl_lc_NP[impr_set] = deepcopy(corr33tl_lc)

        beta = ens.beta
        if impr_set == "1"
            cv_l = cv_loc(beta)
            cv_c = cv_cons(beta)
        elseif impr_set =="1old"
            cv_l = cv_loc_old(beta)
            cv_c = cv_cons_old(beta)
        elseif impr_set =="2"
            cv_l = cv_loc_set2(beta)
            cv_c = cv_cons_set2(beta)
        end

        improve_corr_vkvk!(corr33tl_ll_NP[impr_set], corr33tl_v3s03_ll, 2*cv_l, std=STD_DERIV, treelevel=true)
        improve_corr_vkvk_cons!(corr33tl_lc_NP[impr_set], corr33tl_v3s03_ll, corr33tl_v3s03_lc, cv_l, cv_c, std=STD_DERIV, treelevel=true)
    end

    TMR        = BDIOread_TMR(path_bdio,ens,diag,SU3=false)
    t0[ens.id] = BDIOread_t0(path_bdio,ens)
    factor     = hbarc * sqrt(t0[ens.id])/sqrtt0_ph

    sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1); t = collect(1:sym_points)
    aens = sqrtt0_ph / sqrt(t0[ens.id]); tfm = aens.*(t.-1)

    TMRw = TMR .* Window("SD")(tfm)
    TMRb(Q::Float64) = ((16/(Q/factor)^2)^2 * C4[diag] * (massmu/factor)^2) .* sin.((Q/factor/4) .* (t .- 1)).^4
    TMRsub(Q::Float64) = TMRw .- (Window("SD")(0) .* TMRb(Q))

    HVPQ_ll[ens.id] = []; HVPQ_lc[ens.id] = []
    HVPQ_ll_P[ens.id] = []; HVPQ_lc_P[ens.id] = []
    HVPQ_ll_NP[ens.id] = Dict(); HVPQ_lc_NP[ens.id] = Dict()
    [(HVPQ_ll_NP[ens.id][impr_set] = []; HVPQ_lc_NP[ens.id][impr_set] = []) for impr_set in IMPR_SET]
    for Q in Qlist
        int_ll = value.(corr33tl_ll[t]) .* value.(TMRsub(Q))
        int_lc = value.(corr33tl_lc[t]) .* value.(TMRsub(Q))
        int_ll_P = value.(corr33tl_ll_P[t]) .* value.(TMRsub(Q))
        int_lc_P = value.(corr33tl_lc_P[t]) .* value.(TMRsub(Q))

        amu_ll = (alpha/pi)^exp_diag * sum(int_ll) * 1e10
        amu_lc = (alpha/pi)^exp_diag * sum(int_lc) * 1e10
        amu_ll_P = (alpha/pi)^exp_diag * sum(int_ll_P) * 1e10
        amu_lc_P = (alpha/pi)^exp_diag * sum(int_lc_P) * 1e10

        push!(HVPQ_ll[ens.id],amu_ll/2)
        push!(HVPQ_lc[ens.id],amu_lc/2)
        push!(HVPQ_ll_P[ens.id],amu_ll_P/2)
        push!(HVPQ_lc_P[ens.id],amu_lc_P/2)

        for impr_set in IMPR_SET
            int_ll_NP = value.(corr33tl_ll_NP[impr_set][t]) .* value.(TMRsub(Q))
            int_lc_NP = value.(corr33tl_lc_NP[impr_set][t]) .* value.(TMRsub(Q))

            amu_ll_NP = (alpha/pi)^exp_diag * sum(int_ll_NP) * 1e10
            amu_lc_NP = (alpha/pi)^exp_diag * sum(int_lc_NP) * 1e10

            push!(HVPQ_ll_NP[ens.id][impr_set],amu_ll_NP/2)
            push!(HVPQ_lc_NP[ens.id][impr_set],amu_lc_NP/2)
        end
    end
end

HVP3l0 = compute_HVPtl0(diag,"SDsub",Qlist,path_coef)

##

ens = EnsInfo("H101")
impr_set = "2"

corr33tl_ll, corr33tl_lc = read_tree_level_v33(path_tl, cons=true)
corr33tl_v3s03_ll, corr33tl_v3s03_lc = read_tree_level_v3sig03(path_tl, cons=true, massless=true)

corr33tl_ll = uwreal.(corr33tl_ll)
corr33tl_lc = uwreal.(corr33tl_lc)

beta = ens.beta
if impr_set == "1"
    cv_l = cv_loc(beta)
    cv_c = cv_cons(beta)
elseif impr_set =="1old"
    cv_l = cv_loc_old(beta)
    cv_c = cv_cons_old(beta)
elseif impr_set =="2"
    cv_l = cv_loc_set2(beta)
    cv_c = cv_cons_set2(beta)
end

improve_corr_vkvk!(corr33tl_ll, corr33tl_v3s03_ll, 2*cv_l, std=false, treelevel=true)
improve_corr_vkvk_cons!(corr33tl_lc, corr33tl_v3s03_ll, corr33tl_v3s03_lc, cv_l, cv_c, std=false, treelevel=true)

corr33tl_ll
corr33tl_lc
##

Q = 5.0 # GeV
IMPR_SET = ["1old","2"]


Qarg = Qlist .== Q

color_dict = Dict(
    "1old" => Dict("ll" => "blue",  "lc" => "red"  ),
    "1"    => Dict("ll" => "blue",  "lc" => "red"  ),
    "2"    => Dict("ll" => "green", "lc" => "brown")
) # ; color_dict["1old"] = color_dict["1"]
fmt_dict =  Dict(
    "1old" => Dict("ll" => "^",  "lc" => "o"  ),
    "1"    => Dict("ll" => "^",  "lc" => "o"  ),
    "2"    => Dict("ll" => "s", "lc" => "d")
)

fig = figure(figsize=(8,6))

xdata    = value.([1 / (8*t0[ens.id]) for ens in ensInfo])
arg = sortperm(xdata)

ydata_ll = [HVPQ_ll[ens.id][Qarg][1] for ens in ensInfo]
ydata_lc = [HVPQ_lc[ens.id][Qarg][1] for ens in ensInfo]
ydata_ll_P = [HVPQ_ll_P[ens.id][Qarg][1] for ens in ensInfo]
ydata_lc_P = [HVPQ_lc_P[ens.id][Qarg][1] for ens in ensInfo]

ydata_ll

PyPlot.plot([0.0],HVP3l0[Qarg], color = "black", marker = "o", label = "Continuum")
PyPlot.plot(xdata[arg],ydata_ll[arg], color = "orange", marker = "x", mfc="none", label = "VV")
PyPlot.plot(xdata[arg],ydata_lc[arg], color = "pink", marker = ".", mfc="none", label = "VVc")
PyPlot.plot(xdata[arg],ydata_lc_P[arg], color = "cyan", marker = "v", mfc="none", label = "VVc impr")
for impr_set in IMPR_SET
    ydata_ll = [HVPQ_ll_NP[ens.id][impr_set][Qarg][1] for ens in ensInfo]
    ydata_lc = [HVPQ_lc_NP[ens.id][impr_set][Qarg][1] for ens in ensInfo]

    PyPlot.plot(xdata[arg],ydata_ll[arg], color = color_dict[impr_set]["ll"], marker = fmt_dict[impr_set]["ll"], mfc="none", label = "VV np impr; Set $impr_set")
    PyPlot.plot(xdata[arg],ydata_lc[arg], color = color_dict[impr_set]["lc"], marker = fmt_dict[impr_set]["lc"], mfc="none", label = "VVc np impr; Set $impr_set")
end
legend(loc="best")
xlabel(L"$a^2/8t_0$")
ylabel(latexstring("a_{\\mu}^{\\rm{hvp}}[\\rm{LO}^{CC}_{\\rm{SDsub}}](Q=$Q\\ \\rm{GeV})"))
tight_layout()
display(gcf())
close()



##

sortarg = sortperm(getfield.(ensInfo,:beta))
ensInfo[sortarg]

HVP3l_rat_ll = [HVP3l0./HVPQ_ll[ens.id] for ens in ensInfo[sortarg]]
HVP3l_rat_lc = [HVP3l0./HVPQ_lc[ens.id] for ens in ensInfo[sortarg]]

Q = 3.5

PyPlot.plot(getindex.(HVP3l_rat_ll,findfirst(x -> x == Q, Qlist)))
PyPlot.plot(getindex.(HVP3l_rat_lc,findfirst(x -> x == Q, Qlist)))
display(gcf())
close()

##

ens = EnsInfo("H102")
path_rw = !isempty(filter(x->occursin(ens.id, x), readdir(joinpath(path_rw_,"reweight_deflated"), join=true))) ? joinpath(path_rw_,"reweight_deflated") : path_rw_
corr = get_corr(path_HVP, ens, "light", "V1V1", path_rw=path_rw, frw_bcwd=false, L=1)
corr.obs
v1t10 = get_corr(path_HVP, ens, "light", "V1T10", path_rw=path_rw, frw_bcwd=false, L=1)
v1t10.obs

g33_ll, g33_lc = corr33(path_HVP, ens, path_rw=path_rw, impr=true, impr_set="1old", cons=true, frw_bcwd=true, std=STD_DERIV)

corr33tl_ll, corr33tl_lc = read_tree_level_v33(path_tl, cons=true)
corr33tl_v3s03_ll, corr33tl_v3s03_lc = read_tree_level_v3sig03(path_tl, cons=true, massless=true)

corr33tl_ll
corr33tl_v3s03_ll

cv_loc.(b_values)
cv_loc_old.(b_values)

##==========================> Tree-level computation <==========================##


diag = "LO"  # LO  NLOa  NLOb  NLOa&b
impr_set = "2"

Qlist = [3.5, 4.0, 5.0, 6.0, 7.0, 8.0]

STD_DERIV = false

C4 = Dict("LO" => π^2/9, "NLOa" => NaN, "NLOb" => NaN, "NLOa&b" => NaN) # to be added in constants

exp_diag = diag == "LO" ? 2 : 3

HVPQ_ll = Dict(); HVPQ_lc = Dict()
for ens in ensInfo[[ens ∉ ensNOcharm for ens in getfield.(ensInfo,:id)]]
    corr33tl_ll, corr33tl_lc = read_tree_level_v33(path_tl, cons=true)
    corr33tl_v3s03_ll, corr33tl_v3s03_lc = read_tree_level_v3sig03(path_tl, cons=true, massless=true)

    corr33tl_ll = uwreal.(corr33tl_ll)
    corr33tl_lc = uwreal.(corr33tl_lc)

    # beta = ens.beta
    # if impr_set == "1"
    #     cv_l = cv_loc(beta)
    #     cv_c = cv_cons(beta)
    # elseif impr_set =="1old"
    #     cv_l = cv_loc_old(beta)
    #     cv_c = cv_cons_old(beta)
    # elseif impr_set =="2"
    #     cv_l = cv_loc_set2(beta)
    #     cv_c = cv_cons_set2(beta)
    # end
    cv_l = 0.0
    cv_c = 0.0

    improve_corr_vkvk!(corr33tl_ll, corr33tl_v3s03_ll, 2*cv_l, std=STD_DERIV, treelevel=true)
    improve_corr_vkvk_cons!(corr33tl_lc, corr33tl_v3s03_ll, corr33tl_v3s03_lc, cv_l, cv_c, std=STD_DERIV, treelevel=true)

    TMR    = BDIOread_TMR(path_bdio,ens,diag,SU3=true)
    t0     = t0sym(ens.beta)
    factor = hbarc * sqrt(t0)/sqrtt0_ph

    sym_points = Int64(HVPobs.Data.get_T(ens.id)/2+1); t = collect(1:sym_points)
    aens = sqrtt0_ph / sqrt(t0); tfm = aens.*(t.-1)

    TMRw = TMR .* Window("SD")(tfm)
    TMRb(Q::Float64) = ((16/(Q/factor)^2)^2 * C4[diag] * (massmu/factor)^2) .* sin.((Q/factor/4) .* (t .- 1)).^4
    TMRsub(Q::Float64) = TMRw .- (Window("SD")(0) .* TMRb(Q))

    HVPQ_ll[ens.id] = []; HVPQ_lc[ens.id] = []
    for Q in Qlist
        int_ll = value.(corr33tl_ll[t]) .* value.(TMRsub(Q))
        int_lc = value.(corr33tl_lc[t]) .* value.(TMRsub(Q))
        amu_ll = (alpha/pi)^exp_diag * sum(int_ll) * 1e10
        amu_lc = (alpha/pi)^exp_diag * sum(int_lc) * 1e10

        push!(HVPQ_ll[ens.id],amu_ll/2)
        push!(HVPQ_lc[ens.id],amu_lc/2)
    end
end

HVP3l0 = compute_HVPtl0(diag,"SDsub",Qlist,path_coef)

sortarg = sortperm(getfield.(ensInfo[[ens ∉ ensNOcharm for ens in getfield.(ensInfo,:id)]],:beta))

HVP3l_rat_ll = [HVP3l0./HVPQ_ll[ens.id] for ens in ensInfo[[ens ∉ ensNOcharm for ens in getfield.(ensInfo,:id)]][sortarg]]
HVP3l_rat_lc = [HVP3l0./HVPQ_lc[ens.id] for ens in ensInfo[[ens ∉ ensNOcharm for ens in getfield.(ensInfo,:id)]][sortarg]]

##
Q = 3.5

PyPlot.plot(getindex.(HVP3l_rat_ll,findfirst(x -> x == Q, Qlist)))
PyPlot.plot(getindex.(HVP3l_rat_lc,findfirst(x -> x == Q, Qlist)))
display(gcf())
close()

##

list_ll = [0.757641,0.020128000091040113,0.00562739286362537,0.0018639937243439758,0.0008112608274905469,0.00041611394372048605,0.00023976834105176022]

ens = EnsInfo("H101")

corr33tl_ll, corr33tl_lc = read_tree_level_v33(path_tl, cons=true)
corr33tl_v3s03_ll, corr33tl_v3s03_lc = read_tree_level_v3sig03(path_tl, cons=true, massless=true)

corr33tl_ll = uwreal.(corr33tl_ll); corr33tl_lc = uwreal.(corr33tl_lc)

mycorr33tl_ll = corr33tl_ll[:]

beta = ens.beta
if impr_set == "1"
    cv_l = cv_loc(beta)
    cv_c = cv_cons(beta)
elseif impr_set =="1old"
    cv_l = cv_loc_old(beta)
    cv_c = cv_cons_old(beta)
elseif impr_set =="2"
    cv_l = cv_loc_set2(beta)
    cv_c = cv_cons_set2(beta)
end

improve_corr_vkvk!(mycorr33tl_ll, corr33tl_v3s03_ll, 2*cv_l, std=false, treelevel=true)
improve_corr_vkvk_cons!(corr33tl_lc, corr33tl_v3s03_ll, corr33tl_v3s03_lc, cv_l, cv_c, std=false, treelevel=true)

mycorr33tl_ll
corr33tl_ll
##

improve_corr_vkvk!(g3l_33_ll[k], g3l_v3s03_ll, 2*cv_l, std=false, treelevel=true)
improve_corr_vkvk_cons!(g3l_33_lc[k], g3l_v3s03_ll, g3l_v3s03_lc, cv_l, cv_c, std=false, treelevel=true)