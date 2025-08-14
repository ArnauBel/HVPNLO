# Export packages (including HVPobs)
using HVPobs
using PyPlot
using Plots
using ADerrors

##
# Export Corr33 function from "data_management.jl"
# Export TMR functions from "KernelTMRNLO.jl"
##

include("data_management.jl")
export corr33 

include("KernelTMRNLO.jl")
export Tildef4aInner, Tildef4a, Tildef4bInner, Tildef4b, Tildef4cInner, Tildef4c
export amuHVPNLO

##############################################################################################################
# Correlator export and test :
##
# Path definition and data extraction
##

julia_script_directory = @__DIR__

path_HVP  = "/Users/cesc/Desktop/Physics/JGU - PhD/LatticeData/HVP_data"
path_rw = "/Users/cesc/Desktop/Physics/JGU - PhD/LatticeData/rwf_deflated"

path_HVP  = joinpath(julia_script_directory, "..", "LatticeData", "HVP_data")
path_rw = joinpath(julia_script_directory, "..", "LatticeData", "rwf_deflated")
Ensamble = "N202"
ens = EnsInfo(Ensamble)

mycorr33locallocal, mycorr33localconserved = corr33(path_HVP, ens, path_rw = path_rw, frw_bcwd = true)
uwerr.(mycorr33locallocal.obs)
uwerr.(mycorr33localconserved.obs)

mycorr33locallocal.obs

##
# Plot local - local and local - conserved -> [in lattice units]
##

sym_points = Int64(length(mycorrlocallocal.obs)/2+1)
errorbar(collect(0:sym_points-1), -value.(mycorrlocallocal.obs[1:sym_points]), ADerrors.err.(mycorrlocallocal.obs[1:sym_points]), fmt="s", label="Local-local", color = "green", capsize=2)
errorbar(collect(0:sym_points-1), -value.(mycorrlocalconserved.obs[1:sym_points]), ADerrors.err.(mycorrlocalconserved.obs[1:sym_points]), fmt="s", label="Local-conserved", color = "red", capsize=2)
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title(L" $G^{(33)}(t)$ ; Ens = " * Ensamble)
xlabel("t [a]")
ylabel("G(t)")
ylim((0.5*minimum(filter(x -> x >= 0, -value.(mycorrlocallocal.obs))), 2*maximum(-value.(mycorrlocallocal.obs))))
yscale("log")
legend(["Local-local","Local-conserved"], loc  = "best")
grid("on")
display(gcf())      #display the figure
close()

##
# Plot local - local and local - conserved -> [in femtometers]
##

a_ens =  value(a(ens.beta))
tfm =  a_ens .* collect(0:sym_points-1)

errorbar(tfm, -value.(mycorrlocallocal.obs[1:sym_points]), ADerrors.err.(mycorrlocallocal.obs[1:sym_points]), fmt="s", label="Local-local", color = "green", capsize=2)
errorbar(tfm, -value.(mycorrlocalconserved.obs[1:sym_points]), ADerrors.err.(mycorrlocalconserved.obs[1:sym_points]), fmt="s", label="Local-conserved", color = "red", capsize=2)
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title(L" $G^{(33)}(t)$ ; Ens = " * Ensamble)
xlabel("t [fm]")
ylabel("G(t)")
ylim((0.5*minimum(filter(x -> x >= 0, -value.(mycorrlocallocal.obs))), 2*maximum(-value.(mycorrlocallocal.obs))))
yscale("log")
legend(["Local-local","Local-conserved"], loc  = "best")
grid("on")
display(gcf())      #display the figure
close()

##
# Plot local - local and local - conserved -> [in 1/GeV]
##

hbarc = 0.1973269804        # GeV * fm

tGeV = tfm ./ hbarc

errorbar(tGeV, -value.(mycorrlocallocal.obs[1:sym_points]), ADerrors.err.(mycorrlocallocal.obs[1:sym_points]), fmt="s", label="Local-local", color = "green", capsize=2)
errorbar(tGeV, -value.(mycorrlocalconserved.obs[1:sym_points]), ADerrors.err.(mycorrlocalconserved.obs[1:sym_points]), fmt="s", label="Local-conserved", color = "red", capsize=2)
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title(Ensamble*" correlator")
xlabel(L"$t$ [$GeV^{-1}$]")
ylabel("G(t)")
ylim((0.5*minimum(filter(x -> x >= 0, -value.(mycorrlocallocal.obs))), 2*maximum(-value.(mycorrlocallocal.obs))))
yscale("log")
legend(["Local-local","Local-conserved"], loc  = "best")
grid("on")
display(gcf())      #display the figure
close()

##############################################################################################################
# TMR functions test :
##

# Get the directory of the Julia script and then to the txt files

path_to_coef = joinpath(julia_script_directory, "..", "Coefficients")

tt = collect(range(0, stop=2, length=100))

##
# f[4a]
##

myresultsa = Tildef4a(t,path_to_coef)

figa, axa = subplots(figsize=(8, 6))

axa.plot(t, myresultsa, label="Tildef4a", color = "blue")
axis("tight")
xlabel(L"\hat{t}")
ylabel(L"\tilde{f}(t)")
title(L"$\tilde{f}_4^{(a)}(\hat{t})$")
legend([L"\tilde{f}_4^{(a)}(\hat{t})"], loc  = "best")
grid("on")
display(gcf())      #display the figure
close()

##
# f[4b]
##

myresultsb = Tildef4b(t,path_to_coef)

figb, axb = subplots(figsize=(8, 6))

axb.plot(t, myresultsb, label="Tildef4b", color = "green")
axis("tight")
xlabel(L"\hat{t}")
ylabel(L"\tilde{f}(t)")
title(L"$\tilde{f}_4^{(b)}(\hat{t})$")
legend([L"\tilde{f}_4^{(b)}(\hat{t})"], loc  = "best")
grid("on")
display(gcf())      #display the figure
close()

##
# f[4c]
##

myresultsc = Tildef4c(tt, path_to_coef)

figc, axc = subplots(figsize=(8, 8))
myheatmap = axc.pcolormesh(tt, tt, myresultsc, cmap="Reds")
colorbar(myheatmap, label="TMR value")
xlabel(L"\hat{t}")
ylabel(L"\hat{\tau}")
title(L"$\tilde{f}_2^{(c)}(\hat{t},\hat{\tau})$")
display(gcf())      #display the figure
close()

##
# f[4a] vs f[4b]
##

figab, axab = subplots(figsize=(8, 6))

axab.plot(t, -myresultsa, label=L"\tilde{f}_4^{(a)}(\hat{t})", color = "blue")
axab.plot(t, 2*myresultsb, label=L"\tilde{f}_4^{(b)}(\hat{t})", color = "green")
axis("tight")
xlabel(L"\hat{t}")
ylabel(L"\tilde{f}(t)")
#yscale("log")
title(L"$\tilde{f}_4^{(a,b)}(\hat{t})$")
legend([L"\tilde{f}_4^{(a)}(\hat{t})",L"\tilde{f}_4^{(b)}(\hat{t})"], loc  = "best")
grid("on")
display(gcf())      #display the figure
close()

##############################################################################################################
# Intrgrands and HVP NLO values:
##

## In 1/GeV

alpha = 1/137.035999084
hbarc = 0.1973269804        # GeV * fm

integrand_ll_a_GeV = - mycorrlocallocal.obs[1:sym_points] .* Tildef4a(massmu .* tGeV,path_to_coef)
integrand_lc_a_GeV = - mycorrlocalconserved.obs[1:sym_points] .* Tildef4a(massmu .* tGeV,path_to_coef)
uwerr.(integrand_ll_a_GeV)
uwerr.(integrand_lc_a_GeV)

errorbar(tGeV, value.(integrand_ll_a_GeV), ADerrors.err.(integrand_ll_a_GeV), fmt="s", label="Local-local integrand", color = "green", capsize=2)
errorbar(tGeV, value.(integrand_lc_a_GeV), ADerrors.err.(integrand_lc_a_GeV), fmt="s", label="Local-conserved integrand", color = "red", capsize=2)
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title("Integrand; Ensamble = " * Ensamble)
xlabel(L"$t$ [$GeV^{-1}$]")
ylabel(L"$G(t)$ $\tilde{f}_4^{(a)}(m_\mu t)$")
#ylim((0.5*minimum(filter(x -> x >= 0, -value.(mycorrlocallocal.obs))), 2*maximum(-value.(mycorrlocallocal.obs))))
#yscale("log")
legend(["Local-local integrand","Local-conserved integrand"], loc  = "best")
grid("on")
display(gcf())      #display the figure
close()

## "Manual" computation of the HVP (in 1/GeV)

(alpha/pi)^3 * (a_ens/hbarc) * sum(value.(integrand_ll_a_GeV)) * 10e9
(alpha/pi)^3 * (a_ens/hbarc) * sum(ADerrors.err.(integrand_ll_a_GeV)) * 10e9

(alpha/pi)^3 * (a_ens/hbarc) * sum(value.(integrand_lc_a_GeV)) * 10e9
(alpha/pi)^3 * (a_ens/hbarc) * sum(ADerrors.err.(integrand_lc_a_GeV)) * 10e9

#uwerr((alpha/pi)^3 * (a_ens/hbarc) * sum(integrand_ll_a) * 10e10)

## In fm

alpha = 1/137.035999084
hbarc = 0.1973269804        # GeV * fm

factor = massmu/hbarc

integrand_ll_a = - mycorrlocallocal.obs[1:sym_points] .* Tildef4a(factor .* tfm,path_to_coef)
integrand_lc_a = - mycorrlocalconserved.obs[1:sym_points] .* Tildef4a(factor .* tfm,path_to_coef)
uwerr.(integrand_ll_a)
uwerr.(integrand_lc_a)

errorbar(tfm, value.(integrand_ll_a), ADerrors.err.(integrand_ll_a), fmt="s", label="Local-local integrand", color = "green", capsize=2)
errorbar(tfm, value.(integrand_lc_a), ADerrors.err.(integrand_lc_a), fmt="s", label="Local-conserved integrand", color = "red", capsize=2)
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title("Integrand; Ensamble = " * Ensamble)
xlabel("t [fm]")
ylabel(L"$G(t)$ $\tilde{f}_4^{(a)}(m_\mu t)$")
#ylim((0.5*minimum(filter(x -> x >= 0, -value.(mycorrlocallocal.obs))), 2*maximum(-value.(mycorrlocallocal.obs))))
#yscale("log")
legend(["Local-local integrand","Local-conserved integrand"], loc  = "best")
grid("on")
display(gcf())      #display the figure
close()

## "Manual" computation of the HVP (in fm)

(alpha/pi)^3 * a_ens * sum(value.(integrand_ll_a)) * 10e9 / hbarc
(alpha/pi)^3 * a_ens * sum(ADerrors.err.(integrand_ll_a)) * 10e9 / hbarc

(alpha/pi)^3 * a_ens * sum(value.(integrand_lc_a)) * 10e9 / hbarc
(alpha/pi)^3 * a_ens * sum(ADerrors.err.(integrand_lc_a)) * 10e9 / hbarc

## In lattice units
t0ens = t0(ens.beta)
t0_ph = uwreal([0.1439, 0.006], "sqrtt0 [fm]") 
uwerr(t0_ph)

alpha = 1/137.035999084
hbarc = 0.1973269804        # GeV * fm

integrand_ll_a = - mycorrlocallocal.obs[1:sym_points] .* Tildef4a((massmu/hbarc * value(t0_ph)/sqrt(value(t0ens))) .* collect(0:sym_points-1),path_to_coef)
integrand_lc_a = - mycorrlocalconserved.obs[1:sym_points] .* Tildef4a((massmu/hbarc * value(t0_ph)/sqrt(value(t0ens)))  .* collect(0:sym_points-1),path_to_coef)
uwerr.(integrand_ll_a)
uwerr.(integrand_lc_a)

errorbar(collect(0:sym_points-1), value.(integrand_ll_a), ADerrors.err.(integrand_ll_a), fmt="s", label="Local-local integrand", color = "green", capsize=2)
errorbar(collect(0:sym_points-1), value.(integrand_lc_a), ADerrors.err.(integrand_lc_a), fmt="s", label="Local-conserved integrand", color = "red", capsize=2)
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title("Integrand; Ensamble = " * Ensamble)
xlabel("t [lattice units]")
ylabel(L"$G(t)$ $\tilde{f}_4^{(a)}(m_\mu t)$")
#ylim((0.5*minimum(filter(x -> x >= 0, -value.(mycorrlocallocal.obs))), 2*maximum(-value.(mycorrlocallocal.obs))))
#yscale("log")
legend(["Local-local integrand","Local-conserved integrand"], loc  = "best")
grid("on")
display(gcf())      #display the figure
close()

## "Manual" computation of the HVP (in lattice spacing)
t0_ph
(alpha/pi)^3 * (hbarc * sqrt(value(t0ens))/value(t0_ph))^2 * sum(value.(integrand_ll_a)) * 10e9
(alpha/pi)^3 * (hbarc * sqrt(value(t0ens))/value(t0_ph))^2 * sum(ADerrors.err.(integrand_ll_a)) * 10e9

(alpha/pi)^3 * (hbarc * sqrt(value(t0ens))/value(t0_ph))^2 * sum(value.(integrand_lc_a)) * 10e9
(alpha/pi)^3 * (hbarc * sqrt(value(t0ens))/value(t0_ph))^2 * sum(ADerrors.err.(integrand_lc_a)) * 10e9

sqrt(t0ens)/t0_ph

##############################################################################################################
# NLO HVP computation for diferent ensambles
##
# Ensambles with mpi = mk
##

#amuHVPNLO("a", EnsInfo("A653"), plot=true)     # key missing for all beta = 3.34

amuHVPNLO("a", EnsInfo("B450"), plot=true)
amuHVPNLO("b", EnsInfo("B450"), plot=true)

amuHVPNLO("a", EnsInfo("N202"), plot=true)
amuHVPNLO("b", EnsInfo("N202"), plot=true)

#amuHVPNLO("a", EnsInfo("N300"), plot=true)     # ensamble missing in HVP_data

amuHVPNLO("a", EnsInfo("J500"), plot=true)
amuHVPNLO("b", EnsInfo("J500"), plot=true)

##
# Convergense towards physical mass (a = 0.06426)
##

amuHVPNLO("a", EnsInfo("N202"), plot=true)
amuHVPNLO("b", EnsInfo("N202"), plot=true)

amuHVPNLO("a", EnsInfo("N203"), plot=true)
amuHVPNLO("b", EnsInfo("N203"), plot=true)

amuHVPNLO("a", EnsInfo("N200"), plot=true)
amuHVPNLO("b", EnsInfo("N200"), plot=true)

##############################################################################################################
# Renormalization of the correlators
##

Z3 = get_Z3(ens) 
uwerr(Z3)
Z8 = get_Z8(ens)
uwerr(Z8)
Z08 = get_Z08(ens)
uwerr(Z08)

mycorr88locallocal, mycorr88localconserved = corr88_conn(path_HVP, ens,  mycorr33locallocal, g33_lc = mycorr33localconserved, path_rw = path_rw)
mycorr08locallocal, mycorr08localconserved = corr08_conn(mycorrlocallocal, mycorr88locallocal; g33_lc = mycorrlocalconserved, g88_lc = mycorr88localconserved)

mycorr88locallocal.obs

mycorrRlocallocal = Z3^2 .* mycorrlocallocal.obs .+ (Z8^2/3) .* mycorr88locallocal.obs .+ (2*Z8*Z08/3) .* mycorr08locallocal.obs
uwerr.(mycorrRlocallocal)

mycorrRlocalconserved = Z3 .* mycorrlocalconserved.obs .+ (Z8/3) .* mycorr88localconserved.obs .+ (Z08/3) .* mycorr08localconserved.obs
uwerr.(mycorrRlocalconserved)

##

alpha = 1/137.035999084
hbarc = 0.1973269804        # GeV * fm

factor = massmu/hbarc

integrand_ll_R_a = - mycorrlocallocalR .* Tildef4a(factor .* tfm,path_to_coef)
integrand_lc_R_a = - mycorrlocalconservedR .* Tildef4a(factor .* tfm,path_to_coef)
uwerr.(integrand_ll_R_a)
uwerr.(integrand_lc_R_a)

errorbar(tfm, value.(integrand_ll_R_a), ADerrors.err.(integrand_ll_R_a), fmt="s", label="Local-local integrand", color = "green", capsize=2)
errorbar(tfm, value.(integrand_lc_R_a), ADerrors.err.(integrand_lc_R_a), fmt="s", label="Local-conserved integrand", color = "red", capsize=2)
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title("Integrand; Ensamble = " * Ensamble)
xlabel("t [fm]")
ylabel(L"$G(t)$ $\tilde{f}_4^{(a)}(m_\mu t)$")
#ylim((0.5*minimum(filter(x -> x >= 0, -value.(mycorrlocallocal.obs))), 2*maximum(-value.(mycorrlocallocal.obs))))
#yscale("log")
legend(["Local-local integrand","Local-conserved integrand"], loc  = "best")
grid("on")
display(gcf())      #display the figure
close()


##############################################################################################################
# Pheno Model
##

using DelimitedFiles

# Read the data from the file
data = readdlm("CorrelatorData.txt", '\t', header=true)

# Extract columns and convert to Float64
tPheno = Baseparse.(string.(data[1][:, 1][100:100:300000]))  # Extract and convert the first column
ModelPheno = Baseparse.(string.(data[1][:, 2][100:100:300000]))  # Extract and convert the second column
Tildef4 = Baseparse.(string.(data[1][:, 3][100:100:300000]))  # Extract the third column

Tildef4Pheno = Tildef4a(massmu .* Float64.(tPheno),path_to_coef)

(alpha/pi)^3 * tPheno[1] * sum(ModelPheno .* Tildef4Pheno) * 10e9 

(16 * pi^2/massmu^2) * Tildef4
Tildef4Pheno

##

tPheno[end]
plot(tPheno, ModelPheno)
plot(tGeV, -value.(mycorrlocallocalR))
plot(tGeV, -value.(mycorrlocalconservedR))
yscale("log")
display(gcf())      #display the figure
close()

##

##############################################################################################################
##

function extract_coef(contents::String, list_name::String)
    start_index = findfirst(x -> contains(x, "# $(list_name)"), split(contents, '\n'))
    if start_index === nothing
        error("List of coefficients $(an) has not been found")
    end

    start_index += 1
    end_index = findfirst(x -> isempty(x) || x[1] == '#', split(contents, '\n')[start_index:end])
    end_index = end_index === nothing ? length(contents) : start_index + end_index - 2

    raw_list = string.(split(contents, '\n')[start_index:end_index])
    
    return Baseparse.(string.(raw_list))
end

function extract_coef2(contents::String, list_name::String)
    start_index = findfirst(x -> contains(x, "# $(list_name)"), split(contents, '\n'))
    if start_index === nothing
        error("List of coefficients $(an) has not been found")
    end

    start_index += 1
    end_index = findfirst(x -> isempty(x) || x[1] == '#', split(contents, '\n')[start_index:end])
    end_index = end_index === nothing ? length(contents) : start_index + end_index - 2

    raw_list = string.(split(contents, '\n')[start_index:end_index])
    
    return Baseparse.(string.(raw_list))
end

using DelimitedFiles

@time read(joinpath(path_to_coef, "NLO_diagram4a.txt"), String)
@time contents = readdlm(joinpath(path_to_coef, "NLO_diagram4a.txt"), '\t', '\n', skipstart=2)
contents[2]