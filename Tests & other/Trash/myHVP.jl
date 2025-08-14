using HVPobs
using ADerrors

using Plots
using PyPlot

## Export correlator, fv corr. and other usefull functions from "data_management.jl"
## Export TMR functions from "KernelTMRNLO.jl"
## Export HVP and FV from "amuNLO.jl"

include("data_management.jl")
export get_Z3, get_Z8, get_Z08, corr33, corr88_conn, corr08_conn, corrR  

include("TMRKernel.jl")
export Tildef4aInner, Tildef4a, Tildef4bInner, Tildef4b, Tildef4cInner, Tildef4c

include("amuNLO.jl")
export amuHVPNLO, amu∆G

################## Correlator export and test ##################

## Path definition and data extraction

julia_script_directory = @__DIR__

path_HVP  = joinpath(julia_script_directory, "..", "LatticeData", "HVP_data")
path_rw = joinpath(julia_script_directory, "..", "LatticeData", "rwf_deflated")
Ensamble = "N202"
ens = EnsInfo(Ensamble)

gR_ll, gR_lc = corrR(path_HVP, ens, path_rw = path_rw, frw_bcwd = true)
uwerr.(gR_ll.obs)
uwerr.(gR_lc.obs)

## Scale seting t0

t0_ph = uwreal([0.1439, 0.0006], "sqrtt0 [fm]") 

path_ms = joinpath(julia_script_directory, "..", "LatticeData", "ms_t0_dat")

t0ens = get_t0(path_ms, ens, path_rw = path_rw, pl = true)

aens = t0_ph ./ sqrt.( t0ens )
uwerr.(aens)
aens

## Plot for the renormalized correlator

aens =  value(a(ens.beta))
sym_points = Int64(length(gR_ll.obs)/2+1)
tfm =  aens .* collect(0:sym_points-1)

errorbar(tfm, -value.(gR_ll.obs[1:sym_points]), ADerrors.err.(gR_ll.obs[1:sym_points]), fmt="s", label="Local-local", color = "green", capsize=2)
errorbar(tfm, -value.(gR_lc.obs[1:sym_points]), ADerrors.err.(gR_lc.obs[1:sym_points]), fmt="s", label="Local-conserved", color = "red", capsize=2)
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title(L" $G^{(R)}(t)$; Ens = " * Ensamble)
xlabel("t [fm]")
#ylabel("G(t)")
ylim((0.5*minimum(filter(x -> x >= 0, -value.(gR_ll.obs))), 2*maximum(-value.(gR_lc.obs))))
yscale("log")
legend(["Local-local","Local-conserved"], loc  = "best")
PyPlot.grid("on")
display(gcf())      #display the figure
close()

## TMR Kernel plot

path_coef = joinpath(julia_script_directory, "..", "Coefficients")

massmu = 0.10565837
hbarc = 0.1973269804        # GeV * fm
factor = hbarc * sqrt(t0ens)/t0_ph     # convertion from GeV to 1/a units

TMRa = factor^2 .* Tildef4a((massmu/factor) .* collect(0:sym_points-1),path_coef)
uwerr.(TMRa)
TMRb = factor^2 .* Tildef4b((massmu/factor) .* collect(0:sym_points-1),path_coef)
uwerr.(TMRb)

errorbar(tfm, -value.(TMRa), ADerrors.err.(TMRa), fmt="s", label=L"$-\tilde{f}_4^{(a)}(t)$", color = "blue", capsize=2)
errorbar(tfm, value.(TMRb), ADerrors.err.(TMRb), fmt="s", label=L"$\tilde{f}_4^{(b)}(t)$", color = "orange", capsize=2)
axis("tight")
xlabel("t [fm]")
#ylabel(L"\tilde{f}_4(t)")
title("TMR Kernels")
legend([L"$-\tilde{f}_4^{(a)}(t)$",L"$\tilde{f}_4^{(b)}(t)$"], loc  = "best")
PyPlot.grid("on")
display(gcf())      #display the figure
close()

TMRc = factor^4 .* Tildef4c((massmu/factor) .* collect(0:sym_points-1), path_coef)
uwerr.(TMRc)

figc, axc = subplots(figsize=(8, 8))
myheatmap = axc.pcolormesh(tfm, tfm, value.(TMRc), cmap="Greens")  # Set colormap to Reds for red color
colorbar(myheatmap, label="TMR value")
xlabel("t [fm]")
ylabel(L"$\tau$ [fm]")
title(L"$\tilde{f}_4^{(c)}(t,\tau)$")
display(gcf())      #display the figure
close()

figc, axc = subplots(figsize=(8, 8))
myheatmap = axc.pcolormesh(tfm, tfm, ADerrors.err.(TMRc), cmap="Reds")
colorbar(myheatmap, label="TMR uncertainty")
xlabel("t [fm]")
ylabel(L"$\tau$ [fm]")
title(L"$\delta\tilde{f}_4^{(c)}(t,\tau)$")
display(gcf())      #display the figure
close()


################## NLO HVP computation ##################

## Total HVP

amuHVPNLO("a",ens, pl=true, errmult=10)
amuHVPNLO("b",ens, pl=true, errmult=10)
amuHVPNLO("c",ens, pl=true)

## HVP splitted by contributions

LocalLocal, LocalConserved = corrR(path_HVP, ens, path_rw = path_rw, frw_bcwd = true, split = true)

vec_myHVPNLO, vec_myIntegrand = amuHVPNLO("a", LocalLocal, t0ens, pl=true, int=true)
vec_myHVPNLO, vec_myIntegrand = amuHVPNLO("b", LocalLocal, t0ens, pl=true, int=true)
vec_myHVPNLO, vec_myIntegrand = amuHVPNLO("c", LocalLocal, t0ens, pl=true, int=true)

vec_myHVPNLO, vec_myIntegrand = amuHVPNLO("a", LocalConserved, t0ens, pl=true, int=true)
vec_myHVPNLO, vec_myIntegrand = amuHVPNLO("b", LocalConserved, t0ens, pl=true, int=true)
vec_myHVPNLO, vec_myIntegrand = amuHVPNLO("c", LocalConserved, t0ens, pl=true, int=true)

vec_myHVPNLOa, vec_myIntegranda = amuHVPNLO("a", LocalLocal, TMRa, t0ens=t0ens, pl=true, int=true)
vec_myHVPNLOb, vec_myIntegrandb = amuHVPNLO("b", LocalLocal, TMRb, t0ens=t0ens, pl=true, int=true)
vec_myHVPNLOc, vec_myIntegrandc = amuHVPNLO("c", LocalLocal, TMRc, t0ens=t0ens, pl=true, int=true)

## Finite Volume corrections

path_fvc = joinpath(julia_script_directory, "..", "LatticeData", "JKMPI")

fvc_HP = ∆GHP(path_fvc, ens, nmin=1, nmax=-1)

fvc_myHVPNLOa, fvc_myIntegranda = amu∆G("a", fvc_HP, t0ens, corr=LocalLocal[end], pl=true, int=true)
fvc_myHVPNLOb, fvc_myIntegrandb = amu∆G("b", fvc_HP, t0ens, corr=LocalLocal[end], pl=true, int=true)
fvc_myHVPNLOc, fvc_myIntegrandc = amu∆G("c", fvc_HP, t0ens, corr=LocalLocal[end], pl=true, int=true)

fvc_myHVPNLOa, fvc_myIntegranda = amu∆G("a", fvc_HP, TMRa, corr=LocalLocal[end], t0ens = t0ens, pl=true, int=true)
fvc_myHVPNLOb, fvc_myIntegrandb = amu∆G("b", fvc_HP, TMRb, corr=LocalLocal[end], t0ens = t0ens, pl=true, int=true)
fvc_myHVPNLOc, fvc_myIntegrandc = amu∆G("c", fvc_HP, TMRc, corr=LocalLocal[end], t0ens = t0ens, pl=true, int=true)

################## Bounding method ##################

## Upper and improved lowerr bound

t = collect(0:sym_points-1)
tcut = 20

E_eff=Eeff(tcut, LocalLocal[end].obs)
mpi = m_ens[LocalLocal[end].id]["m_pi"] * (1e-3/factor) #* (t0_ph / sqrt(t0ens))

upperbound = GUB(t, tcut, LocalLocal[end], mpi=mpi)
uwerr.(upperbound)
lowerbound = GUB(t, tcut, LocalLocal[end], Eeff=E_eff)
uwerr.(lowerbound)

errorbar(t, -value.(LocalLocal[end].obs[1:sym_points]), ADerrors.err.(LocalLocal[end].obs[1:sym_points]), fmt="s", label="Local-local", color = "green", capsize=2)
errorbar(t[tcut+1:end], value.(lowerbound), ADerrors.err.(lowerbound), fmt="s", label="Lower-bound", color = "gray", capsize=2)
errorbar(t[tcut+1:end], value.(upperbound), ADerrors.err.(upperbound), fmt="s", label="Upper-bound", color = "black", capsize=2)
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title(L" $G^{(R)}(t)$; Ens = " * Ensamble)
xlabel("t/a")
#ylabel("G(t)")
ylim((0.5*minimum(filter(x -> x >= 0, -value.(LocalLocal[end].obs))), 2*maximum(-value.(LocalLocal[end].obs))))
yscale("log")
legend(["Local-local","Lower-bound","Upper-bound"], loc  = "best")
PyPlot.grid("on")
display(gcf())      #display the figure
close()

## Application of the Bounding Method

amuNLO_a, lb_a, ub_a = amu_BM("a", LocalLocal[end], TMRa, mpi, t_step=1, t0ens=t0ens, pl=true)
amuNLO_b, lb_b, ub_b = amu_BM("b", LocalLocal[end], TMRb, mpi, t_step=1, t0ens=t0ens, pl=true)

amuNLO_b, lb_b, lb_impr, ub_b = amu_BM("b", LocalLocal[end], TMRb, mpi, impr=true, t_step=1, t0ens=t0ens, pl=true)












##
####################################
####################################
##













using SpecialFunctions
using QuadGK



m_pi = 412
L = 3.1
n = 1
yArray = collect(1:1e-3:20-1e-3)
integrand = besselk.(0, (mpi * sqrt(L^2 * n^2 + 4 * t^2)) .* yArray) .* sinh.((L * mpi * n) .* (yArray .- 1))
integrand

sum(integrand)
function sumnInner(t::Float64, mpi::Float64, L::Float64, n::Int64)

    term1 = besselk(2, mpi * sqrt(L^2 * n^2 + 4 * t^2)) / (mpi^2 * (L^2 * n^2 + 4 * t^2))

    yArray = collect(1:1e-3:20-1e-3)
    integrand = besselk(0, mpi * sqrt(L^2 * n^2 + 4 * t^2) .* yArray) * sinh(L * mpi * n .* (yArray .- 1))
    term2 = - 1 / (L * mpi * n) * 1e-3 * sum(integrand)

    return (mpi^4 * t) / (3 * π^2) * (term1 + term2)
end

factor = hbarc * sqrt(t0ens)/t0_ph     # convertion from GeV to 1/a units
n = 1
mpi = m_pi * factor * 1e-3
L = 3.1/aens[1]
quadgk(y -> besselk(0, mpi * sqrt(L^2 * n^2 + 4 * t^2) * y) * sinh(L * mpi * n * (y - 1)), 1, 10)

t = 10

intfunc(y, p) = besselk(2, mpi * sqrt(p[1]^2 * n^2 + 4 * t^2) * y) * sinh(L * p[2] * n * (y - 1))

p = [L,mpi]

int_error(intfunc, 1, 4, p)


##

using ADerrors

# Average values and covariance from
# https://inspirehep.net/literature/1477411
# here we try to reproduce the result
# Scale ratio = 21.86(42) Eq.5.6
avg = [16.26, 0.12, -0.0038]
Mcov = [ 0.478071  -0.176116   0.0135305
        -0.176116   0.0696489 -0.00554431
         0.0135305 -0.00554431 0.000454180]

p = cobs(avg, Mcov, "Beta function fit parameters")
g1s = uwreal([2.6723, 0.0064], 4)
g2s = 11.31

fint(x, p) = - (p[1] + p[2]*x^2 + p[3]*x^4)/x^3

int_error(fint, 1, 10, p)