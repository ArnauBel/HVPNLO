# Import packages

using Revise

include("../../HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs

using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2
using DelimitedFiles

using SpecialFunctions

using Plots
using PyPlot
using Colors

using ProgressBars
using Suppressor

# include uwreal constants

include("../../HVPtool/uwConst.jl")

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "..", "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "..", "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "..", "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_coef = joinpath(julia_script_directory, "..", "..", "..", "KernelCoeff")

path_spec = joinpath(julia_script_directory, "..", "..", "..", "HVPData", "spectroscopy")

pFVC_MLL  = joinpath(julia_script_directory, "..", "..", "..", "HVPData", "FSE_MLL")


# path_plot = joinpath(julia_script_directory, "..", "..", "..", "Slides & Plots","Plots")

# We do not have charm or disconnected data for some of the ensembles

ensNOcharm = ["J501","N451","D150","D451","J304","C102","D251","D201","J306","J307","F300","H200","N452"]
ensNOdisc  = ["D251","J306","J307","F300","D450"]

ensSPECdata = ["D200","E250"]  # J303

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Kernel piece by piece NLOa")

SAVE     = false
OVERSAVE = false

that = collect(range(0.0,10.0,10000))

that_small = collect(range(0.0,7.0,10000))
that_asymp = collect(range(2.0,10.0,10000))

contents4a = read(joinpath(path_coef, "NLO_diagram4a.txt"), String)
an, bn, cn, dn, anb11, anb12, anb21, anb22, anb23 = [Float64.(vec) for vec in HVPtool.TMR.extract_coef.(contents4a, ["an","bn","cn","dn","anb11","anb12","anb21","anb22","anb23"])]

idcsSmall = collect(4:2:Integer(length(an)*2+2))
factorials = factorial.(big.(idcsSmall))
idcsLarge = collect(0:length(anb11)-1)

tildef_small = []
for t in that_small
    push!(tildef_small, - sum((an .+  bn .* pi^2 .+ cn .* (custom_log(t)+GammaEuler) .+ dn .* ((custom_log(t)+GammaEuler)^2)) .* ((t .^ idcsSmall) ./ Float64.(factorials,RoundUp))))
end

tildef_asymp = []
for t in that_asymp
    push!(tildef_asymp, - (t^2 * (197/144 + pi^2/12 - pi^2*log(2)/2 + 3*zeta(3)/4)/2 - pi * t /8 + (log(t) + GammaEuler) * (1 - 5 / (12 * t^2)) + 653/216 - 127 * pi^2/144 - 7 * zeta(3)/4 + 7 * pi^2 * log(2)/6 -1.472467138/t + 1.1589872337/t^2))
end

fig = figure(figsize=(6,4.5))

PyPlot.plot(that_small, tildef_small, color = "blue" , linestyle="--", label=latexstring("t \\ll 1/m_\\mu"))
PyPlot.plot(that_asymp, tildef_asymp, color = "green", linestyle="--", label=latexstring("t \\gg 1/m_\\mu"))
# PyPlot.plot(that, -Tildef4a(that,path_coef), color = "orange", linestyle=":")

xlabel(latexstring("\\hat{t}"))
ylabel(latexstring("-\\frac{m_\\mu^2}{16\\pi^2}\\tilde{f}^{(\\rm{NLOa})}(\\hat{t})"))

legend()
tight_layout()
display(gcf())
if SAVE
    p = create_path(julia_script_directory,["KenrelNLOa_piecebypiece.pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p)
end
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Kernel piece by piece NLOb")

SAVE     = false
OVERSAVE = false

that = collect(range(0.0,20.0,10000))

that_small = collect(range(0.0,8.0,10000))
that_asymp = collect(range(8.0,20.0,10000))

contents4b = read(joinpath(path_coef, "NLO_diagram4b.txt"), String)
an, bn, cn, dn = HVPtool.TMR.extract_coef.(contents4b, ["an","bn","cn","dn"])

idcs = collect(4:2:Integer(length(an)*2+2))
factorials = factorial.(big.(idcs))

tildef_small = []
for t in that_small
    first_term = an .+  bn .* pi^2 .+ cn .* (custom_log(t)+GammaEuler) .+ dn .* ((custom_log(t)+GammaEuler)^2)
    second_term = (t .^ idcs) ./ Float64.(factorials,RoundUp)
    push!(tildef_small, sum(first_term .* second_term))
end

M = masse/massmu
tildef_asymp = []
for t in that_asymp
    push!(tildef_asymp, -((3*pi^2)/(64*M)) +
              (1/216) * (428 + 24*pi^2 - 75*t^2) +
              (1/32) * M * pi^2 * (-15 + 4*t^2) +
              (1/6) * M^2 * (5 + 9*t^2) -
              (5/192) * M^3 * pi^2 * (-7 + 24*t^2) +
              (M^6 * (911930975 + 102977616*t^2)) / 8467200 +
              (1/900) * M^4 * (-537 + 50*(44 + 3*pi^2)*t^2) +
              ( (13/9) - t^2/6 +
                M^4*(2/5 - (7*t^2)/3) +
                2*M^2*(-1 + t^2) +
                M^6*(5653/63 + (329*t^2)/40)
              ) * log(M) +
              (2/3 + M^4*t^2) * log(M)^2)
end

fig = figure(figsize=(6,4.5))

PyPlot.plot(that_small, tildef_small, color = "blue" , linestyle="--", label=latexstring("t \\ll 1/m_\\mu"))
PyPlot.plot(that_asymp, tildef_asymp, color = "green", linestyle="--", label=latexstring("t \\gg 1/m_e"))
PyPlot.plot(that, (massmu^2/(16π^2)).*Tildef4b_num.(that), color = "gray", linestyle=":", label="Exact")

xlabel(latexstring("\\hat{t}"))
ylabel(latexstring("\\frac{m_\\mu^2}{16\\pi^2}\\tilde{f}^{(\\rm{NLOb})}(\\hat{t})"))

legend()
tight_layout()
display(gcf())
if SAVE
    p = create_path(julia_script_directory,["KenrelNLOb_piecebypiece.pdf"],OVERWRITE=OVERSAVE)
    PyPlot.savefig(p)
end
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

@info("Pie-chart")
