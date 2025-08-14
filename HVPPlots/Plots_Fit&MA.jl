# Import packages

using Revise

include("../HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs

using ALPHAio, ADerrors
import ADerrors: err

using Statistics
using Distributions

using BDIO
using JLD2

using Plots
using PyPlot
using Colors

using ProgressBars
using Suppressor


# include uwreal constants

include("../HVPtool/uwConst.jl")

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_plot = joinpath(julia_script_directory,"..","..","Slides & Plots","Plots")

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

# Set plot parameters

wind = "ID"  #  NW  SD  SDsub  ID  LD  ILD
diag = "NLOa&b"  #  LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)
comp = "g33"  #  g33  g88  gCCconn  gCCdisc  gC8disc  g3333  g8888  gCCCC  g3388  g33CC  g88CC  ∆ls_amu  ∆lc_b

Q = 5.0  # virtuality for SDsub

# model_var_list = Function[a3,a2phi2,phi2sqr,phi2log,phi2inv,logphi2]
model_var_list = Function[a3,a2y,ysqr,ylog,yinv,logy]
MultFunc = nothing  # nothing  deltaphi

readIMPR_SET = ["1","2"] # ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

BLIND = false

FitCUT = "None"  # "None"  "beta"  "mass"  "beta&mass"

STD_DERIV  = false
tl_IMPR    = false
VREF       = true
RESC       = true

SimpleBase = false

path_bdio = path_bdio_dict["local"]


# Data reading and definitions

if wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]
    @info("Reading obsBDIO data:\n - wind: $wind\n - diag: $diag\n - comp: $comp\n - Fit type: $FitCUT\n - IMPR SET: $readIMPR_SET\n - Q: $Q")
else
    @info("Reading obsBDIO data:\n - wind: $wind\n - diag: $diag\n - comp: $comp\n - Fit type: $FitCUT\n - IMPR SET: $readIMPR_SET")
end


DictComptoKey = Dict{String,Vector{String}}(
    "g33"      => ["g33_ll","g33_lc"],
    "g88"      => ["g88_ll","g88_lc"],
    # only interested in the lc (local-conserved) discr. for the cc conn
    "gCCconn"  => ["gCCconn_SU3_lc"], # ["gCCconn_ll","gCCconn_lc"]
    # "gCCconn"  => ["gCCconn_SU3_ll","gCCconn_SU3_lc"], # ["gCCconn_SU3_ll","gCCconn_SU3_lc"]

    "∆ls_amu" => ["∆ls_amu_ll","∆ls_amu_lc"],
    "∆lc_b"   => ["∆lc_b_ll","∆lc_b_lc"],

    # only interested in the cc (conserved-conserved) discr. for the cc disc  & c8 disc
    "gCCdisc" => ["gCCdisc_cc"],
    "gC8disc" => ["gC8disc_cc"],

    "g3333"    => ["g3333_ll","g3333_lc"],
    "g8888"    => ["g8888_ll","g8888_lc"],
    "gCCCC"    => ["gCCCC_ll","gCCCC_lc"],
    "g3388"    => ["g3388_ll","g3388_lc"],
    "g33CC"    => ["g33CC_ll","g33CC_lc"],
    "g88CC"    => ["g88CC_ll","g88CC_lc"]
)

mykeys = DictComptoKey[comp]
key_len    = length(mykeys)

xdata = Dict(); ydata = Dict()
fitres = Dict()
res =  Dict(); param = Dict()
modelinfo = Dict()
MA = Dict(); info = Dict()
for impr_set in readIMPR_SET
    println("- Reading fit & MA for impr. set $impr_set")

    println("   - Reading X & Y data..")

    xdata, ydata[impr_set] = BDIOread_XYdata(path_bdio,diag,wind,comp,model_var_list,FitCUT,impr_set,resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

    println("   - Reading FitRes...")

    fitres[impr_set] = JDL2read_FitRes(path_bdio,diag,wind,comp,model_var_list,FitCUT,impr_set,resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

    println("   - Reading results and parameters...")

    res[impr_set], param[impr_set] = BDIOread_res(path_bdio,diag,wind,comp,model_var_list,FitCUT,impr_set,resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,Q=Q,BLIND=BLIND,param=true)

    println("   - Reading model information...")

    modelinfo = JDL2read_ModelInfo(path_bdio,diag,wind,comp,model_var_list,FitCUT,impr_set,resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)

    println("- Reading MA...")

    MA[impr_set], info[impr_set] = BDIOread_MA(path_bdio,diag,wind,comp,model_var_list,FitCUT,impr_set,resc=RESC,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=STD_DERIV,tlImpr=tl_IMPR,Vref=VREF,BLIND=BLIND,Q=Q)
end

xydata = Dict("xdata" => xdata, "ydata" => ydata)

model_len = modelinfo["length"]
a2RESCAL  = modelinfo["a2Rescaling"]
nens      = modelinfo["nens"]
ensInfo   = EnsInfo.(modelinfo["ensList"])

f_tot_isov, n_par_tot_isov, label_tot_isov = call_models(model_var_list,ensInfo,4,SimpleBase=SimpleBase,a2resc=a2RESCAL,MultFunc=MultFunc,fPiresc=RESC,na_max=modelinfo["na_max"],nmPi_max=modelinfo["nmPi_max"],nmK_max=modelinfo["nmK_max"])


@info("All data ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##


@info("Continiuum 'magic' plot")

IMPR_SET = readIMPR_SET  #  readIMPR_SET  ["1"]  ["2"]  ["1","2"]  ["1old","2"]  ["1","1old","2"]

npoints = 50

PROJpoints = true

println("   - Computing 'x' and 'y' plot points \n       Projected points ..... $PROJpoints")

noRESCALarg = !RESCAL ? collect(1:model_len) : collect(1:2:model_len)

xarr = [Float64.(range(1e-5, 1.5*maximum(value.(xydata["xdata"][:,1])), length=npoints)) fill(value(phi2_ph), npoints) fill(value(phi4_ph), npoints)]; yarr = Dict()
if PROJpoints
    xproj = deepcopy(xydata["xdata"])
    xproj[:,2] = fill(phi2_ph, length(xydata["xdata"][:,1]))::Vector{uwreal}
    xproj[:,3] = fill(phi4_ph, length(xydata["xdata"][:,1]))::Vector{uwreal}
end
yproj = Dict(); yproj_syst = Dict()
warg = Dict()
for (j,impr_set) in enumerate(IMPR_SET)

    yarr[impr_set] = Dict()
    yproj[impr_set] = Dict(); yproj_syst[impr_set] = Dict()
    warg[impr_set] = Dict()
    for (k,key) in enumerate(mykeys)

        yarr[impr_set][key] = Vector{Vector{uwreal}}()
        println("- Impr. set = $impr_set; Comp = $key:")
        println("   - [Computing bands...]")

        valid_indices = filter(i -> label_tot_isov[i][1] ∉ ["baseResc","baseSimpResc"], eachindex(info[impr_set]["weight"][key]))
        w = info[impr_set]["weight"][key][valid_indices]
        w = RESCAL ? w./sum(w) : w

        warg[impr_set][key] = sortperm(info[impr_set]["weight"][key][valid_indices]; rev=true) .|> i -> valid_indices[i]
        warg[impr_set][key] = warg[impr_set][key][info[impr_set]["weight"][key][warg[impr_set][key]] .> info[impr_set]["weight"][key][warg[impr_set][key][1]]/10000]

        for i in ProgressBar(warg[impr_set][key])
            my = f_tot_isov[i](xarr, param[impr_set][key][i]); uwerr.(my)
            push!(yarr[impr_set][key], my)
        end
        if PROJpoints
            println("   - [Computing projection points...]")
            yproj_list = Vector{Vector{uwreal}}()
            # w_ = Vector{Float64}()
            for i in ProgressBar(valid_indices)
                push!(yproj_list, f_tot_isov[i](xproj,value.(param[impr_set][key][i])))
            end
            
            yproj[impr_set][key] = Vector{uwreal}()
            yproj_syst[impr_set][key] = Vector{Float64}()
            for arg=1:nens
                y_ = [element[arg] for element in yproj_list]
                val, syst = model_average(y_,w)
                push!(yproj[impr_set][key],val[1])
                push!(yproj_syst[impr_set][key],syst)
            end
            uwerr.(yproj[impr_set][key])
        end
    end
end

#-- Stop for only plot compilation

# IMPR_SET = readIMPR_SET

# mykeys = [""]

wPen = 0.2

SHOWRES = false

# color_list = comp ∉ ["cc conn","cc disc","c8 disc"] ? ["blue","red","green","brown"] : ["blue","green"]
color_dict = Dict(
    ""     => Dict("cc" => "gold"),
    "1"    => Dict("ll" => "blue",  "lc" => "red"  ),
    "2"    => Dict("ll" => "green", "lc" => "brown")
); color_dict["1old"] = color_dict["1"]
fmt_list = ["^","o","s","d"]

factor = diag=="LO" ? 1 : 10

fig = figure(figsize=(10,7.5))
for (j,impr_set) in enumerate(IMPR_SET)
    for (k,key) in enumerate(mykeys)
        for i=1:length(yarr[impr_set][key])
            fill_between(xarr[:,1], value.(yarr[impr_set][key][i]).-err.(yarr[impr_set][key][i]), value.(yarr[impr_set][key][i]).+err.(yarr[impr_set][key][i]), alpha=(!RESCAL ? 1 : 2)*info[impr_set]["weight"][key][warg[impr_set][key][i]]*wPen, color=color_dict[impr_set]["$(key[end-1:end])"])
        end
        if PROJpoints
            errorbar(value.(xproj[:,1]), value.(yproj[impr_set][key]), err.(yproj[impr_set][key]), fmt=fmt_list[key_len*(j-1)+k], capsize=2, color=color_dict[impr_set]["$(key[end-1:end])"], mfc="none")
            errorbar(value.(xproj[:,1]), value.(yproj[impr_set][key]), sqrt.(err.(yproj[impr_set][key]).^2 .+ yproj_syst[impr_set][key].^2), fmt=fmt_list[key_len*(j-1)+k], capsize=2, color=color_dict[impr_set]["$(key[end-1:end])"], mfc="none", label="discr. $(key[end-1:end]); set $impr_set")
        end
    end
end
# compute final estimation from combining IMPR_SET and KEYS
weight = vcat([vcat([info[impr_set]["weight"][key] for key in mykeys]...) for impr_set in IMPR_SET]...)./(length(mykeys)*length(IMPR_SET))
res_tot = vcat([vcat([res[impr_set][key] for key in mykeys]...) for impr_set in IMPR_SET]...)
amu, syst = model_average(res_tot, weight); uwerr.(amu)

errorbar([0.0],value.(amu),err.(amu),fmt="o",mfc="none",color="black", ms=5, capsize=3)
errorbar([0.0],value.(amu),sqrt.(err.(amu).^2 .+ (syst)^2),fmt="o",color="black", ms=5, capsize=3)
digits = comp in ["gCCdisc","gC8disc"] ? 7 : 3
res_str = SHOWRES ? "  [$(round(value(amu[1]),digits=digits))($(round(err(amu[1]),digits=digits)))($(round(syst,digits=digits)))]" : ""
PyPlot.title("Projection to continuum extrapolation$res_str")
axvline(ls="dashed", color="black", lw=0.2, alpha=0.7)
xlabel(L"$a^2/8t_0$")
if diag != "LO" # we multiply the y axis by a factor 10
    formatter(x, pos) = string(round(10 * x, digits=2))  # Round to 2 decimal places
    ax = gca()
    ax.yaxis.set_major_formatter(PyPlot.matplotlib.ticker.FuncFormatter(formatter))
end
xlim(right=0.065)
# ylim(bottom=26.5/10)
# ylim([1.8/10,2.9/10])
diag_str = diag == "LO" ? "\\rm{LO}" : (diag == "NLOa&b" ? "\\rm{NLO}_{\\rm{a}\\&\\rm{b}}" : "\\rm{NLO}_{\\rm{$(diag[end])}}")
comp_str = comp[1] == '∆' ? (comp[end] != 'b' ? "\\Delta_{ls}(a_{\\mu})" : "\\Delta_{lc}(b)") : "\\mathrm{$(comp[2]),$(comp[3])}"
Q_str    = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "($Q\\ \\rm{GeV})" : ""
fact_str = diag == "LO" ? "\\times10^{10}" : "\\times10^{11}"
if wind == "NW"
    ylabel(latexstring("a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]$Q_str$fact_str"))
else
    ylabel(latexstring("\\left(a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]\\right)^{\\rm{$wind}}$Q_str$fact_str"))
end
PROJpoints ? legend(loc="lower center") : nothing
tight_layout()
display(gcf())
close()

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##


@info("Chiral projection plot")

SAVE     = true
OVERSAVE = false

ARGPLOT = 1  #  set to 1 for best plot

IMPR_SET = readIMPR_SET  #  readIMPR_SET  ["1"]  ["2"]  ["1","2"]  ["1old","2"]  ["1","1old","2"]

color = ["orange","red","purple","blue","green","brown"]
fmt = ["o","^","s","d","h","v"]

Pi_ph = !RESC ? phi2_ph : y_ph
K_ph  = !RESC ? phi4_ph : z_ph

for (i,impr_set) in enumerate(IMPR_SET)
    for (j,key) in enumerate(mykeys)
        label = []
        # firstarg, lastarg = [(key_len*(i-1)+j-1)*model_len+1,(key_len*(i-1)+j)*model_len]
        argSort = sortperm(eachindex(info[impr_set]["weight"][key]), 
                   by=i -> label_tot_isov[i][1] in ["baseResc", "baseSimpResc"] ? -Inf : info[impr_set]["weight"][key][i], 
                   rev=true)
        argPlot = argSort[ARGPLOT]
        discr   = key[end-1:end]
        fit     = fitres[impr_set][key][argPlot]
        model   = label_tot_isov[argPlot]
        myres   = res[impr_set][key][argPlot]; uwerr(myres)
        myparam = param[impr_set][key][argPlot]; uwerr.(myparam)

        xproj_ch = deepcopy(xydata["xdata"])
        xproj_ch[:,3] = fill(K_ph, length(xydata["xdata"][:,1]))::Vector{uwreal}
        ydata = f_tot_isov[argPlot](xproj_ch,value.(myparam)); uwerr.(ydata)

        println("Chosen model (diag=$diag, set=$impr_set, discr=$discr): $(label_tot_isov[argPlot]) (n: $argPlot)")
        println("- chi2/dof     = $(fit.chi2/fit.dof)")
        println("- chi2/chi2exp = $(fit.chi2/fit.chi2exp)")

        for (k,b) in  enumerate(sort(unique(getfield.(ensInfo, :beta))))
            push!(label,L"$\beta=$"*"$b")
            n_ = findall(x->x.beta == b, ensInfo)
            a2_aux = mean(value.(xydata["xdata"][:,1][n_]))
            errorbar(value.(xydata["xdata"][n_,2]), value.(ydata[n_]), err.(ydata[n_]), fmt=fmt[k], capsize=2, color=color[k], mfc="none", label=label[k])
            xxx = !RESC ? [fill(a2_aux,100) Float64.(range(0.02, 0.8, length=100)) fill(value(phi4_ph), 100)] : [fill(a2_aux,100) Float64.(range(0.02, 0.4, length=100)) fill(value(z_ph), 100)]
            yyy = f_tot_isov[argPlot](xxx, myparam)
            PyPlot.plot(xxx[:,2],value.(yyy),ls="--",color=color[k],lw=0.5)
        end
        xxx_ph = !RESC ? [fill(0.0,100) Float64.(range(0.02, 0.8, length=100)) fill(value(phi4_ph), 100)] : [fill(0.0,100) Float64.(range(0.02, 0.4, length=100)) fill(value(z_ph), 100)]
        yyy_ph = f_tot_isov[argPlot](xxx_ph, myparam); uwerr.(yyy_ph)
        
        errorbar(value(Pi_ph), value(myres), err(myres), fmt="^", capsize=2, color="black",label="ph.")
        push!(label,"ph.")
        PyPlot.plot(xxx_ph[:,2],value.(yyy_ph),ls="--",color="gray",lw=0.5)
        fill_between(xxx_ph[:,2], value.(yyy_ph)-err.(yyy_ph), value.(yyy_ph)+err.(yyy_ph), alpha=0.2, color="gray")
        PyPlot.title("Chiral and continuum extrapolation (Discr. $(key[end-1:end]); impr. set $impr_set)")
        axvline(x=value(Pi_ph), ls="dashed", color="black", lw=0.2, alpha=0.7) 
        if !RESC
            xlabel(L"$\Phi_2$")
        else
            xlabel(L"$y$")
        end
        if diag != "LO" # we multiply the y axis by a factor 10
            formatter(x, pos) = string(round(10 * x, digits=2))  # Round to 2 decimal places
            ax = gca()
            ax.yaxis.set_major_formatter(PyPlot.matplotlib.ticker.FuncFormatter(formatter))
        end
        # ylim([3,12.5])
        diag_str = diag == "LO" ? "\\rm{LO}" : (diag == "NLOa&b" ? "\\rm{NLO}_{\\rm{a}\\&\\rm{b}}" : "\\rm{NLO}_{\\rm{$(diag[end])}}")
        comp_str = comp[1] == '∆' ? (comp[end] != 'b' ? "\\Delta_{ls}(a_{\\mu})" : "\\Delta_{lc}(b)") : "\\mathrm{$(comp[2]),$(comp[3])}"
        Q_str    = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "($Q\\ GeV)" : ""
        fact_str = diag == "LO" ? "\\times10^{10}" : "\\times10^{11}"
        if wind == "NW"
            ylabel(latexstring("a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]$Q_str$fact_str"))
        else
            ylabel(latexstring("\\left(a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]\\right)^{\\rm{$wind}}$Q_str$fact_str"))
        end
        legend(label,loc="best")
        tight_layout()
        display(gcf())
        if SAVE
            RESCstr = !RESC ? "" : "_resc"
            p = create_path(path_plot,[diag,wind,"ChExtr_$(key)_$(impr_set)$(RESCstr).pdf"],OVERWRITE=OVERSAVE)
            PyPlot.savefig(p)
        end
        close()
    end
end

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##


@info("Hit plot")

IMPR_SET = readIMPR_SET #  readIMPR_SET  ["1"]  ["2"]  ["1","2"]  ["1old","2"]  ["1","1old","2"]

ACC = false

color_dict = Dict(
    ""     => Dict("cc" => "cyan"),
    "1"    => Dict("ll" => "blue",  "lc" => "red"  ),
    "2"    => Dict("ll" => "green", "lc" => "brown")
); color_dict["1old"] = color_dict["1"]


weight = vcat([vcat([info[impr_set]["weight"][key] for key in mykeys]...) for impr_set in IMPR_SET]...)./(length(mykeys)*length(IMPR_SET))
res_tot = vcat([vcat([res[impr_set][key] for key in mykeys]...) for impr_set in IMPR_SET]...); uwerr.(res_tot)

amu, syst = model_average(res_tot, weight); uwerr.(amu)

gca().axvspan(value(amu[1])-err(amu[1]), value(amu[1])+err(amu[1]), color="limegreen", alpha=0.3)
gca().axvspan(value(amu[1])-sqrt(err(amu[1])^2+syst^2), value(amu[1])+sqrt(err(amu[1])^2+syst^2), color="limegreen", alpha=0.3, label="Quoted result")

xmin = value(amu[1])-5*sqrt(err(amu[1])^2+syst^2)
xmax = value(amu[1])+5*sqrt(err(amu[1])^2+syst^2)
x = Float64.(range(xmin, xmax, length=1000))

Fpdf = zeros(1000)
for impr_set in IMPR_SET
    for key in mykeys
        tr_valerr = [value.(res[impr_set][key]),err.(res[impr_set][key])]
        valerr    = [[tr_valerr[1][i],tr_valerr[2][i]] for i=1:length(tr_valerr[1])]
        dist = [Normal(val, err) for (val,err) in valerr]
        Npdf = pdf.(dist, Ref(x))
        if ACC
            Fpdfpre = Fpdf[:]
            for n=1:length(Npdf)
                Fpdf .+= (info[impr_set]["weight"][key][n]/(length(mykeys)*length(IMPR_SET))) .* Npdf[n] ./ (length(IMPR_SET) + length(mykeys))
            end

            fill_between(x,Fpdf,Fpdfpre,alpha=1.0,color=color_dict[impr_set][key[end-1:end]],label="discr. $(key[end-1:end]); set $impr_set")
        else
            Fpdf = zeros(1000)
            for n=1:length(Npdf)
                Fpdf .+= (info[impr_set]["weight"][key][n]/(length(mykeys)*length(IMPR_SET))) .* Npdf[n]
            end

            PyPlot.plot(x,Fpdf,color=color_dict[impr_set][key[end-1:end]],label="discr. $(key[end-1:end]); set $impr_set")
        end
    end
end
if diag != "LO" # we multiply the y axis by a factor 10
    formatter(x, pos) = string(round(10 * x, digits=2))  # Round to 2 decimal places
    ax = gca()
    ax.xaxis.set_major_formatter(PyPlot.matplotlib.ticker.FuncFormatter(formatter))
end
PyPlot.title("Gaussian hit distribution; impr. sets [$(paste_str(IMPR_SET))]")
diag_str = diag == "LO" ? "\\rm{LO}" : (diag == "NLOa&b" ? "\\rm{NLO}_{\\rm{a}\\&\\rm{b}}" : "\\rm{NLO}_{\\rm{$(diag[end])}}")
comp_str = comp[1] == '∆' ? (comp[end] != 'b' ? "\\Delta_{ls}(a_{\\mu})" : "\\Delta_{lc}(b)") : "\\mathrm{$(comp[2]),$(comp[3])}"
Q_str    = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "($Q\\ GeV)" : ""
fact_str = diag == "LO" ? "\\times10^{10}" : "\\times10^{11}"
if wind == "NW"
    xlabel(latexstring("a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]$Q_str$fact_str"))
else
    xlabel(latexstring("\\left(a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]\\right)^{\\rm{$wind}}$Q_str$fact_str"))
end
ylabel("Normalized hits")
legend(loc="best")
tight_layout()
current_ylim = ylim()  # Get current y-limits
# ylim(0, current_ylim[2])  # Set lower bound and keep upper bound
# xlim(110,150)
# ylim(0,0.025)
display(gcf())
close()


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##


@info("Model average plot")

IMPR_SET = readIMPR_SET #  readIMPR_SET  ["1"]  ["2"]  ["1","2"]  ["1old","2"]  ["1","1old","2"]

Norm2DOF = false

mykeys = DictComptoKey[comp]

color_dict = Dict(
    ""     => Dict("cc" => "cyan"),
    "1"    => Dict("ll" => "blue",  "lc" => "red"  ),
    "2"    => Dict("ll" => "green", "lc" => "brown")
); color_dict["1old"] = color_dict["1"]


PVAL = sum(isnothing.(vcat([vcat([getfield.(fitres[impr_set][key],:pval) for key in mykeys]...) for impr_set in IMPR_SET]...))) < 1

weight = vcat([vcat([info[impr_set]["weight"][key] for key in mykeys]...) for impr_set in IMPR_SET]...)./(length(mykeys)*length(IMPR_SET))
res_tot = vcat([vcat([res[impr_set][key] for key in mykeys]...) for impr_set in IMPR_SET]...)

amu, syst = model_average(res_tot, weight); uwerr.(amu)

v =  value(amu[1])
e = err(amu[1]); esyst = sqrt(err(amu[1])^2+syst^2)
p1 = res_tot[:]; uwerr.(p1)
w = weight[:]

chi2chi2exp=[]; chi2dof=[]; pval=[]; mods=[]
for impr_set in IMPR_SET
    for key in mykeys
        push!(chi2chi2exp, getfield.(fitres[impr_set][key],:chi2)./getfield.(fitres[impr_set][key],:chi2exp))
        push!(chi2dof    , getfield.(fitres[impr_set][key],:chi2)./getfield.(fitres[impr_set][key],:dof    ))
        push!(pval       , getfield.(fitres[impr_set][key],:pval))
        push!(mods, modelinfo["label_tot_isov"])
    end
end

chi2chi2exp = vcat(chi2chi2exp...)
chi2dof     = vcat(chi2dof...)
pval        = vcat(pval...)
mods        = vcat(mods...)

fig = figure(figsize=(15,10))

subplots_adjust(hspace=0.1) 
subplot(411)    
ax1 = gca()                
x = 1:length(p1)
y = value.(p1)
dy = err.(p1)

# PyPlot.title("Model Average (diag. $diag; comp. $comp)")
PyPlot.title("Model Average; impr. sets [$(paste_str(IMPR_SET))]")

if comp in ["gCCdisc","gC8disc"]
    fill_between(x[1:Int64(length(x)/2)], v-e, v+e, color="limegreen", alpha=0.2)
    fill_between(x[1:Int64(length(x)/2)], v-esyst, v+esyst, color="limegreen", alpha=0.2)
    IMPR_SET_ = IMPR_SET[1]
else
    fill_between(x, v-e, v+e, color="limegreen", alpha=0.2)
    fill_between(x, v-esyst, v+esyst, color="limegreen", alpha=0.2)
    # color_list = ["blue","purple","red","orange","brown","gray"]
    IMPR_SET_ = IMPR_SET
end
for (i,impr_set) in enumerate(IMPR_SET_)
    for (j,key) in enumerate(mykeys)
        o,f=[(key_len*(i-1)+j-1)*model_len+1,(key_len*(i-1)+j)*model_len]
        errorbar(x[o:f], y[o:f], dy[o:f], fmt="o", mfc="none", color=color_dict[impr_set]["$(key[end-1:end])"], ms=10, capsize=2)
    end
end
setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels
if diag != "LO" # we multiply the y axis by a factor 10
    formatter(x, pos) = string(round(10 * x, digits=2))  # Round to 2 decimal places
    ax = gca()
    ax.yaxis.set_major_formatter(PyPlot.matplotlib.ticker.FuncFormatter(formatter))
end
diag_str = diag == "LO" ? "\\rm{LO}" : (diag == "NLOa&b" ? "\\rm{NLO}_{\\rm{a}\\&\\rm{b}}" : "\\rm{NLO}_{\\rm{$(diag[end])}}")
comp_str = comp[1] == '∆' ? (comp[end] != 'b' ? "\\Delta_{ls}(a_{\\mu})" : "\\Delta_{lc}(b)") : "\\mathrm{$(comp[2]),$(comp[3])}"
Q_str    = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "($Q)" : ""
fact_str = diag == "LO" ? "\\times10^{10}" : "\\times10^{11}"
if wind == "NW"
    ylabel(latexstring("a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]$Q_str$fact_str"))
else
    ylabel(latexstring("\\left(a_{\\mu}^{$comp_str}[\\rm{$(diag_str)}]\\right)^{\\rm{$wind}}$Q_str$fact_str"))
end
ylim([v-6*e,v+6*e])

subplot(412)
ax2=gca()
for (i,impr_set) in enumerate(IMPR_SET_)
    for (j,key) in enumerate(mykeys)
        o,f=[(key_len*(i-1)+j-1)*model_len+1,(key_len*(i-1)+j)*model_len]
        PyPlot.bar(x[o:f], w[o:f], alpha=0.4, color=color_dict[impr_set]["$(key[end-1:end])"], edgecolor="black", linewidth=1.0)
    end
end
# PyPlot.bar(x, w, alpha=0.4, color="royalblue", edgecolor="blue", linewidth=1.5)
setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels
# errorbar(mods, weight_model, 0*dy, color="green")
ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))

subplot(413)
ax3=gca()
for (i,impr_set) in enumerate(IMPR_SET_)
    for (j,key) in enumerate(mykeys)
        o,f=[(key_len*(i-1)+j-1)*model_len+1,(key_len*(i-1)+j)*model_len]
        if !Norm2DOF
            PyPlot.bar(x[o:f], chi2chi2exp[o:f], alpha=0.4, color=color_dict[impr_set]["$(key[end-1:end])"], edgecolor="black", linewidth=1.0)
        else
            PyPlot.bar(x[o:f], chi2dof[o:f], alpha=0.4, color=color_dict[impr_set]["$(key[end-1:end])"], edgecolor="black", linewidth=1.0)
        end
    end
end
setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels
!Norm2DOF ? ylabel(L"$\chi^2/\chi^2_{\mathrm{exp}}$") : ylabel(L"$\chi^2/\mathrm{d.o.f.}$")
ylim([0.0,2.0])

if typeof(pval) != Vector{Nothing}
    subplot(414)
    for (i,impr_set) in enumerate(IMPR_SET_)
        for (j,key) in enumerate(mykeys)
            o,f=[(key_len*(i-1)+j-1)*model_len+1,(key_len*(i-1)+j)*model_len]
            PyPlot.bar(x[o:f], pval[o:f], alpha=0.4, color=color_dict[impr_set]["$(key[end-1:end])"], edgecolor="black", linewidth=1.0)
        end
    end
    ylabel(L"$\rm{p-values}$")
    ylim([0.0,1.0])
end

# PyPlot.xticks(x, mods, rotation=90)
xlabel("Models [n = $model_len]")
tight_layout()
display(fig)
close("all")
