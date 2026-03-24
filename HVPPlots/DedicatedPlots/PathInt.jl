
using PyPlot
using Colors

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18
##-----------------------------------------------

# Oscillatory path integral

phi = collect(range(-10,10,10000))

obs = phi.^3 .+ 5*phi.^2 .- 40*phi .- 50
exp_S_re = 100*cos.(phi.^2)

plot(phi, exp_S_re, color="orange", label=latexstring("\\mathcal{R}\\{\\exp (i\\,S[\\Phi])\\}"))
plot(phi, obs, color="blue", label=latexstring("\\mathcal{O}[\\Phi]"))
xlim(-10,10)
ylim(-200,200)
xlabel(latexstring("\\Phi-\\Phi_{\\rm cl}"))
legend(loc="best")
ax = gca()
axvline(x=0.0, color="black", lw=0.2, alpha=0.7) 
ax.set_xticks([], [])
ax.set_yticks([], [])
display(gcf())
close()


##-----------------------------------------------

# Exponential path integral

phi = collect(range(-10,10,10000))

obs = phi.^3 .+ 5*phi.^2 .- 40*phi .- 50
exp_S = 150*exp.(-0.05*phi.^2)

plot(phi, exp_S, color="orange", label=latexstring("\\exp (-\\,S_E[\\Phi])"))
plot(phi, obs, color="blue", label=latexstring("\\mathcal{O}[\\Phi]"))
xlim(-10,10)
ylim(-200,200)
xlabel(latexstring("\\Phi-\\Phi_{\\rm cl}"))
legend(loc="best")
ax = gca()
axvline(x=0.0, color="black", lw=0.2, alpha=0.7) 
ax.set_xticks([], [])
ax.set_yticks([], [])
display(gcf())
close()
