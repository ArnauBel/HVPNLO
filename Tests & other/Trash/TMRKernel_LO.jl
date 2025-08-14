######################################################################################
# This file was created by Arnau Beltran
# Here we define the Kernel functions of the LO and NLO diagrams for the HVP
# This Kernels are writen in a polynomial expansion for t << 1 and t >> 1
########################################################################################

##
#  LO example
##

GammaEuler = 0.57721566490153286060651209008240243

an = [1/9 , -169/5400 , -401/88200 , -787/2916000 , -7353/768398400]
bn = [0 , 120/5400 , 210/88200 , 360/2916000 , 3080/768398400]

function Tildef2(t::Vector{Float64}, an::Vector{Float64}, bn::Vector{Float64})
    if length(an) != length(bn)
        error("Length of an and bn are not the same")
    else
        return [sum((an[i] + bn[i] * (log(ti)+GammaEuler)) * ti^(2*i+2) for i in 1:length(an)) for ti in t]
    end
end

##

function inner(t::Float64, an::Vector{Float64}, bn::Vector{Float64})::Float64
    idcs = collect(4:2:12)
    return sum((an .+ bn .* (log(t)+GammaEuler)) .* (t .^ idcs))
end

function Tildef22(t::Vector{Float64}, an::Vector{Float64}, bn::Vector{Float64})
    if length(an) != length(bn)
        error("Length of an and bn are not the same")
    else
        inner2 = t0 -> inner(t0, an, bn)
        return inner2.(t)
    end
end

##
# Function test -> Plot
##

using LaTeXStrings
using PyPlot

t = collect(range(0.0001, stop=0.5, length=1000))

#@time Tildef = Tildef2(t, an, bn)
Tildef = Tildef22(t, an, bn)

plot(t, Tildef, label="Tildef2", color = "green")
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title("TMR Kernel function")
xlabel(L"\hat{t}")
ylabel(L"\tilde{f}(t)")
yscale("linear")
#yticks(range(0, stop=0.004, length=8))
legend([L"\tilde{f}_2(\hat{t})"], loc  = "best")
grid("on")
display(gcf())      #display the figure
close()

##