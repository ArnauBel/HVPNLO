using ADerrors
using HVPobs

using Plots

include("TMRKernel.jl")
export Tildef4aInner, Tildef4a, Tildef4bInner, Tildef4b, Tildef4cInner

const hbarc = 0.1973269804        # GeV * fm

julia_script_directory = @__DIR__
path_coef = joinpath(julia_script_directory, "..", "Coefficients")

##== Tildef4c definition:

function Tildef4cInner(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal}, a2n::Union{Vector{Rational{BigInt}},Vector{Float64}}, b2n::Union{Vector{Rational{BigInt}},Vector{Float64}}, a3n::Union{Vector{Rational{BigInt}},Vector{Float64}}, b3n::Union{Vector{Rational{BigInt}},Vector{Float64}}, aidcs::Vector{Int64}, bidcs::Vector{Int64})
        
    aterm1 = sum(a2n[3:end] .* (t^2 .* (tau .^aidcs) .+ tau^2 .* (t .^aidcs)))
    aterm2 = sum(b2n[3:end] .* (t^2 .* (custom_log(tau)+GammaEuler) .* (tau .^aidcs) .+ tau^2 .* (custom_log(t)+GammaEuler) .* (t .^aidcs)))
    
    bterm1 = 2 * sum(a3n[4:end] .* ((t + tau) .^bidcs .+ (t - tau) .^bidcs .- 2 .* t .^bidcs .- 2 .* tau .^bidcs)) 
    bterm2 = 2 * sum(b3n[4:end] .* ((custom_log(t+tau)+GammaEuler) .* (t + tau) .^bidcs .+ (custom_log(t-tau)+GammaEuler) .* (t - tau) .^bidcs .- 2 .* (custom_log(t)+GammaEuler) .* t .^bidcs .- 2 .* (custom_log(tau)+GammaEuler) .* tau .^bidcs))
    
    #return (32 * pi^4 / massmu^4) * (aterm1 + aterm2 + bterm1 + bterm2)
    return aterm1 + aterm2 + bterm1 + bterm2
end

function Tildef4cSUPa(t::Union{Int64,Float64,uwreal}, cn::Vector{Float64}, cidcs::Vector{Int64})
    t0 = 2.2
    mySUP = exp(-2*t)/sqrt(t) * sum(cn .* (t0/t -  1) .^cidcs)
    #return (32 * pi^4 / massmu^4)*sqrt(pi) * mySUP
    return sqrt(pi) * mySUP
end

function Tildef4cDOM_blue(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal})
    myDOM = -11/6 - t^2/2 + t^2/tau^2 - tau^2/2 + tau^2/t^2 + (t^2*tau^2)/4 - 1/48*pi*(t + tau)*(-45 + 4*(t+tau)^2) + 1/48*pi*(t-tau)*(-45 + 4*(t-tau)^2) + 2*(t^2 + tau^2 - 1)*(GammaEuler + custom_log(t*tau)) + (1 - (t+tau)^2)*custom_log(t+tau) + (1 - (t-tau)^2)*custom_log(t-tau)
    #return (32 * pi^4 / massmu^4) * myDOM
    return myDOM
end
function Tildef4cSUPP_blue(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal}, c2n::Vector{Float64}, c3n::Vector{Float64}, cidcs::Vector{Int64})
    mySUPP = 0.25*(tau^2*Tildef4cSUPa(t,c2n,cidcs) + t^2*Tildef4cSUPa(tau,c2n,cidcs)) + 2*(Tildef4cSUPa(t+tau,c3n,cidcs) + Tildef4cSUPa(t-tau,c3n,cidcs) - 2*Tildef4cSUPa(t,c3n,cidcs) - 2*Tildef4cSUPa(tau,c3n,cidcs))
    #return (32 * pi^4 / massmu^4) * mySUPP
    return mySUPP
end

function Tildef4cDOM_red(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal})
    myDOM = -11/4 - t^2/2 + t^2/tau^2 - tau^2/2 + tau^2/t^2 + (t^2*tau^2)/4 - 1/48*pi*(t + tau)*(-45 + 4*(t+tau)^2) + 2*(t^2 + tau^2 - 1)*(GammaEuler + custom_log(t*tau)) + (1 - (t+tau)^2)*custom_log(t+tau) - GammaEuler
    #return (32 * pi^4 / massmu^4) * myDOM
    return myDOM
end
function Tildef4cSUPP_red(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal}, c2n::Vector{Float64}, c3n::Vector{Float64}, cidcs::Vector{Int64})
    mySUPP = 0.25*(tau^2*Tildef4cSUPa(t,c2n,cidcs) + t^2*Tildef4cSUPa(tau,c2n,cidcs)) + 2*(Tildef4cSUPa(t+tau,c3n,cidcs) - 2*Tildef4cSUPa(t,c3n,cidcs) - 2*Tildef4cSUPa(tau,c3n,cidcs))
    #return (32 * pi^4 / massmu^4) * mySUPP
    return mySUPP
end

function Tildef4cDOM_yellow(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal})
    myDOM =  tau^2/t^2+ 2*(t^2 + tau^2 - 1)custom_log(t) + (1 - (t+tau)^2)*custom_log(t+tau) + (1 - (t-tau)^2)*custom_log(t-tau)    #return (32 * pi^4 / massmu^4) * myDOM
    return myDOM
end
function Tildef4cSUPP_yellow(t::Union{Int64,Float64,uwreal}, tau::Union{Int64,Float64,uwreal}, c2n::Vector{Float64}, c3n::Vector{Float64}, cidcs::Vector{Int64})
    mySUPP = 0.25*tau^2*Tildef4cSUPa(t,c2n,cidcs) + 2*(Tildef4cSUPa(t+tau,c3n,cidcs) + Tildef4cSUPa(t-tau,c3n,cidcs) - 2*Tildef4cSUPa(t,c3n,cidcs))
    #return (32 * pi^4 / massmu^4) * mySUPP
    return mySUPP
end

function Tildef4c(t::Union{Vector{Int64},Vector{Float64},Vector{uwreal}}, path_to_coef::String; sym::Bool=false)
    
    contents4c = read(joinpath(path_to_coef, "NLO_diagram4c.txt"), String)
    a2n, b2n, a3n, b3n, c2n, c3n = [Float64.(vec) for vec in extract_coef.(contents4c, ["a2n","b2n","a3n","b3n","c2n","c3n"])]
    
    length(a2n) != length(b2n) ? error("The lengths of a2n and b2n should be the same") : nothing
    length(a3n) != length(b3n) ? error("The lengths of a3n and b3n should be the same") : nothing
    length(c2n) != length(c3n) ? error("The lengths of c2n and c3n should be the same") : nothing

    aidcs  = collect(4:2:Int64(length(a2n)*2-2))
    bidcs  = collect(6:2:Int64(length(a3n)*2-2))
    cidcs  = collect(0:1:Int64(length(c2n)-1))

    n = length(t)
    result_matrix = Array{Any}(undef, n, n)

    for i in 1:n
        for j in 1:i
            ti = t[i]; tj = t[j]
            if ti + tj <= 3.8
                result_matrix[i, j] = Tildef4cInner(ti,tj,a2n,b2n,a3n,b3n,aidcs,bidcs)
            elseif ti-tj <= 3 && tj >= 1.3
                result_matrix[i, j] = Tildef4cDOM_red(ti, tj) + Tildef4cSUPP_red(ti, tj, c2n, c3n, cidcs) + (GammaEuler-7/4)*(ti-tj)^2 + 2*sum((a3n[3:end].+b3n[3:end].*(GammaEuler+custom_log(ti-tj))).*(ti-tj).^aidcs)
            elseif tj <= 2.4
                result_matrix[i, j] = Tildef4cDOM_yellow(ti, tj) + Tildef4cSUPP_yellow(ti, tj, c2n, c3n, cidcs) + 3*tj^2 + ti^2*sum((a2n[3:end].+b2n[3:end].*(GammaEuler+custom_log(tj))).* tj.^aidcs) - 4*sum((a3n[3:end].+b3n[3:end].*(GammaEuler+custom_log(tj))).* tj.^aidcs)
            else
                result_matrix[i, j] = Tildef4cDOM_blue(ti, tj) + Tildef4cSUPP_blue(ti, tj, c2n, c3n, cidcs)
            end
        end
    end
    if sym
        for i in 1:n
            for j in 1:i-1
                result_matrix[j,i] = result_matrix[i,j]
            end
        end
    end
    #return (32 * pi^4 / massmu^4) .* result_matrix
    return result_matrix
end

##== Analytical estimation and numerical evaluation

using QuadGK
using ProgressBars

t = collect(1:0.005:4)
# t = collect(1.1:0.05:3)

myKernelc = Tildef4c(t,path_coef)./(32 * pi^4 / massmu^4)

Tildef4c_num(u,t,tau) = u/(2*(u-1)^4) * (t^2*(u-1)^2 - 2*u + 2*u*cos((t*(u-1))/sqrt(u))) * (tau^2*(u-1)^2 - 2*u + 2*u*cos((tau*(u-1))/sqrt(u)))

print(quadgk(u -> Tildef4c_num(u, 1, 1), 0, 1))
print(myKernelc[2,2])

n = length(t)
myKernelc_num = Array{Any}(undef, n, n)
myKernelc_err = Array{Any}(undef, n, n)

for i in ProgressBar(1:n)
    for j in 1:i
        myKernelc_num[i,j], myKernelc_err[i,j] = quadgk(u -> Tildef4c_num(u, t[i], t[j]), 0, 1, atol=1e-13)
        myKernelc_num[j,i], myKernelc_err[j,i] = myKernelc_num[i,j], myKernelc_err[i,j]
    end
end

##== Max error 

data = abs.(myKernelc_num-myKernelc)
#myKernelc_err

max_value = maximum(data)
max_index = argmax(data)
max_point = [t[max_index[1]],t[max_index[2]]]

##== Plot

using Plots
gr() # Use the GR backend for plotting

heatmap(t,t,log10.(data), c=:inferno, xlabel=L"$\hat{t}$", ylabel=L"$\hat{\tau}$", title="Error budget in 'log' scale (zoom)", size=(1000, 1000), titlefontsize=20, xguidefontsize=20, yguidefontsize=20, zguidefontsize=20, xtickfont=12, ytickfont=12)

Plots.savefig("heatmap.png")
