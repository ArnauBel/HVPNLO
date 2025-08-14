using Revise
using HVPobs
using DelimitedFiles
using ADerrors, PyPlot

using Statistics

##

#1st - Data is read and the correlator computed
path = "/Users/cesc/Desktop/NLO_lattice_implementation/HVP_data/N202/light/raw_data/data_HVP_V1V1"      #path to data
cdata = read_hvp_data(path, "N202");        #read the data
corr = corr_obs(cdata, rw=nothing)      #compute average points
uwerr.(corr.obs)        #to add uncertainties into "corr"

#2nd - The correlator is symmetrized
n = length(corr.obs)
m = div(n, 2)       #calculate midpoint

corr_sym_values = Float64[]     #initialize value and error vectors
corr_sym_errors = Float64[]

# Map "n right -> (n-1) left" and average
push!(corr_sym_values, value(corr.obs[1]))
push!(corr_sym_errors, ADerrors.err(corr.obs[1]))
for i in 2:(m-1)
    sym_value = (value(corr.obs[i]) + value(corr.obs[end - i + 2])) / 2
    sym_error = (ADerrors.err(corr.obs[i]) + ADerrors.err(corr.obs[end - i + 2])) / 2 
    push!(corr_sym_values, sym_value)
    push!(corr_sym_errors, sym_error)
end
push!(corr_sym_values, value(corr.obs[end - m + 2]))
push!(corr_sym_errors, ADerrors.err(corr.obs[end - m + 2]))

#3rd - Plot
errorbar(collect(1:length(corr.obs)), -value.(corr.obs), ADerrors.err.(corr.obs), fmt="s", label="Raw data", capsize=2)     #plot of raw data
errorbar(collect(1:length(corr_sym_values)), -corr_sym_values, corr_sym_errors, fmt="s", label="Correlator (sym)", capsize=2)       #plot of sym data
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title("N202 correlator")
xlabel("t [a]")
ylabel("-G(t)")
ylim((0.5*minimum(filter(x -> x >= 0, -value.(corr.obs))), 2*maximum(-value.(corr.obs))))
yscale("log")
legend(["Raw data","Correlator (sym)"], loc  = "best")
grid("on")
#clf()       #clear figure; comment to stack new plot with the old plots
display(gcf())      #display the figure
close()     #clear figure

## --------------------------------------------------------------------------

#1st - Data is read and the correlator computed
basicPath = "/Users/cesc/Desktop/NLO_lattice_implementation/HVP_data/N202/light/raw_data/"      #path to data
directionPath = ["data_HVP_V1V1","data_HVP_V2V2","data_HVP_V3V3"]

for path_it in 1:3
    path =  basicPath*directionPath[path_it]
    cdata = read_hvp_data(path, "N202");        #read the data
    corr = corr_obs(cdata, rw=nothing)      #compute average points
    uwerr.(corr.obs)        #to add uncertainties into "corr"

    #2nd - The correlator is symmetrized
    n = length(corr.obs)
    m = div(n, 2)       #calculate midpoint

    corr_sym_values = [Float64[],Float64[],Float64[]]     #initialize value and error vectors
    corr_sym_errors = [Float64[],Float64[],Float64[]]

    # Map "n right -> (n-1) left" and average
    push!(corr_sym_values[path_it], value(corr.obs[1]))
    push!(corr_sym_errors[path_it], ADerrors.err(corr.obs[1]))
    for i in 2:(m-1)
        sym_value = (value(corr.obs[i]) + value(corr.obs[end - i + 2])) / 2
        sym_error = (ADerrors.err(corr.obs[i]) + ADerrors.err(corr.obs[end - i + 2])) / 2 
        push!(corr_sym_values[path_it], sym_value)
        push!(corr_sym_errors[path_it], sym_error)
    end
    push!(corr_sym_values[path_it], value(corr.obs[end - m + 2]))
    push!(corr_sym_errors[path_it], ADerrors.err(corr.obs[end - m + 2]))
end

##

function symmetrize_corr(corr_obs::Vector{uwreal})
    n = length(corr_obs)
    m = div(n, 2)

    sym_values = []
    sym_errors = []

    # Map "n right -> (n-1) left" and average
    push!(sym_values, value(corr_obs[1]))
    push!(sym_errors, ADerrors.err(corr_obs[1]))

    for i in 2:(m-1)
        sym_value = (value(corr_obs[i]) + value(corr_obs[end - i + 2])) / 2
        sym_error = (ADerrors.err(corr_obs[i]) + ADerrors.err(corr_obs[end - i + 2])) / 2 
        push!(sym_values, sym_value)
        push!(sym_errors, sym_error)
    end

    push!(sym_values, value(corr_obs[end - m + 2]))
    push!(sym_errors, ADerrors.err(corr_obs[end - m + 2]))

    return sym_values, sym_errors
end


basicPath = "/Users/cesc/Desktop/NLO_lattice_implementation/HVP_data/N202/light/raw_data/data_HVP_"
directionPath = ["V1V1", "V2V2", "V3V3"]

corr_sym_values = Vector{Vector{Float64}}(undef, 3)
corr_sym_errors = Vector{Vector{Float64}}(undef, 3)

for path_it in 1:3
    path = basicPath * directionPath[path_it]
    cdata = read_hvp_data(path, "N202")
    corr = corr_obs(cdata, rw=nothing)
    uwerr.(corr.obs)

    corr_sym_values[path_it], corr_sym_errors[path_it] = symmetrize_corr(corr.obs)
end

#3rd - Plot
for path_it in 1:3
    errorbar(collect(1:length(corr_sym_values[path_it])), -corr_sym_values[path_it], corr_sym_errors[path_it], fmt="s", label="Correlator"*directionPath[path_it], capsize=2)       #plot of sym data
end
axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title("N202 correlator")
xlabel("t [a]")
ylabel("-G(t)")
ylim((0.5*minimum(filter(x -> x >= 0, -value.(corr.obs))), 2*maximum(-value.(corr.obs))))
yscale("log")
legend(directionPath, loc  = "best")
grid("on")
#clf()       #clear figure; comment to stack new plot with the old plots
display(gcf())      #display the figure
close()     #clear figure

##

corr_sym_values_average = mean(corr_sym_values)
corr_sym_errors_average = std(corr_sym_errors)

errorbar(collect(1:length(corr_sym_values_average)), -corr_sym_values_average, corr_sym_errors_average, fmt="s", label="Averaged correlator", capsize=2)       #plot of sym data

axis("tight")
ax = gca()      # get the handle of the current axis (not really used here)
PyPlot.title("N202 correlator")
xlabel("t [a]")
ylabel("-G(t)")
ylim((0.1*minimum(filter(x -> x >= 0, -corr_sym_values_average)), 2*maximum(-value.(corr.obs))))
yscale("log")
legend(["Averaged correlator"], loc  = "best")
grid("on")
#clf()       #clear figure; comment to stack new plot with the old plots
display(gcf())      #display the figure
close()     #clear figure


corr_sym_values_average[end]
corr_sym_errors_average[end]

corr_sym_values_average[end-3]
corr_sym_errors_average[end-3]