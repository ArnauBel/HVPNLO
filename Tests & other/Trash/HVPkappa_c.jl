# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

using Plots
using PyPlot
using Colors

# Data treatment 

IMPR      = true
RENORM    = true
STD_DERIV = false

# Path definition

julia_script_directory = @__DIR__

path_HVP   = joinpath(julia_script_directory, "..", "LMEData", "HVP_data")
path_rw    = joinpath(julia_script_directory, "..", "LMEData", "rwf_deflated")
path_ms    = joinpath(julia_script_directory, "..", "LMEData", "ms_t0_data")
path_fvcPI = joinpath(julia_script_directory, "..", "LMEData", "FSE", "JKMPI")
path_fvcK  = joinpath(julia_script_directory, "..", "LMEData", "FSE", "JKMK")

path_coef = joinpath(julia_script_directory, "..", "Coefficients")

if STD_DERIV
    path_bdio = joinpath(julia_script_directory, "..", "ObsBDIOstd")
else
    path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")
end

# Plot parameters

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] =  true
rcParams["mathtext.fontset"]  = "cm"
rcParams["font.size"] = 13
rcParams["axes.labelsize"] = 22
rcParams["axes.titlesize"] = 18

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302"]
ensInfo = EnsInfo.(ensList)

# kappaC dict

kcd_in = Dict(
    "A653" => Dict("kappaC" => 0.119743, "kappaC_err" => 0.000017, "kappaC_sim" => 0.119743, "kappaC_sim_plus" => 0.119793),
    "A654" => Dict("kappaC" => 0.120079, "kappaC_err" => 0.000025, "kappaC_sim" => 0.120177, "kappaC_sim_plus" => 0.120227),
    "H101" => Dict("kappaC" => 0.122897, "kappaC_err" => 0.000018, "kappaC_sim" => 0.122908, "kappaC_sim_plus" => 0.122938),
    "H102" => Dict("kappaC" => 0.123041, "kappaC_err" => 0.000026, "kappaC_sim" => 0.123050, "kappaC_sim_plus" => 0.123080),
    "U101" => Dict("kappaC" => 0.123244, "kappaC_err" => 0.000019, "kappaC_sim" => 0.123251, "kappaC_sim_plus" => 0.123281),
    "H105" => Dict("kappaC" => 0.123244, "kappaC_err" => 0.000019, "kappaC_sim" => 0.123251, "kappaC_sim_plus" => 0.123281),
    "N101" => Dict("kappaC" => 0.123244, "kappaC_err" => 0.000019, "kappaC_sim" => 0.123251, "kappaC_sim_plus" => 0.123281),
    "C101" => Dict("kappaC" => 0.123362, "kappaC_err" => 0.000015, "kappaC_sim" => 0.123367, "kappaC_sim_plus" => 0.123397),
    "B450" => Dict("kappaC" => 0.125093, "kappaC_err" => 0.000017, "kappaC_sim" => 0.125089, "kappaC_sim_plus" => 0.125129),
    "S400" => Dict("kappaC" => 0.125252, "kappaC_err" => 0.000020, "kappaC_sim" => 0.125267, "kappaC_sim_plus" => 0.125317),
    "N401" => Dict("kappaC" => 0.125439, "kappaC_err" => 0.000015, "kappaC_sim" => 0.125447, "kappaC_sim_plus" => 0.125477),
    "N451" => Dict("kappaC" => 0.125439, "kappaC_err" => 0.000015, "kappaC_sim" => 0.125447, "kappaC_sim_plus" => 0.125477),
    "D450" => Dict("kappaC" => 0.125585, "kappaC_err" => 0.000007, "kappaC_sim" => 0.125585, "kappaC_sim_plus" => 0.125635),
    "D452" => Dict("kappaC" => 0.125645, "kappaC_err" => 0.000005, "kappaC_sim" => 0.125640, "kappaC_sim_plus" => 0.125690),
    "H200" => Dict("kappaC" => 0.127579, "kappaC_err" => 0.000016, "kappaC_sim" => 0.127626, "kappaC_sim_plus" => 0.127666),
    "N202" => Dict("kappaC" => 0.127579, "kappaC_err" => 0.000016, "kappaC_sim" => 0.127626, "kappaC_sim_plus" => 0.127666),
    "N203" => Dict("kappaC" => 0.127714, "kappaC_err" => 0.000011, "kappaC_sim" => 0.127713, "kappaC_sim_plus" => 0.127733),
    "N200" => Dict("kappaC" => 0.127858, "kappaC_err" => 0.000007, "kappaC_sim" => 0.127859, "kappaC_sim_plus" => 0.127879),
    "D200" => Dict("kappaC" => 0.127986, "kappaC_err" => 0.000006, "kappaC_sim" => 0.127986, "kappaC_sim_plus" => 0.127956),
    "E250" => Dict("kappaC" => 0.128052, "kappaC_err" => 0.000005, "kappaC_sim" => 0.128054, "kappaC_sim_plus" => 0.128064),
    "N300" => Dict("kappaC" => 0.130099, "kappaC_err" => 0.000018, "kappaC_sim" => 0.130099, "kappaC_sim_plus" => 0.130149),
    "N302" => Dict("kappaC" => 0.130247, "kappaC_err" => 0.000009, "kappaC_sim" => 0.130243, "kappaC_sim_plus" => 0.130263),
    "J303" => Dict("kappaC" => 0.130362, "kappaC_err" => 0.000009, "kappaC_sim" => 0.130362, "kappaC_sim_plus" => 0.130382),
    "E300" => Dict("kappaC" => 0.130432, "kappaC_err" => 0.000010, "kappaC_sim" => 0.130421, "kappaC_sim_plus" => 0.130400),
    "J500" => Dict("kappaC" => 0.131663, "kappaC_err" => 0.000016, "kappaC_sim" => 0.131644, "kappaC_sim_plus" => 0.131600)
)

@info("Ready")

##
par = []
for ens in ensInfo

    println("- Reading HVP")

    pdata = joinpath(path_bdio,"HVP&FVC","NW",ens.id,"HVP")

    fb = BDIO_open(joinpath(pdata,"$(ens.id)_HVPNLOa&b_set1"),"r")
    val = Dict{String, Dict{String, uwreal}}()
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        nobs = d["nobs"]
        dims = d["dimensions"]
        ks = collect(d["keys"])
        val["HVP"] =  ALPHAdobs_read_next(fb, keys=ks)
    end
    BDIO_close!(fb)
    info = load(joinpath(pdata,"$(ens.id)_HVPNLOa&b_BMinfo_set1.jld2"), "HVPinfo")
    HVP = merge(val,info) 

    obs = [HVP["HVP"]["gcc_lc_conn"],HVP["HVP"]["gcc_lc_conn_p"]]
    kappaC = [kcd_in[ens.id]["kappaC_sim"],kcd_in[ens.id]["kappaC_sim_plus"]]
    kappaC_tar = uwreal([kcd_in[ens.id]["kappaC"],kcd_in[ens.id]["kappaC_err"]], "kappaC target")

    @. lin_model(x,p) = p[1] + p[2] * x

    fit  = fit_routine(lin_model, kappaC, obs, 2)
    par = fit.param
    obs_tar = lin_model(kappaC_tar, par)[1]

    uwerr(kappaC_tar)
    uwerr(obs_tar)
    obs_val = value.(obs)
    obs_err = err.(obs)

    minimum(vcat(value(kappaC_tar),kappaC...))

    figure()
    PyPlot.title("$(ens.id)")
    errorbar(kappaC, obs_val, obs_err, fmt="d", capsize=2, mfc="none", color="black")
    errorbar(value(kappaC_tar), value(obs_tar), xerr = err(kappaC_tar), yerr=err(obs_tar), capsize=2, color="red", fmt="d")
    kappaC_arr = vcat(value(kappaC_tar)+err(kappaC_tar),value(kappaC_tar)-err(kappaC_tar),kappaC...)
    xarr = Float64.(range(minimum(kappaC_arr), maximum(kappaC_arr), length=100))
    yarr = lin_model(xarr, par); uwerr.(yarr)
    fill_between(xarr, value.(yarr) .- err.(yarr), value.(yarr) .+ err.(yarr), color="royalblue", alpha=0.2 )
    ylabel(L"$a_\mu^{\rm{hvp}}[\rm{NLO}^{\rm{cc-conn}}]$")
    xlabel(L"$\kappa_c$")
    tight_layout()
    display(gcf())
    PyPlot.savefig(joinpath(julia_script_directory,"..","Slides & Plots","Plots kappaC corr",ens.id))
    close()
end

##

