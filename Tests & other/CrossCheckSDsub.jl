## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> 2206.06582v2: Intermediate window paper (Tables VIII and IX) <==========================##


using ALPHAio, ADerrors
import ADerrors: err

using BDIO

using Printf

julia_script_directory = @__DIR__

ensList = ["A654","B450","C101","C102","D150","D200","D201","D450","D451","D452","E250","E300","F300","H101","H102","H200","J303","J304","J306","J307","J500","J501","N101","N200","N202","N203","N300","N302","N451","N452","S400"]
ensInfo = EnsInfo.(ensList)

# Substracted iso-vector contributions (Table VI)

IsoVecRes = Dict(
    "A653" => Dict("1old" => Dict("ll" => uwreal([39.138, 0.075], "result"), "lc" => uwreal([40.010, 0.045], "result")), "2" => Dict("ll" => uwreal([28.547, 0.047], "result"), "lc" => uwreal([34.913, 0.031], "result"))),
    "A654" => Dict("1old" => Dict("ll" => uwreal([39.171, 0.082], "result"), "lc" => uwreal([40.102, 0.056], "result")), "2" => Dict("ll" => uwreal([29.198, 0.063], "result"), "lc" => uwreal([35.484, 0.044], "result"))),
    "H101" => Dict("1old" => Dict("ll" => uwreal([38.650, 0.041], "result"), "lc" => uwreal([39.252, 0.030], "result")), "2" => Dict("ll" => uwreal([31.712, 0.047], "result"), "lc" => uwreal([35.985, 0.036], "result"))),
    "H102" => Dict("1old" => Dict("ll" => uwreal([38.668, 0.040], "result"), "lc" => uwreal([39.293, 0.029], "result")), "2" => Dict("ll" => uwreal([31.712, 0.047], "result"), "lc" => uwreal([35.985, 0.036], "result"))),
    "H105" => Dict("1old" => Dict("ll" => uwreal([38.695, 0.063], "result"), "lc" => uwreal([39.357, 0.058], "result")), "2" => Dict("ll" => uwreal([32.018, 0.058], "result"), "lc" => uwreal([36.279, 0.056], "result"))),
    "N101" => Dict("1old" => Dict("ll" => uwreal([38.823, 0.034], "result"), "lc" => uwreal([39.495, 0.022], "result")), "2" => Dict("ll" => uwreal([32.144, 0.017], "result"), "lc" => uwreal([36.412, 0.022], "result"))),
    "C101" => Dict("1old" => Dict("ll" => uwreal([38.782, 0.038], "result"), "lc" => uwreal([39.461, 0.026], "result")), "2" => Dict("ll" => uwreal([32.366, 0.024], "result"), "lc" => uwreal([36.488, 0.027], "result"))),
    "C102" => Dict("1old" => Dict("ll" => uwreal([38.834, 0.039], "result"), "lc" => uwreal([39.526, 0.024], "result")), "2" => Dict("ll" => uwreal([32.411, 0.028], "result"), "lc" => uwreal([36.626, 0.026], "result"))),
    "D150" => Dict("1old" => Dict("ll" => uwreal([38.738, 0.035], "result"), "lc" => uwreal([39.427, 0.022], "result")), "2" => Dict("ll" => uwreal([32.660, 0.036], "result"), "lc" => uwreal([36.832, 0.014], "result"))),
    "B450" => Dict("1old" => Dict("ll" => uwreal([37.916, 0.031], "result"), "lc" => uwreal([38.351, 0.030], "result")), "2" => Dict("ll" => uwreal([32.874, 0.031], "result"), "lc" => uwreal([35.804, 0.029], "result"))),
    "S400" => Dict("1old" => Dict("ll" => uwreal([39.093, 0.038], "result"), "lc" => uwreal([38.445, 0.037], "result")), "2" => Dict("ll" => uwreal([32.144, 0.022], "result"), "lc" => uwreal([36.412, 0.022], "result"))),
    "N451" => Dict("1old" => Dict("ll" => uwreal([38.202, 0.015], "result"), "lc" => uwreal([38.678, 0.013], "result")), "2" => Dict("ll" => uwreal([33.582, 0.016], "result"), "lc" => uwreal([36.485, 0.013], "result"))),
    "D450" => Dict("1old" => Dict("ll" => uwreal([38.232, 0.011], "result"), "lc" => uwreal([38.722, 0.008], "result")), "2" => Dict("ll" => uwreal([33.817, 0.012], "result"), "lc" => uwreal([36.708, 0.012], "result"))),
    "D451" => Dict("1old" => Dict("ll" => uwreal([38.265, 0.012], "result"), "lc" => uwreal([38.762, 0.008], "result")), "2" => Dict("ll" => uwreal([33.806, 0.013], "result"), "lc" => uwreal([36.702, 0.014], "result"))),
    "D452" => Dict("1old" => Dict("ll" => uwreal([38.218, 0.015], "result"), "lc" => uwreal([38.713, 0.013], "result")), "2" => Dict("ll" => uwreal([33.556, 0.015], "result"), "lc" => uwreal([36.853, 0.013], "result"))),
    "H200" => Dict("1old" => Dict("ll" => uwreal([36.992, 0.082], "result"), "lc" => uwreal([37.229, 0.082], "result")), "2" => Dict("ll" => uwreal([34.122, 0.081], "result"), "lc" => uwreal([35.725, 0.080], "result"))),
    "N202" => Dict("1old" => Dict("ll" => uwreal([37.169, 0.032], "result"), "lc" => uwreal([37.409, 0.032], "result")), "2" => Dict("ll" => uwreal([34.274, 0.034], "result"), "lc" => uwreal([35.878, 0.033], "result"))),
    "J500" => Dict("1old" => Dict("ll" => uwreal([35.529, 0.032], "result"), "lc" => uwreal([35.517, 0.030], "result")), "2" => Dict("ll" => uwreal([35.008, 0.029], "result"), "lc" => uwreal([35.200, 0.029], "result"))),
    "J501" => Dict("1old" => Dict("ll" => uwreal([35.704, 0.042], "result"), "lc" => uwreal([35.696, 0.045], "result")), "2" => Dict("ll" => uwreal([35.227, 0.040], "result"), "lc" => uwreal([35.418, 0.044], "result")))
)



∆lsRes = Dict(
    "A654" => Dict("1old" => Dict("ll" => uwreal([0.048, 0.014], "result"), "lc" => uwreal([0.105, 0.014], "result")), "2" => Dict("ll" => uwreal([0.435, 0.015], "result"), "lc" => uwreal([0.401, 0.014], "result"))),
    "H102" => Dict("1old" => Dict("ll" => uwreal([0.044, 0.010], "result"), "lc" => uwreal([0.074, 0.010], "result")), "2" => Dict("ll" => uwreal([0.257, 0.010], "result"), "lc" => uwreal([0.244, 0.010], "result"))),
    "H105" => Dict("1old" => Dict("ll" => uwreal([0.126, 0.020], "result"), "lc" => uwreal([0.187, 0.020], "result")), "2" => Dict("ll" => uwreal([0.546, 0.019], "result"), "lc" => uwreal([0.520, 0.019], "result"))),
    "N101" => Dict("1old" => Dict("ll" => uwreal([0.147, 0.005], "result"), "lc" => uwreal([0.208, 0.005], "result")), "2" => Dict("ll" => uwreal([0.562, 0.005], "result"), "lc" => uwreal([0.535, 0.005], "result"))),
    "C101" => Dict("1old" => Dict("ll" => uwreal([0.182, 0.012], "result"), "lc" => uwreal([0.264, 0.011], "result")), "2" => Dict("ll" => uwreal([0.720, 0.011], "result"), "lc" => uwreal([0.691, 0.010], "result"))),
    "C102" => Dict("1old" => Dict("ll" => uwreal([0.223, 0.008], "result"), "lc" => uwreal([0.318, 0.006], "result")), "2" => Dict("ll" => uwreal([0.823, 0.006], "result"), "lc" => uwreal([0.790, 0.006], "result"))),
    "D150" => Dict("1old" => Dict("ll" => uwreal([0.215, 0.010], "result"), "lc" => uwreal([0.315, 0.007], "result")), "2" => Dict("ll" => uwreal([0.883, 0.007], "result"), "lc" => uwreal([0.848, 0.006], "result"))),
    "S400" => Dict("1old" => Dict("ll" => uwreal([0.073, 0.012], "result"), "lc" => uwreal([0.097, 0.012], "result")), "2" => Dict("ll" => uwreal([0.239, 0.012], "result"), "lc" => uwreal([0.234, 0.012], "result"))),
    "N451" => Dict("1old" => Dict("ll" => uwreal([0.172, 0.004], "result"), "lc" => uwreal([0.217, 0.004], "result")), "2" => Dict("ll" => uwreal([0.468, 0.004], "result"), "lc" => uwreal([0.458, 0.004], "result"))),
    "D450" => Dict("1old" => Dict("ll" => uwreal([0.248, 0.004], "result"), "lc" => uwreal([0.309, 0.003], "result")), "2" => Dict("ll" => uwreal([0.656, 0.003], "result"), "lc" => uwreal([0.643, 0.003], "result"))),
    "D451" => Dict("1old" => Dict("ll" => uwreal([0.299, 0.005], "result"), "lc" => uwreal([0.368, 0.003], "result")), "2" => Dict("ll" => uwreal([0.754, 0.003], "result"), "lc" => uwreal([0.738, 0.003], "result"))),
    "D452" => Dict("1old" => Dict("ll" => uwreal([0.280, 0.006], "result"), "lc" => uwreal([0.350, 0.005], "result")), "2" => Dict("ll" => uwreal([0.763, 0.005], "result"), "lc" => uwreal([0.746, 0.005], "result"))),
    "N203" => Dict("1old" => Dict("ll" => uwreal([0.092, 0.012], "result"), "lc" => uwreal([0.110, 0.012], "result")), "2" => Dict("ll" => uwreal([0.202, 0.012], "result"), "lc" => uwreal([0.204, 0.011], "result"))),
    "N200" => Dict("1old" => Dict("ll" => uwreal([0.212, 0.013], "result"), "lc" => uwreal([0.240, 0.012], "result")), "2" => Dict("ll" => uwreal([0.409, 0.012], "result"), "lc" => uwreal([0.407, 0.012], "result"))),
    "D200" => Dict("1old" => Dict("ll" => uwreal([0.303, 0.011], "result"), "lc" => uwreal([0.343, 0.011], "result")), "2" => Dict("ll" => uwreal([0.581, 0.011], "result"), "lc" => uwreal([0.578, 0.010], "result"))),
    "D201" => Dict("1old" => Dict("ll" => uwreal([0.356, 0.006], "result"), "lc" => uwreal([0.400, 0.005], "result")), "2" => Dict("ll" => uwreal([0.660, 0.005], "result"), "lc" => uwreal([0.656, 0.005], "result"))),
    "E250" => Dict("1old" => Dict("ll" => uwreal([0.355, 0.013], "result"), "lc" => uwreal([0.402, 0.012], "result")), "2" => Dict("ll" => uwreal([0.688, 0.012], "result"), "lc" => uwreal([0.685, 0.011], "result"))),
    "N302" => Dict("1old" => Dict("ll" => uwreal([0.141, 0.015], "result"), "lc" => uwreal([0.148, 0.014], "result")), "2" => Dict("ll" => uwreal([0.201, 0.014], "result"), "lc" => uwreal([0.200, 0.014], "result"))),
    "J303" => Dict("1old" => Dict("ll" => uwreal([0.301, 0.012], "result"), "lc" => uwreal([0.318, 0.011], "result")), "2" => Dict("ll" => uwreal([0.416, 0.012], "result"), "lc" => uwreal([0.418, 0.011], "result"))),
    "J304" => Dict("1old" => Dict("ll" => uwreal([0.415, 0.008], "result"), "lc" => uwreal([0.435, 0.008], "result")), "2" => Dict("ll" => uwreal([0.558, 0.008], "result"), "lc" => uwreal([0.558, 0.007], "result"))),
    "E300" => Dict("1old" => Dict("ll" => uwreal([0.419, 0.005], "result"), "lc" => uwreal([0.442, 0.004], "result")), "2" => Dict("ll" => uwreal([0.573, 0.005], "result"), "lc" => uwreal([0.576, 0.004], "result"))),
    "J501" => Dict("1old" => Dict("ll" => uwreal([0.163, 0.012], "result"), "lc" => uwreal([0.163, 0.013], "result")), "2" => Dict("ll" => uwreal([0.196, 0.012], "result"), "lc" => uwreal([0.192, 0.013], "result")))
)



Qlist = [3.5, 4.0, 5.0, 6.0, 7.0, 8.0]


##------ 

diag = "LO"
comp = "∆ls_amu"

Q = 5.0

STD_DERIV = false
FVCbool   = true
t0SHIFT   = true

IMPR_SET = ["1old","2"]

if comp == "33"
    paperRes = IsoVecRes
    mult = 1
elseif comp == "∆ls_amu"
    paperRes = ∆lsRes
    mult = -1/3
end

qarg = findfirst(x -> x == Q, Qlist)

ensBool = [ensid in keys(paperRes) for ensid in getfield.(ensInfo,:id)]

if STD_DERIV
    path_bdio = joinpath(julia_script_directory, "..", "..", "ObsBDIOstd")
else
    path_bdio = joinpath(julia_script_directory, "..", "..", "ObsBDIO")
end

Res = Dict()
for ens in ensInfo[ensBool]
    println("- Reading data ensemble: $(ens.id)")

    println("   - Reading t0...")

    t0 = BDIOread_t0(path_bdio, ens)

    println("   - Reading HVP...    [applying systematics]")

    HVP = Dict()
    for impr_set in IMPR_SET
        hvp, info = BDIOread_HVPens(path_bdio,diag,"SDsub",ens,impr_set,info=true)
    
        HVP[impr_set] = apply_syst_HVP(hvp,info["HVPsyst"],diag,"SDsub",ens.id)
    end

    if FVCbool
        println("   - Reading FVC...    [applying systematics]")

        fvc = BDIOread_FVCens(path_bdio,diag,"SDsub",ens)
        
        FVC = apply_syst_FVC(fvc,diag,"SDsub",ens.id,IMPR_SET=IMPR_SET)

        println("   - aµ = HVP + FVC")

        Res_ = HVP_VolCorrect(HVP,FVC,diag,IMPR_SET=IMPR_SET)
    else
        Res_ = HVP
    end
    
    if t0SHIFT
        println("   - Performing t0 shift...")

        for impr_set in IMPR_SET
            KEY = comp == "33" ? ["g$(comp)_ll","g$(comp)_lc"] : ["$(comp)_ll","$(comp)_lc"]
            for key in KEY
                uwerr.(Res_[impr_set][key]); der = [mchist(Res_[impr_set][key][i], "sqrtt0 [fm]")[1] for i=1:length(Res_[impr_set][key])] ./ artificial_err
                Res_[impr_set][key] = Res_[impr_set][key] .+ value(0.1443 - sqrtt0_ph_Regensburg) .* der
            end
        end
    end

    Res[ens.id] = Res_
end

println("\n\n\n")

##

# Initialize the LaTeX table string
latex_table = """
\\begin{table}[h]
\\centering
\\begin{tabular}{|c|c|c|c|c|}
\\hline
ens & Set 1 ll & Set 1 lc & Set 2 ll & Set 2 lc \\\\
\\hline
\\hline

"""

# Iterate over each ens.id
for ens in ensInfo[ensBool]
    row = ens.id * " & "
    for impr_set in IMPR_SET
        for discr in ["ll","lc"]
            if comp == "33"
                res_value = mult*Res[ens.id][impr_set]["g$(comp)_$(discr)"][qarg]; uwerr(res_value)
            elseif comp == "∆ls_amu"
                res_value = mult*Res[ens.id][impr_set]["$(comp)_$(discr)"][1]; uwerr(res_value)
            end
            if FVCbool
                pap_value = paperRes[ens.id][impr_set][discr]; uwerr(pap_value)
            else
                # if comp == "33"
                #     if ens.kappa_l == ens.kappa_s
                #         pap_value = paperRes[ens.id][impr_set][discr] - (FVCRes[ens.id]["HP(t<t*)"]+FVCRes[ens.id]["HP(t>t*)"]); uwerr(pap_value)
                #     else
                #         pap_value = paperRes[ens.id][impr_set][discr] - (FVCRes[ens.id]["HP(t<t*)"]+FVCRes[ens.id]["HP(t>t*)"]+FVCRes[ens.id]["Kaon loop"]); uwerr(pap_value)
                #     end
                # elseif comp == "88"
                #     if ens.kappa_l == ens.kappa_s
                #         pap_value = paperRes[ens.id][impr_set][discr] - (1/3) * (FVCRes[ens.id]["HP(t<t*)"]+FVCRes[ens.id]["HP(t>t*)"]); uwerr(pap_value)
                #     else
                #         pap_value = paperRes[ens.id][impr_set][discr] - 3/2 * (FVCRes[ens.id]["Kaon loop"]); uwerr(pap_value) # 3/2*
                #     end
                # end
            end
            difference = res_value - pap_value; uwerr(difference)
            sigma = abs(value(difference))/err(difference)
            cell = "\\begin{tabular}{@{}c@{}} $(round(value(pap_value),digits=3)) \$\\pm\$ $(round(err(pap_value),digits=3)) \\\\ $(round(value(res_value),digits=3)) \$\\pm\$ $(round(err(res_value),digits=3)) \\\\ $(round(sigma,digits=2)) \\end{tabular}"
            row *= cell * " & "
        end
    end
    # Remove the last " & " and add a line break
    row = row[1:end-3] * " \\\\\n \\hline \n"
    latex_table *= row
end

# Close the LaTeX table
if comp == "33"
    latex_table *= """
    \\hline
    \\end{tabular}
    \\caption{Comparison iso-vector}
    \\end{table}
    """
elseif comp == "∆ls_amu"
    latex_table *= """
    \\hline
    \\end{tabular}
    \\caption{Comparison iso-scalar-iso-vector}
    \\end{table}
    """
end

# Print or save the LaTeX table
println(latex_table)



## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> 1904.03120v1: 2019 (g-2)µ paper (Tables VIII and IX) <==========================##


using ALPHAio, ADerrors
import ADerrors: err

using BDIO


julia_script_directory = @__DIR__
path_bdio = joinpath(julia_script_directory, "..", "..", "ObsBDIO")

ensList = ["H101", "H102", "N101", "C101", "B450", "S400", "N202", "N203", "N200", "D200", "E250", "N300", "N302", "J303"]
ensInfo = EnsInfo.(ensList)

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

data_dict = Dict(
    "H101" => Dict(
        "interpolation" => 0.8615,
        "simulated" => [
            uwreal([0.8513, 0.0006], "sim H101"),
            uwreal([0.8614, 0.0006], "sim H101"),
            uwreal([0.8714, 0.0006], "sim H101"),
            uwreal([0.8813, 0.0006], "sim H101")
        ]
    ),
    "H102" => Dict(
        "interpolation" => 0.8615,
        "simulated" => [
            uwreal([0.8528, 0.0009], "sim H102"),
            uwreal([0.8629, 0.0009], "sim H102"),
            uwreal([0.8729, 0.0009], "sim H102"),
            uwreal([0.8828, 0.0009], "sim H102")
        ]
    ),
    "N101" => Dict(
        "interpolation" => 0.8615,
        "simulated" => [
            uwreal([0.8495, 0.0006], "sim N101"),
            uwreal([0.8563, 0.0006], "sim N101"),
            uwreal([0.8630, 0.0006], "sim N101"),
            uwreal([0.8697, 0.0006], "sim N101")
        ]
    ),
    "C101" => Dict(
        "interpolation" => 0.8615,
        "simulated" => [
            uwreal([0.8467, 0.0004], "sim C101"),
            uwreal([0.8534, 0.0004], "sim C101"),
            uwreal([0.8602, 0.0004], "sim C101"),
            uwreal([0.8669, 0.0004], "sim C101")
        ]
    ),
    "B450" => Dict(
        "interpolation" => 0.7615,
        "simulated" => [
            uwreal([0.7543, 0.0008], "sim B450"),
            uwreal([0.7614, 0.0008], "sim B450"),
            uwreal([0.7683, 0.0008], "sim B450"),
            uwreal([0.7752, 0.0008], "sim B450")
        ]
    ),
    "S400" => Dict(
        "interpolation" => 0.7615,
        "simulated" => [
            uwreal([0.7457, 0.0007], "sim S400"),
            uwreal([0.7528, 0.0007], "sim S400"),
            uwreal([0.7599, 0.0007], "sim S400"),
            uwreal([0.7669, 0.0007], "sim S400")
        ]
    ),
    "N202" => Dict(
        "interpolation" => 0.6410,
        "simulated" => [
            uwreal([0.6347, 0.0006], "sim N202"),
            uwreal([0.6384, 0.0006], "sim N202"),
            uwreal([0.6421, 0.0006], "sim N202"),
            uwreal([0.6458, 0.0006], "sim N202")
        ]
    ),
    "N203" => Dict(
        "interpolation" => 0.6410,
        "simulated" => [
            uwreal([0.6341, 0.0004], "sim N203"),
            uwreal([0.6416, 0.0004], "sim N203"),
            uwreal([0.6490, 0.0004], "sim N203"),
            uwreal([0.6563, 0.0004], "sim N203")
        ]
    ),
    "N200" => Dict(
        "interpolation" => 0.6410,
        "simulated" => [
            uwreal([0.6320, 0.0003], "sim N200"),
            uwreal([0.6395, 0.0003], "sim N200"),
            uwreal([0.6469, 0.0003], "sim N200"),
            uwreal([0.6543, 0.0003], "sim N200")
        ]
    ),
    "D200" => Dict(
        "interpolation" => 0.6410,
        "simulated" => [
            uwreal([0.6330, 0.0002], "sim D200"),
            uwreal([0.6367, 0.0002], "sim D200"),
            uwreal([0.6405, 0.0002], "sim D200"),
            uwreal([0.6442, 0.0002], "sim D200")
        ]
    ),
    "E250" => Dict(
        "interpolation" => 0.6410,
        "simulated" => [
            uwreal([0.6317, 0.0002], "sim E250"),
            uwreal([0.6392, 0.0002], "sim E250"),
            uwreal([0.6467, 0.0002], "sim E250"),
            uwreal([0.6541, 0.0002], "sim E250")
        ]
    ),
    "N300" => Dict(
        "interpolation" => 0.4969,
        "simulated" => [
            uwreal([0.5126, 0.0007], "sim N300"),
            uwreal([0.5008, 0.0007], "sim N300"),
            uwreal([0.4888, 0.0007], "sim N300"),
            uwreal([0.4766, 0.0007], "sim N300")
        ]
    ),
    "N302" => Dict(
        "interpolation" => 0.4969,
        "simulated" => [
            uwreal([0.5066, 0.0003], "sim N302"),
            uwreal([0.4947, 0.0003], "sim N302"),
            uwreal([0.4827, 0.0003], "sim N302"),
            uwreal([0.4704, 0.0003], "sim N302")
        ]
    ),
    "J303" => Dict(
        "interpolation" => 0.4969,
        "simulated" => [
            uwreal([0.5033, 0.0004], "sim J303"),
            uwreal([0.4954, 0.0004], "sim J303"),
            uwreal([0.4873, 0.0004], "sim J303"),
            uwreal([0.4792, 0.0003], "sim J303")
        ]
    )
)


fb = BDIO_open(joinpath(path_bdio,"kappaC_tar","mDs_prime"),"r")
amDs_SU3 = Dict{String,uwreal}()
i = 0
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    i=i+1
    if i==1
        mDs_ph = ALPHAdobs_read_next(fb)
    else
        ks = collect(d["keys"])
        amDs_SU3 = ALPHAdobs_read_next(fb, keys=ks)
    end
end


KappaC_tar = Dict{String,uwreal}(); myKappaC_tar = Dict{String,uwreal}()
mDs_sim = Dict{String,Vector{uwreal}}(); mymDs_sim = Dict{String,Vector{uwreal}}()
for ens in ensInfo
    fb = BDIO_open(joinpath(path_bdio,"kappaC_tar",ens.id,"$(ens.id)_kappaC"),"r")
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        nobs = d["nobs"]
        dims = d["dimensions"]
        myKappaC_tar[ens.id] = ALPHAdobs_read_next(fb)
    end
    uwerr(myKappaC_tar[ens.id])

    fb = BDIO_open(joinpath(path_bdio,"kappaC_tar",ens.id,"$(ens.id)_mDs"),"r")
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        nobs = d["nobs"]
        dims = d["dimensions"]
        sz = tuple(d["size"]...)
        mymDs_sim[ens.id] = ALPHAdobs_read_next(fb, size=sz)
    end
    uwerr.(mymDs_sim[ens.id])

    KappaC_tar[ens.id] = uwreal([kcd_in[ens.id]["kappaC"],kcd_in[ens.id]["kappaC_err"]],"kappaC_int"); uwerr(KappaC_tar[ens.id])
    mDs_sim[ens.id] = data_dict[ens.id]["simulated"]; uwerr.(mDs_sim[ens.id])
end


##


# Initialize the LaTeX table string
latex_table = """
\\begin{table}[h]
\\centering
\\begin{tabular}{|c|c|c c c c|c|}
\\hline
ens & a mDs (int) & a mDs (sim) &  &  &  & kappaC (sim) \\\\
\\hline
\\hline

"""

# Iterate over each ens.id
for ens in ensInfo
    row = ens.id * " & "

    cell = "\\begin{tabular}{@{}c@{}} $(data_dict[ens.id]["interpolation"]) \\\\ $(round(value(amDs_SU3["$(ens.beta)"]),digits=4)) \\\\ ... \\end{tabular}"
    row *= cell * " & "

    difference = mDs_sim[ens.id] .- mymDs_sim[ens.id]; uwerr.(difference)
    sigma = abs.(value.(difference))./err.(difference)
    # cell = "\\begin{tabular}{@{}c@{}} $(round(value(mDs_sim[ens.id][1]),digits=4)) +/- $(round(err(mDs_sim[ens.id][1]),digits=4)) & $(round(value(mDs_sim[ens.id][2]),digits=4)) +/- $(round(err(mDs_sim[ens.id][2]),digits=4)) & $(round(value(mDs_sim[ens.id][3]),digits=4)) +/- $(round(err(mDs_sim[ens.id][3]),digits=4)) & $(round(value(mDs_sim[ens.id][4]),digits=4)) +/- $(round(err(mDs_sim[ens.id][4]),digits=4)) \\\\ $(round(value(mymDs_sim[ens.id][1]),digits=4)) +/- $(round(err(mymDs_sim[ens.id][1]),digits=4)) & $(round(value(mymDs_sim[ens.id][2]),digits=4)) +/- $(round(err(mymDs_sim[ens.id][2]),digits=4)) & $(round(value(mymDs_sim[ens.id][3]),digits=4)) +/- $(round(err(mymDs_sim[ens.id][3]),digits=4)) & $(round(value(mymDs_sim[ens.id][4]),digits=4)) +/- $(round(err(mymDs_sim[ens.id][4]),digits=4)) \\\\ $(round(sigma[1],digits=2)) & $(round(sigma[2],digits=2)) & $(round(sigma[3],digits=2)) & $(round(sigma[4],digits=2)) \\end{tabular}"
    
    for i=collect(1:4)
        cell = "\\begin{tabular}{@{}c@{}} $(round(value(mDs_sim[ens.id][i]),digits=4)) +/- $(round(err(mDs_sim[ens.id][i]),digits=4)) \\\\ $(round(value(mymDs_sim[ens.id][i]),digits=4)) +/- $(round(err(mymDs_sim[ens.id][i]),digits=4)) \\\\ $(round(sigma[i],digits=2)) \\end{tabular}"
        row *= cell * " & "
    end

    difference = KappaC_tar[ens.id] - myKappaC_tar[ens.id]; uwerr(difference)
    sigma = abs(value(difference))/err(difference)
    cell = "\\begin{tabular}{@{}c@{}} $(round(value(KappaC_tar[ens.id]),digits=6)) +/- $(round(err(KappaC_tar[ens.id]),digits=6)) \\\\ $(round(value(myKappaC_tar[ens.id]),digits=6)) +/- $(round(err(myKappaC_tar[ens.id]),digits=6)) \\\\ $(round(sigma,digits=2)) \\end{tabular}"
    row *= cell

    row = row * " \\\\\n \\hline \n"
    latex_table *= row
end

# Close the LaTeX table
latex_table *= """
\\end{tabular}
\\caption{Comparison with 2019: results strange-heavy data for kappaC}
\\end{table}
"""


# Print or save the LaTeX table
println(latex_table)



## Data from window paper


amu = Dict()
for ensid in keys(IsoVecRes)
    amu[ensid] = Dict()
    for impr_set = ["1old","2"]
        amu[ensid][impr_set] = Dict()
        for discr in ["ll","lc"]
            amu[ensid][impr_set]["g33_$discr"] = IsoVecRes[ensid][impr_set][discr]
            if ensid != "J501"
                amu[ensid][impr_set]["g88_$discr"] = ∆lsRes[ensid][impr_set][discr]
            end
        end
    end
end

ensList = []
for ensid in keys(IsoVecRes)
    ensid != "H200" ? push!(ensList,ensid) : nothing
end
ensInfo = EnsInfo.(ensList)

