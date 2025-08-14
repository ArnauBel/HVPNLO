using ALPHAio, ADerrors
import ADerrors: err

const hbarc = 0.1973269804  # GeV fm

const sqrtt0_ph_Regensburg = uwreal([0.1449, 0.0007], "sqrtt0 [fm] (Regensburg)")  # Regensburg
const sqrtt0_ph_Madrid     = uwreal([0.1439, 0.0006], "sqrtt0 [fm] (Madrid)")
const sqrtt0_ph_Bruno      = uwreal([0.1460,0.0019], "sqrtt0 [fm] (Bruno)")
# const sqrtt0_ph_CLS        = uwreal([0.1443, 0.0007], "sqrtt0 [fm] (CLS)")  # CLS
const sqrtt0_ph_CLS = uwreal([0.1443, 0.0014], "sqrtt0 [fm] (CLS)")  # 2206.06582v2 (ID window paper)

artificial_err = 1e-8 # 1e-10
const sqrtt0_ph = uwreal([value(sqrtt0_ph_Regensburg),artificial_err], "sqrtt0 [fm]")  # Regengsburg scale choice
# const sqrtt0_ph = uwreal([value(sqrtt0_ph_CLS),artificial_err], "sqrtt0 [fm]")  # CLS 2023 Window paper scale choice
const MD_ph     = uwreal([1968.47*1e-3,artificial_err], "MD_ph [GeV]")


const b_values = [3.34, 3.40, 3.46, 3.55, 3.70, 3.85]

#2211.03744
const t0_ = [
    uwreal([2.204,0.005], "t0sym/a2"),
    uwreal([2.872,0.010], "t0sym/a2"),
    uwreal([3.682,0.012], "t0sym/a2"),
    uwreal([5.162,0.016], "t0sym/a2"),
    uwreal([8.613,0.025], "t0sym/a2"),
    uwreal([14.011,0.039], "t0sym/a2")]
const a_ = sqrtt0_ph ./ sqrt.( t0_)

# no Bruno data for beta = 3.34, 3.85 => values provided by Simon (Antoine)
const t0_Bruno = [
    uwreal([2.1715,0.0073], "t0sym/a2"),
    uwreal([2.860,0.011], "t0sym/a2"),
    uwreal([3.659,0.016], "t0sym/a2"),
    uwreal([5.164,0.018], "t0sym/a2"),
    uwreal([8.595,0.025], "t0sym/a2"),
    uwreal([13.99,0.067], "t0sym/a2")]
const a_Bruno = sqrtt0_ph_Bruno ./ sqrt.( t0_Bruno)

function t0sym(beta::Float64; Bruno::Bool=false)
    if Bruno
        return t0_Bruno[b_values .== beta][1]
    else
        return t0_[b_values .== beta][1]
    end
end
function asym(beta::Float64; Bruno::Bool=false)
    if Bruno
        return a_Bruno[b_values .== beta][1]
    else
        return a_[b_values .== beta][1]
    end
end

const alpha_s_beta = [
    uwreal([0.29472,0.00019], "BLO scal"),
    uwreal([0.27571,0.00018], "BLO scal"),
    uwreal([0.26072,0.00017], "BLO scal"),
    uwreal([0.24279,0.00014], "BLO scal"),
    uwreal([0.22000,0.00010], "BLO scal"),
    uwreal([0.202337,0.000081], "BLO scal")
]

function a2_rescaling(beta::Float64,Gamma::Float64=0.395;ERR::Bool=false)
    a2Res = alpha_s_beta[b_values .== beta][1]^Gamma
    ERR ? (return [value(a2Res),err(a2Res)]) : (return value(a2Res))
end

const masse  = 0.5109989461e-3     # electron mass (GeV) from PDG
const massmu = 0.10565837         # muon mass (GeV) from PDG

const alpha = 1/137.035999084     #em coupling q = 0
const GammaEuler = 0.57721566490153286060651209008240243

const mpi_ph = uwreal([134.9768,0.0005],"mpi_ph")
const mK_ph  = uwreal([495.011,0.010],"mK_ph")

const phi2_ph = 8*(sqrtt0_ph*(mpi_ph*1e-3)/hbarc)^2
const phi4_ph = 8*sqrtt0_ph^2*(((mK_ph*1e-3)/hbarc)^2 + 0.5*((mpi_ph*1e-3)/hbarc)^2)


const Qlist = [3.5, 4.0, 5.0, 6.0, 7.0, 8.0]

M = masse/massmu
C4 = Dict(
    "LO"     => (t -> 1/9),
    "NLOa"   => (t -> 317/324 + (23)/27 * (GammaEuler + custom_log(t)) - (2*π^2)/9),
    "NLOb"   => (t -> -1/27 * (1 + 4*log(M)) + M^2 * (2/3) + M^3 * (-4*π^2/9) + M^4 * (4/9 * (6 + π^2 + 6*log(M)^2)) + M^6 * (-(4/81) * (23 + 6*π^2 - 42*log(M) + 36*log(M)^2))),
    "NLOa&b" => (t -> C4["NLOa"](t) + C4["NLOb"](t))
)

const Zvc_l = Dict(
"A653" => uwreal([1.32281, 0.00072], "Zvc"),
"A654" => uwreal([1.30495, 0.00106], "Zvc"),

"H101" => uwreal([1.20324, 0.00071], "Zvc"),
"H102" => uwreal([1.19743, 0.00100], "Zvc"),
"H105" => uwreal([1.18964, 0.00075], "Zvc"),
"N101" => uwreal([1.18964, 0.00075], "Zvc"),
"C101" => uwreal([1.18500, 0.00044], "Zvc"),
"C102" => uwreal([1.18500, 0.00044], "Zvc"), # TO BE UPDATED!
"D150" => uwreal([1.18500, 0.00044], "Zvc"), # TO BE UPDATED!

"B450" => uwreal([1.12972, 0.00083], "Zvc"),
"S400" => uwreal([1.11159, 0.00089], "Zvc"),
"N451" => uwreal([1.11412, 0.00058], "Zvc"),
"D450" => uwreal([1.10790, 0.00026], "Zvc"),
"D451" => uwreal([1.10790, 0.00026], "Zvc"), # TO BE UPDATED!
"D452" => uwreal([1.10790, 0.00023], "Zvc"),

"H200" => uwreal([1.04843, 0.00085], "Zvc"), 
"N202" => uwreal([1.04843, 0.00085], "Zvc"), 
"N203" => uwreal([1.04534, 0.00039], "Zvc"),
"N200" => uwreal([1.04012, 0.00025], "Zvc"),
"D251" => uwreal([1.04012, 0.00045], "Zvc"), # TO BE UPDATED! 
"D200" => uwreal([1.03587, 0.00022], "Zvc"),
"D201" => uwreal([1.03587, 0.00022], "Zvc"), # TO BE UPDATED!
"E250" => uwreal([1.03310, 0.00011], "Zvc"),

"N300" => uwreal([0.97722, 0.00060], "Zvc"),
"N302" => uwreal([0.97241, 0.00030], "Zvc"),
"J303" => uwreal([0.96037, 0.00039], "Zvc"),
"J304" => uwreal([0.96037, 0.00039], "Zvc"), # TO BE UPDATED!
"E300" => uwreal([0.96639, 0.00026], "Zvc"),

"J500" => uwreal([0.93412, 0.00051], "Zvc"),
"J501" => uwreal([0.93412, 0.00051], "Zvc") # TO BE UPDATED!
)

# pion and kaon masses extracted from 2401.11895v1
# rho masses extracted from 2203.08676v2
# m_ens = Dict(
#     "A653" => Dict("m_pi" => uwreal([0.21184, 0.00105], "mpi A653"), "m_K" => uwreal([0.21184, 0.00105], "mK A653"), "m_rho" => uwreal([NaN, NaN]   , "mrho A653")),
#     "A654" => Dict("m_pi" => uwreal([0.16633, 0.00131], "mpi A654"), "m_K" => uwreal([0.22727, 0.00112], "mK A654"), "m_rho" => uwreal([NaN, NaN]   , "mrho A654")),
#     "H101" => Dict("m_pi" => uwreal([0.18250, 0.00105], "mpi H101"), "m_K" => uwreal([0.18250, 0.00105], "mK H101"), "m_rho" => uwreal([0.375,0.002], "mrho H101")),
#     "H102" => Dict("m_pi" => uwreal([0.15383, 0.00080], "mpi H102"), "m_K" => uwreal([0.19135, 0.00071], "mK H102"), "m_rho" => uwreal([0.358,0.003], "mrho H102")),
#     "H105" => Dict("m_pi" => uwreal([0.12155, 0.00115], "mpi H105"), "m_K" => uwreal([0.20223, 0.00085], "mK H105"), "m_rho" => uwreal([0.338,0.011], "mrho H105")),
#     "N101" => Dict("m_pi" => uwreal([0.12120, 0.00056], "mpi N101"), "m_K" => uwreal([0.20146, 0.00035], "mK N101"), "m_rho" => uwreal([0.340,0.004], "mrho N101")),
#     "C101" => Dict("m_pi" => uwreal([0.09570, 0.00078], "mpi C101"), "m_K" => uwreal([0.20584, 0.00044], "mK C101"), "m_rho" => uwreal([0.326,0.003], "mrho C101")),
#     "C102" => Dict("m_pi" => uwreal([0.09640, 0.00087], "mpi C102"), "m_K" => uwreal([0.21766, 0.00050], "mK C102"), "m_rho" => uwreal([NaN, NaN]   , "mrho C102")),
#     "D150" => Dict("m_pi" => uwreal([0.05654, 0.00094], "mpi D150"), "m_K" => uwreal([0.20835, 0.00035], "mK D150"), "m_rho" => uwreal([NaN, NaN]   , "mrho D150")),
#     "B450" => Dict("m_pi" => uwreal([0.16081, 0.00050], "mpi B450"), "m_K" => uwreal([0.16081, 0.00050], "mK B450"), "m_rho" => uwreal([0.337,0.001], "mrho B450")),
#     "S400" => Dict("m_pi" => uwreal([0.13503, 0.00046], "mpi S400"), "m_K" => uwreal([0.17022, 0.00041], "mK S400"), "m_rho" => uwreal([0.312,0.004], "mrho S400")),
#     "N452" => Dict("m_pi" => uwreal([0.13546, 0.00030], "mpi N452"), "m_K" => uwreal([0.17031, 0.00026], "mK N452"), "m_rho" => uwreal([NaN, NaN]   , "mrho N452")),
#     "N451" => Dict("m_pi" => uwreal([0.11064, 0.00045], "mpi N451"), "m_K" => uwreal([0.17822, 0.00026], "mK N451"), "m_rho" => uwreal([0.302,0.004], "mrho N451")),
#     "D450" => Dict("m_pi" => uwreal([0.08346, 0.00051], "mpi D450"), "m_K" => uwreal([0.18393, 0.00026], "mK D450"), "m_rho" => uwreal([0.303,0.008], "mrho D450")),
#     "D451" => Dict("m_pi" => uwreal([0.08338, 0.00035], "mpi D451"), "m_K" => uwreal([0.19382, 0.00016], "mK D451"), "m_rho" => uwreal([NaN, NaN]   , "mrho D451")),
#     "D452" => Dict("m_pi" => uwreal([0.05932, 0.00059], "mpi D452"), "m_K" => uwreal([0.18645, 0.00018], "mK D452"), "m_rho" => uwreal([NaN, NaN]   , "mrho D452")),
#     "H200" => Dict("m_pi" => uwreal([0.13625, 0.00064], "mpi H200"), "m_K" => uwreal([0.13625, 0.00064], "mK H200"), "m_rho" => uwreal([0.286,0.003], "mrho H200")),
#     "N202" => Dict("m_pi" => uwreal([0.13436, 0.00032], "mpi N202"), "m_K" => uwreal([0.13436, 0.00032], "mK N202"), "m_rho" => uwreal([0.280,0.003], "mrho N202")),
#     "N203" => Dict("m_pi" => uwreal([0.11249, 0.00027], "mpi N203"), "m_K" => uwreal([0.14395, 0.00023], "mK N203"), "m_rho" => uwreal([0.268,0.001], "mrho N203")),
#     "N200" => Dict("m_pi" => uwreal([0.09221, 0.00029], "mpi N200"), "m_K" => uwreal([0.15065, 0.00024], "mK N200"), "m_rho" => uwreal([0.252,0.002], "mrho N200")),
#     "D251" => Dict("m_pi" => uwreal([0.09203, 0.00016], "mpi D251"), "m_K" => uwreal([0.15041, 0.00012], "mK D251"), "m_rho" => uwreal([NaN, NaN]   , "mrho D251")),
#     "D200" => Dict("m_pi" => uwreal([0.06502, 0.00028], "mpi D200"), "m_K" => uwreal([0.15630, 0.00017], "mK D200"), "m_rho" => uwreal([0.250,0.002], "mrho D200")),
#     "D201" => Dict("m_pi" => uwreal([0.06498, 0.00043], "mpi D201"), "m_K" => uwreal([0.16308, 0.00024], "mK D201"), "m_rho" => uwreal([NaN, NaN]   , "mrho D201")),
#     "E250" => Dict("m_pi" => uwreal([0.04232, 0.00033], "mpi E250"), "m_K" => uwreal([0.15936, 0.00008], "mK E250"), "m_rho" => uwreal([0.251,0.004], "mrho E250")),
#     "N300" => Dict("m_pi" => uwreal([0.10574, 0.00030], "mpi N300"), "m_K" => uwreal([0.10574, 0.00030], "mK N300"), "m_rho" => uwreal([0.222,0.003], "mrho N300")),
#     "J307" => Dict("m_pi" => uwreal([0.10547, 0.00042], "mpi J307"), "m_K" => uwreal([0.10547, 0.00042], "mK J307"), "m_rho" => uwreal([NaN, NaN]   , "mrho J307")),
#     "N302" => Dict("m_pi" => uwreal([0.08707, 0.00054], "mpi N302"), "m_K" => uwreal([0.11363, 0.00046], "mK N302"), "m_rho" => uwreal([0.216,0.003], "mrho N302")),
#     "J306" => Dict("m_pi" => uwreal([0.08690, 0.00019], "mpi J306"), "m_K" => uwreal([0.11335, 0.00019], "mK J306"), "m_rho" => uwreal([NaN, NaN]   , "mrho J306")),
#     "J303" => Dict("m_pi" => uwreal([0.06467, 0.00022], "mpi J303"), "m_K" => uwreal([0.11963, 0.00019], "mK J303"), "m_rho" => uwreal([0.200,0.002], "mrho J303")),
#     "J304" => Dict("m_pi" => uwreal([0.06561, 0.00020], "mpi J304"), "m_K" => uwreal([0.13187, 0.00017], "mK J304"), "m_rho" => uwreal([NaN, NaN]   , "mrho J304")),
#     "E300" => Dict("m_pi" => uwreal([0.04399, 0.00012], "mpi E300"), "m_K" => uwreal([0.12402, 0.00009], "mK E300"), "m_rho" => uwreal([0.198,0.002], "mrho E300")),
#     "F300" => Dict("m_pi" => uwreal([0.03381, 0.00023], "mpi F300"), "m_K" => uwreal([0.12358, 0.00017], "mK F300"), "m_rho" => uwreal([NaN, NaN]   , "mrho F300")),
#     "J500" => Dict("m_pi" => uwreal([0.08157, 0.00017], "mpi J500"), "m_K" => uwreal([0.08157, 0.00017], "mK J500"), "m_rho" => uwreal([NaN, NaN]   , "mrho J500")),
#     "J501" => Dict("m_pi" => uwreal([0.06590, 0.00023], "mpi J501"), "m_K" => uwreal([0.08796, 0.00024], "mK J501"), "m_rho" => uwreal([NaN, NaN]   , "mrho J501"))
# )

    # "F300" => Dict("m_pi" => uwreal([0.03381, 0.00023], "mpi F300"), "m_K" => uwreal([0.12358, 0.00017], "mK F300"), "m_rho" => uwreal([NaN, NaN]   , "mrho F300")),
    # "J306" => Dict("m_pi" => uwreal([0.08690, 0.00019], "mpi J306"), "m_K" => uwreal([0.11335, 0.00019], "mK J306"), "m_rho" => uwreal([NaN, NaN]   , "mrho J306")),
    # "J307" => Dict("m_pi" => uwreal([0.10547, 0.00042], "mpi J307"), "m_K" => uwreal([0.10547, 0.00042], "mK J307"), "m_rho" => uwreal([NaN, NaN]   , "mrho J307")),


m_ens = Dict(
    # pion and kaon masses extracted from 2206.06582v2
    # rho masses taken fron Dalibor (available also at 2203.08676v2)
    "A653" => Dict("m_pi" => uwreal([0.21193,0.00091], "mpi A653"), "m_K" => uwreal([0.21193,0.00091], "mK A653"), "m_rho" => uwreal([0.4240, 0.0088], "mrh A653o")), 
    "A654" => Dict("m_pi" => uwreal([0.16647,0.00121], "mpi A654"), "m_K" => uwreal([0.22712,0.00089], "mK A654"), "m_rho" => uwreal([0.3988, 0.0019], "mrh A654o")), 
                   
    "H101" => Dict("m_pi" => uwreal([0.18217,0.00062], "mpi H101"), "m_K" => uwreal([0.18217,0.00062], "mK H101"), "m_rho" => uwreal([0.3709,0.0018], "mrho H101")),
    "H102" => Dict("m_pi" => uwreal([0.15395,0.00071], "mpi H102"), "m_K" => uwreal([0.19144,0.00057], "mK H102"), "m_rho" => uwreal([0.3559,0.0036], "mrho H102")),
    "H105" => Dict("m_pi" => uwreal([0.12136,0.00124], "mpi H105"), "m_K" => uwreal([0.20230,0.00061], "mK H105"), "m_rho" => uwreal([0.3468,0.0025], "mrho H105")),
    "N101" => Dict("m_pi" => uwreal([0.12150,0.00055], "mpi N101"), "m_K" => uwreal([0.20158,0.00031], "mK N101"), "m_rho" => uwreal([0.3368,0.0045], "mrho N101")),
    "C101" => Dict("m_pi" => uwreal([0.09569,0.00073], "mpi C101"), "m_K" => uwreal([0.20579,0.00034], "mK C101"), "m_rho" => uwreal([0.3262,0.0022], "mrho C101")),
    "C102" => Dict("m_pi" => uwreal([0.09671,0.00078], "mpi C102"), "m_K" => uwreal([0.21761,0.00047], "mK C102"), "m_rho" => uwreal([0.3308,0.0038], "mrho C102")),   
    "D150" => Dict("m_pi" => uwreal([0.05654,0.00094], "mpi D150"), "m_K" => uwreal([0.20835,0.00035], "mK D150"), "m_rho" => uwreal([0.3198,0.0026], "mrho D150")), 
                   
    "B450" => Dict("m_pi" => uwreal([0.16063,0.00045], "mpi B450"), "m_K" => uwreal([0.16063,0.00045], "mK B450"), "m_rho" => uwreal([0.3336,0.0013], "mrho B450")),
    "S400" => Dict("m_pi" => uwreal([0.13506,0.00044], "mpi S400"), "m_K" => uwreal([0.17022,0.00039], "mK S400"), "m_rho" => uwreal([0.3094,0.0021], "mrho S400")),
    "N451" => Dict("m_pi" => uwreal([0.11072,0.00029], "mpi N451"), "m_K" => uwreal([0.17824,0.00018], "mK N451"), "m_rho" => uwreal([0.3071,0.0019], "mrho N451")),
    "N452" => Dict("m_pi" => uwreal([0.13546,0.00030], "mpi N452"), "m_K" => uwreal([0.17031,0.00026], "mK N452"), "m_rho" => uwreal([NaN,NaN]      , "mrho N452")), # from 2203.08676v2
    "D450" => Dict("m_pi" => uwreal([0.08329,0.00043], "mpi D450"), "m_K" => uwreal([0.18384,0.00018], "mK D450"), "m_rho" => uwreal([0.2934,0.0017], "mrho D450")),
    "D451" => Dict("m_pi" => uwreal([0.08359,0.00030], "mpi D451"), "m_K" => uwreal([0.19402,0.00014], "mK D451"), "m_rho" => uwreal([0.2942,0.0017], "mrho D451")),  
    "D452" => Dict("m_pi" => uwreal([0.05941,0.00055], "mpi D452"), "m_K" => uwreal([0.18651,0.00015], "mK D452"), "m_rho" => uwreal([0.2855,0.0016], "mrho D452")), 
                   
    "H200" => Dict("m_pi" => uwreal([0.13535,0.00060], "mpi H200"), "m_K" => uwreal([0.13535,0.00060], "mK H200"), "m_rho" => uwreal([0.2878,0.0019], "mrho H200")),
    "N202" => Dict("m_pi" => uwreal([0.13424,0.00031], "mpi N202"), "m_K" => uwreal([0.13424,0.00031], "mK N202"), "m_rho" => uwreal([0.2808,0.0015], "mrho N202")),
    "N203" => Dict("m_pi" => uwreal([0.11254,0.00024], "mpi N203"), "m_K" => uwreal([0.14402,0.00020], "mK N203"), "m_rho" => uwreal([0.2684,0.0014], "mrho N203")),
    "N200" => Dict("m_pi" => uwreal([0.09234,0.00031], "mpi N200"), "m_K" => uwreal([0.15071,0.00023], "mK N200"), "m_rho" => uwreal([0.2591,0.0023], "mrho N200")),
    "D251" => Dict("m_pi" => uwreal([0.09293,0.00016], "mpi D251"), "m_K" => uwreal([0.15041,0.00012], "mK D251"), "m_rho" => uwreal([0.2584,0.0015], "mrho D251")), 
    "D200" => Dict("m_pi" => uwreal([0.06507,0.00028], "mpi D200"), "m_K" => uwreal([0.15630,0.00015], "mK D200"), "m_rho" => uwreal([0.2392,0.0051], "mrho D200")), 
    "D201" => Dict("m_pi" => uwreal([0.06499,0.00043], "mpi D201"), "m_K" => uwreal([0.16309,0.00024], "mK D201"), "m_rho" => uwreal([0.2477,0.0036], "mrho D201")), 
    "E250" => Dict("m_pi" => uwreal([0.04170,0.00041], "mpi E250"), "m_K" => uwreal([0.15924,0.00009], "mK E250"), "m_rho" => uwreal([0.2362,0.0019], "mrho E250")), 
                   
    "N300" => Dict("m_pi" => uwreal([0.10569,0.00023], "mpi N300"), "m_K" => uwreal([0.10569,0.00023], "mK N300"), "m_rho" => uwreal([0.2302,0.0022], "mrho N300")),
    "J307" => Dict("m_pi" => uwreal([0.10547,0.00042], "mpi J307"), "m_K" => uwreal([0.10547,0.00042], "mK J307"), "m_rho" => uwreal([0.2156,0.0023], "mrho J307")), # from 2203.08676v2
    "N302" => Dict("m_pi" => uwreal([0.08690,0.00034], "mpi N302"), "m_K" => uwreal([0.11358,0.00028], "mK N302"), "m_rho" => uwreal([0.2172,0.0011], "mrho N302")),
    "J306" => Dict("m_pi" => uwreal([0.08690,0.00019], "mpi J306"), "m_K" => uwreal([0.11335,0.00019], "mK J306"), "m_rho" => uwreal([0.2107,0.0022], "mrho J306")), # from 2203.08676v2
    "J303" => Dict("m_pi" => uwreal([0.06475,0.00018], "mpi J303"), "m_K" => uwreal([0.11963,0.00016], "mK J303"), "m_rho" => uwreal([0.2016,0.0015], "mrho J303")),
    "J304" => Dict("m_pi" => uwreal([0.06550,0.00020], "mpi J304"), "m_K" => uwreal([0.13187,0.00017], "mK J304"), "m_rho" => uwreal([0.2020,0.0015], "mrho J304")), 
    "E300" => Dict("m_pi" => uwreal([0.04393,0.00016], "mpi E300"), "m_K" => uwreal([0.12372,0.00010], "mK E300"), "m_rho" => uwreal([0.1927,0.0010], "mrho E300")),
    "F300" => Dict("m_pi" => uwreal([0.03381,0.00023], "mpi F300"), "m_K" => uwreal([0.03381,0.00023], "mK F300"), "m_rho" => uwreal([0.1849,0.0038], "mrho F300")),
                   
    "J500" => Dict("m_pi" => uwreal([0.08153,0.00019], "mpi J500"), "m_K" => uwreal([0.08153,0.00019], "mK J500"), "m_rho" => uwreal([0.1701,0.0011], "mrho J500")), 
    "J501" => Dict("m_pi" => uwreal([0.06582,0.00023], "mpi J501"), "m_K" => uwreal([0.08794,0.00022], "mK J501"), "m_rho" => uwreal([0.1646,0.0011], "mrho J501"))  
)


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
# the first two rows are no longer used but extracted with the mDs' interpolation