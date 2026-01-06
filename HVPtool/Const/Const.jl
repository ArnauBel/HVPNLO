using ADerrors

function custom_log(t::Union{Int64,Float64,uwreal})
    if typeof(t) == uwreal
        value(t)  == 0 ? (0) : (log(t))
    else
        t == 0 ? (0) : (log(abs(t)))
    end
end

const hbarc = 0.1973269804  # GeV fm

const masse   = 0.5109989461e-3     # electron mass (GeV) from PDG
const massmu  = 0.10565837          # muon mass (GeV) from PDG
const masstau = 1.77693             # tauon mass (GeV) from PDG

const alpha = 1/137.035999084     # em coupling q = 0
const GammaEuler = 0.57721566490153286060651209008240243

const Qlist = [3.5, 4.0, 5.0, 6.0, 7.0, 8.0]

M = masse/massmu
const C4 = Dict(
    "LO"     => (t -> 1/9),
    "NLOa"   => (t -> 317/324 + (23)/27 * (GammaEuler + custom_log(t)) - (2*π^2)/9),
    "NLOb"   => (t -> -1/27 * (1 + 4*log(M)) + M^2 * (2/3) + M^3 * (-4*π^2/9) + M^4 * (4/9 * (6 + π^2 + 6*log(M)^2)) + M^6 * (-(4/81) * (23 + 6*π^2 - 42*log(M) + 36*log(M)^2))),
    "NLOa&b" => (t -> C4["NLOa"](t) + C4["NLOb"](t))
)

const NLOab_zero = 3.5826239251814007

const kcd_in = Dict(
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

const p0_smallpbc_dict = Dict(
    "A653_set1_g33_ll" => 15:20,
    "A653_set1_g33_lc" => 15:20,
    "A653_set2_g33_ll" => 22:22,
    "A653_set2_g33_lc" => 15:16,

    "A654_set1_g33_ll" => 16:20,
    "A654_set1_g33_lc" => 15:17,
    "A654_set2_g33_ll" => 13:15,
    "A654_set2_g33_lc" => 14:16,

    "B450_set1_g33_ll" => 20:24,
    "B450_set1_g33_lc" => 20:24,
    "B450_set2_g33_ll" => 20:24,
    "B450_set2_g33_lc" => 20:24,
)

# the first two rows are no longer used but extracted with the mDs' interpolation

wpmm = Dict{String, Vector{Float64}}()
wpmm["H101"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["H102r002"] = [5.0, -2.0, -1.0, -1.0]
wpmm["H400"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["N202"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N200"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N203"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["N300"]     = [5.0, -1.5, -1.0, -1.0]
wpmm["J303"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J304"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["F300"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J306"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J307"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["E300"]     = [5.0, -1.0, -1.0, -1.0]
wpmm["F300"]     = [10.0, -1.0, -1.0, -1.0]