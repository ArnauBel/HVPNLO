
using ADerrors
import ADerrors: err
using HVPobs

include("Const/Const.jl")


const sqrtt0_ph_Madrid     = uwreal([0.1442, 0.0007], "sqrtt0 [fm] (Madrid)")
const sqrtt0_ph_Regensburg = uwreal([0.1449, 0.0007], "sqrtt0 [fm] (Regensburg)")  # Regensburg
const sqrtt0_ph_DESY       = uwreal([0.1443, 0.0014], "sqrtt0 [fm] (CLS)")  # 2206.06582v2 (ID window paper)
const sqrtt0_ph_Bruno      = uwreal([0.1460, 0.0019], "sqrtt0 [fm] (Bruno)")

const fPi_ph_PDGFLAG = uwreal([130.56,0.14]*1e-3, "fPi [GeV] PDG-FLAG")  # PDG and FLAG value

const artificial_err = 1e-8 # 1e-10
const fPi_ph    = uwreal([fPi_ph_PDGFLAG.mean,artificial_err], "fPi [GeV]")  # PDG and FLAG value
const sqrtt0_ph = uwreal([sqrtt0_ph_Madrid.mean,artificial_err], "sqrtt0 [fm]")  # Regengsburg scale choice
# const sqrtt0_ph = uwreal([sqrtt0_ph_Regensburg.mean,artificial_err], "sqrtt0 [fm]")  # Regengsburg scale choice
# const sqrtt0_ph = uwreal([sqrtt0_ph_DESY.mean,artificial_err], "sqrtt0 [fm]")  # CLS 2023 Window paper scale choice
const MD_ph     = uwreal([1968.47*1e-3,artificial_err], "MD_ph [GeV]")

const fK_ph     = uwreal([157.2,0.5]*1e-3, "fK [GeV]")  # PDG and FLAG value

const mpi_ph = uwreal([134.9768,0.0005],"mpi_ph")
const mK_ph  = uwreal([495.011,0.010],"mK_ph")

const phi2_ph = 8*(sqrtt0_ph*(mpi_ph*1e-3)/hbarc)^2
const phi4_ph = 8*sqrtt0_ph^2*(((mK_ph*1e-3)/hbarc)^2 + 0.5*((mpi_ph*1e-3)/hbarc)^2)

const y_ph = (mpi_ph*1e-3)^2/(8π*fPi_ph^2)
const z_ph = ((mK_ph*1e-3)^2+0.5*(mpi_ph*1e-3)^2)/(8π*(2/3*(fK_ph+0.5*fPi_ph))^2)


const b_values = [3.34, 3.40, 3.46, 3.55, 3.70, 3.85]

#2211.03744
const t0_ = [
    uwreal([2.204,0.005], "t0sym/a2 b=3.34"),
    uwreal([2.872,0.010], "t0sym/a2 b=3.40"),
    uwreal([3.682,0.012], "t0sym/a2 b=3.46"),
    uwreal([5.162,0.016], "t0sym/a2 b=3.55"),
    uwreal([8.613,0.025], "t0sym/a2 b=3.70"),
    uwreal([14.011,0.039],"t0sym/a2 b=3.85")]
const a_ = sqrtt0_ph ./ sqrt.(t0_)

# no Bruno data for beta = 3.34, 3.85 => values provided by Simon (Antoine)
const t0_Bruno = [
    uwreal([2.1715,0.0073], "t0sym/a2 b=3.34"),
    uwreal([2.860,0.011],   "t0sym/a2 b=3.40"),
    uwreal([3.659,0.016],   "t0sym/a2 b=3.46"),
    uwreal([5.164,0.018],   "t0sym/a2 b=3.55"),
    uwreal([8.595,0.025],   "t0sym/a2 b=3.70"),
    uwreal([13.99,0.067],   "t0sym/a2 b=3.85")]
const a_Bruno = sqrtt0_ph_Bruno ./ sqrt.(t0_Bruno)

const fPi_ = [
    uwreal([0.05894,0.00032],   "afPiph b=3.34")
    uwreal([0.05271,0.00023],   "afPiph b=3.40")
    uwreal([0.04741,0.00017],   "afPiph b=3.46")
    uwreal([0.04075,0.00013],   "afPiph b=3.55")
    uwreal([0.032035,0.000090], "afPiph b=3.70")
    uwreal([0.025374,0.000069], "afPiph b=3.85")
]

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


# const sqrtt0sym_over_sqrtt0_ph = uwreal([0.99974,0.00004],"sqrt(t0*)/sqrt(t0ph)")  # arXiv:2211.03744v2 [hep-lat] 4 May 2023 : Eq(3.42)


function fPiph(beta::Float64)
    fPi_[b_values .== beta][1]
end

const alpha_s_beta = [
    uwreal([0.29472,0.00019],   "BLO scal b=3.34"),
    uwreal([0.27571,0.00018],   "BLO scal b=3.40"),
    uwreal([0.26072,0.00017],   "BLO scal b=3.46"),
    uwreal([0.24279,0.00014],   "BLO scal b=3.55"),
    uwreal([0.22000,0.00010],   "BLO scal b=3.70"),
    uwreal([0.202337,0.000081], "BLO scal b=3.85")
]

function a2_rescaling(beta::Float64,Gamma::Float64=0.395;ERR::Bool=false)
    a2Res = alpha_s_beta[b_values .== beta][1]^Gamma
    ERR ? (return [value(a2Res),err(a2Res)]) : (return value(a2Res))
end

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

const m_ens = Dict(
    # pion and kaon masses extracted from 2206.06582v2
    # rho masses taken fron Dalibor (available also at 2203.08676v2)
    "A653" => Dict("mPi" => uwreal([0.21193,0.00091], "mPi"), "mK" => uwreal([0.21193,0.00091], "mK"), "mRho" => uwreal([0.4240, 0.0088], "mRho")), 
    "A654" => Dict("mPi" => uwreal([0.16647,0.00121], "mPi"), "mK" => uwreal([0.22712,0.00089], "mK"), "mRho" => uwreal([0.3988, 0.0019], "mRho")), 
    # H650 is new, computed by myself (mRho=uwreal([0.3846, 0.0200], "mRho") from interpolation from A654 and A653)
    "H650" => Dict("mPi" => uwreal([0.12721,0.00177], "mPi"), "mK" => uwreal([0.23098,0.00083], "mK"), "mRho" => uwreal([10.0,0.0], "mRho")), 
    
    "H101" => Dict("mPi" => uwreal([0.18217,0.00062], "mPi"), "mK" => uwreal([0.18217,0.00062], "mK"), "mRho" => uwreal([0.3709,0.0018], "mRho")),
    "H102" => Dict("mPi" => uwreal([0.15395,0.00071], "mPi"), "mK" => uwreal([0.19144,0.00057], "mK"), "mRho" => uwreal([0.3559,0.0036], "mRho")),
    "H105" => Dict("mPi" => uwreal([0.12136,0.00124], "mPi"), "mK" => uwreal([0.20230,0.00061], "mK"), "mRho" => uwreal([0.3468,0.0025], "mRho")),
    "N101" => Dict("mPi" => uwreal([0.12150,0.00055], "mPi"), "mK" => uwreal([0.20158,0.00031], "mK"), "mRho" => uwreal([0.3368,0.0045], "mRho")),
    "C101" => Dict("mPi" => uwreal([0.09569,0.00073], "mPi"), "mK" => uwreal([0.20579,0.00034], "mK"), "mRho" => uwreal([0.3262,0.0022], "mRho")),
    "C102" => Dict("mPi" => uwreal([0.09671,0.00078], "mPi"), "mK" => uwreal([0.21761,0.00047], "mK"), "mRho" => uwreal([0.3308,0.0038], "mRho")),   
    "D150" => Dict("mPi" => uwreal([0.05654,0.00094], "mPi"), "mK" => uwreal([0.20835,0.00035], "mK"), "mRho" => uwreal([0.3198,0.0026], "mRho")), 
    
    "B450" => Dict("mPi" => uwreal([0.16063,0.00045], "mPi"), "mK" => uwreal([0.16063,0.00045], "mK"), "mRho" => uwreal([0.3336,0.0013], "mRho")),
    "S400" => Dict("mPi" => uwreal([0.13506,0.00044], "mPi"), "mK" => uwreal([0.17022,0.00039], "mK"), "mRho" => uwreal([0.3094,0.0021], "mRho")),
    "N452" => Dict("mPi" => uwreal([0.13546,0.00030], "mPi"), "mK" => uwreal([0.17031,0.00026], "mK"), "mRho" => uwreal([0.3163,0.0023], "mRho")),
    "N451" => Dict("mPi" => uwreal([0.11072,0.00029], "mPi"), "mK" => uwreal([0.17824,0.00018], "mK"), "mRho" => uwreal([0.3071,0.0019], "mRho")),
    "D450" => Dict("mPi" => uwreal([0.08329,0.00043], "mPi"), "mK" => uwreal([0.18384,0.00018], "mK"), "mRho" => uwreal([0.2934,0.0017], "mRho")),
    "D451" => Dict("mPi" => uwreal([0.08359,0.00030], "mPi"), "mK" => uwreal([0.19402,0.00014], "mK"), "mRho" => uwreal([0.2942,0.0017], "mRho")),  
    "D452" => Dict("mPi" => uwreal([0.05941,0.00055], "mPi"), "mK" => uwreal([0.18651,0.00015], "mK"), "mRho" => uwreal([0.2855,0.0016], "mRho")), 
    
    "H200" => Dict("mPi" => uwreal([0.13535,0.00060], "mPi"), "mK" => uwreal([0.13535,0.00060], "mK"), "mRho" => uwreal([0.2878,0.0019], "mRho")),
    "N202" => Dict("mPi" => uwreal([0.13424,0.00031], "mPi"), "mK" => uwreal([0.13424,0.00031], "mK"), "mRho" => uwreal([0.2808,0.0015], "mRho")),
    "N203" => Dict("mPi" => uwreal([0.11254,0.00024], "mPi"), "mK" => uwreal([0.14402,0.00020], "mK"), "mRho" => uwreal([0.2684,0.0014], "mRho")),
    "N200" => Dict("mPi" => uwreal([0.09234,0.00031], "mPi"), "mK" => uwreal([0.15071,0.00023], "mK"), "mRho" => uwreal([0.2591,0.0023], "mRho")),
    "D251" => Dict("mPi" => uwreal([0.09293,0.00016], "mPi"), "mK" => uwreal([0.15041,0.00012], "mK"), "mRho" => uwreal([0.2584,0.0015], "mRho")), 
    "D200" => Dict("mPi" => uwreal([0.06507,0.00028], "mPi"), "mK" => uwreal([0.15630,0.00015], "mK"), "mRho" => uwreal([0.2392,0.0051], "mRho")), 
    "D201" => Dict("mPi" => uwreal([0.06499,0.00043], "mPi"), "mK" => uwreal([0.16309,0.00024], "mK"), "mRho" => uwreal([0.2477,0.0036], "mRho")), 
    "E250" => Dict("mPi" => uwreal([0.04170,0.00041], "mPi"), "mK" => uwreal([0.15924,0.00009], "mK"), "mRho" => uwreal([0.2362,0.0019], "mRho")), 
    
    "N300" => Dict("mPi" => uwreal([0.10569,0.00023], "mPi"), "mK" => uwreal([0.10569,0.00023], "mK"), "mRho" => uwreal([0.2302,0.0022], "mRho")),
    "J307" => Dict("mPi" => uwreal([0.10547,0.00042], "mPi"), "mK" => uwreal([0.10547,0.00042], "mK"), "mRho" => uwreal([0.2156,0.0023], "mRho")),
    "N302" => Dict("mPi" => uwreal([0.08690,0.00034], "mPi"), "mK" => uwreal([0.11358,0.00028], "mK"), "mRho" => uwreal([0.2172,0.0011], "mRho")),
    "J306" => Dict("mPi" => uwreal([0.08690,0.00019], "mPi"), "mK" => uwreal([0.11335,0.00019], "mK"), "mRho" => uwreal([0.2107,0.0022], "mRho")),
    "J303" => Dict("mPi" => uwreal([0.06475,0.00018], "mPi"), "mK" => uwreal([0.11963,0.00016], "mK"), "mRho" => uwreal([0.2016,0.0015], "mRho")),
    "J304" => Dict("mPi" => uwreal([0.06550,0.00020], "mPi"), "mK" => uwreal([0.13187,0.00017], "mK"), "mRho" => uwreal([0.2020,0.0015], "mRho")), 
    "E300" => Dict("mPi" => uwreal([0.04393,0.00016], "mPi"), "mK" => uwreal([0.12372,0.00010], "mK"), "mRho" => uwreal([0.1927,0.0010], "mRho")),
    "F300" => Dict("mPi" => uwreal([0.03381,0.00023], "mPi"), "mK" => uwreal([0.12358,0.00017], "mK"), "mRho" => uwreal([0.1849,0.0038], "mRho")),
    
    "J500" => Dict("mPi" => uwreal([0.08153,0.00019], "mPi"), "mK" => uwreal([0.08153,0.00019], "mK"), "mRho" => uwreal([0.1701,0.0011], "mRho")), 
    "J501" => Dict("mPi" => uwreal([0.06582,0.00023], "mPi"), "mK" => uwreal([0.08794,0.00022], "mK"), "mRho" => uwreal([0.1664,0.0011], "mRho"))  
)

const meson_ens = Dict(
    # extracted from 2411.07969v1 through ChatGPT
    "A653" => Dict("mPi" => uwreal([0.21184,0.00105], "mPi"), "mK" => uwreal([0.21184,0.00105], "mK"), "fPi" => uwreal([0.07144,0.00025], "fPi")),
    "A654" => Dict("mPi" => uwreal([0.16633,0.00131], "mPi"), "mK" => uwreal([0.22727,0.00112], "mK"), "fPi" => uwreal([0.06725,0.00025], "fPi")),
    # H650 is new, computed by myself
    "H650" => Dict("mPi" => uwreal([0.12721,0.00177], "mPi"), "mK" => uwreal([0.23098,0.00083], "mK"), "fPi" => uwreal([0.06396,0.00058], "fPi")),
    
    "H101" => Dict("mPi" => uwreal([0.18250,0.00071], "mPi"), "mK" => uwreal([0.18250,0.00071], "mK"), "fPi" => uwreal([0.06364,0.00030], "fPi")),
    "H102" => Dict("mPi" => uwreal([0.15383,0.00080], "mPi"), "mK" => uwreal([0.19135,0.00071], "mK"), "fPi" => uwreal([0.06044,0.00031], "fPi")),
    "H105" => Dict("mPi" => uwreal([0.12155,0.00115], "mPi"), "mK" => uwreal([0.20223,0.00085], "mK"), "fPi" => uwreal([0.05781,0.00083], "fPi")),
    "N101" => Dict("mPi" => uwreal([0.12120,0.00056], "mPi"), "mK" => uwreal([0.20146,0.00035], "mK"), "fPi" => uwreal([0.05773,0.00035], "fPi")),
    "C101" => Dict("mPi" => uwreal([0.09570,0.00078], "mPi"), "mK" => uwreal([0.20584,0.00044], "mK"), "fPi" => uwreal([0.05511,0.00036], "fPi")),
    "C102" => Dict("mPi" => uwreal([0.09640,0.00087], "mPi"), "mK" => uwreal([0.21766,0.00050], "mK"), "fPi" => uwreal([0.05507,0.00046], "fPi")),
    "D150" => Dict("mPi" => uwreal([0.05654,0.00094], "mPi"), "mK" => uwreal([0.20835,0.00035], "mK"), "fPi" => uwreal([0.05222,0.00030], "fPi")),
    
    "B450" => Dict("mPi" => uwreal([0.16081,0.00050], "mPi"), "mK" => uwreal([0.16081,0.00050], "mK"), "fPi" => uwreal([0.05685,0.00020], "fPi")),
    "S400" => Dict("mPi" => uwreal([0.13503,0.00046], "mPi"), "mK" => uwreal([0.17022,0.00041], "mK"), "fPi" => uwreal([0.05399,0.00034], "fPi")),
    "N452" => Dict("mPi" => uwreal([0.13546,0.00030], "mPi"), "mK" => uwreal([0.17031,0.00026], "mK"), "fPi" => uwreal([0.05462,0.00015], "fPi")),
    "N451" => Dict("mPi" => uwreal([0.11064,0.00045], "mPi"), "mK" => uwreal([0.17822,0.00026], "mK"), "fPi" => uwreal([0.05229,0.00015], "fPi")),
    "D450" => Dict("mPi" => uwreal([0.08346,0.00051], "mPi"), "mK" => uwreal([0.18393,0.00026], "mK"), "fPi" => uwreal([0.04977,0.00014], "fPi")),
    "D451" => Dict("mPi" => uwreal([0.08338,0.00035], "mPi"), "mK" => uwreal([0.19382,0.00016], "mK"), "fPi" => uwreal([0.05000,0.00024], "fPi")),
    "D452" => Dict("mPi" => uwreal([0.05932,0.00059], "mPi"), "mK" => uwreal([0.18645,0.00018], "mK"), "fPi" => uwreal([0.04758,0.00015], "fPi")),
    
    "H200" => Dict("mPi" => uwreal([0.13625,0.00064], "mPi"), "mK" => uwreal([0.13625,0.00064], "mK"), "fPi" => uwreal([0.04775,0.00034], "fPi")),
    "N202" => Dict("mPi" => uwreal([0.13436,0.00032], "mPi"), "mK" => uwreal([0.13436,0.00032], "mK"), "fPi" => uwreal([0.04846,0.00013], "fPi")),
    "N203" => Dict("mPi" => uwreal([0.11249,0.00027], "mPi"), "mK" => uwreal([0.14395,0.00023], "mK"), "fPi" => uwreal([0.04643,0.00016], "fPi")),
    "N200" => Dict("mPi" => uwreal([0.09221,0.00029], "mPi"), "mK" => uwreal([0.15065,0.00024], "mK"), "fPi" => uwreal([0.04420,0.00018], "fPi")),
    "D251" => Dict("mPi" => uwreal([0.09203,0.00016], "mPi"), "mK" => uwreal([0.15041,0.00012], "mK"), "fPi" => uwreal([0.04461,0.00010], "fPi")),
    "D200" => Dict("mPi" => uwreal([0.06502,0.00028], "mPi"), "mK" => uwreal([0.15630,0.00017], "mK"), "fPi" => uwreal([0.04237,0.00020], "fPi")),
    "D201" => Dict("mPi" => uwreal([0.06498,0.00043], "mPi"), "mK" => uwreal([0.16308,0.00024], "mK"), "fPi" => uwreal([0.04263,0.00025], "fPi")),
    "E250" => Dict("mPi" => uwreal([0.04232,0.00023], "mPi"), "mK" => uwreal([0.15936,0.00008], "mK"), "fPi" => uwreal([0.04018,0.00012], "fPi")),
    
    "N300" => Dict("mPi" => uwreal([0.10574,0.00030], "mPi"), "mK" => uwreal([0.10574,0.00030], "mK"), "fPi" => uwreal([0.03817,0.00016], "fPi")),
    "J307" => Dict("mPi" => uwreal([0.10547,0.00042], "mPi"), "mK" => uwreal([0.10547,0.00042], "mK"), "fPi" => uwreal([0.03785,0.00017], "fPi")),
    "N302" => Dict("mPi" => uwreal([0.08707,0.00054], "mPi"), "mK" => uwreal([0.11363,0.00046], "mK"), "fPi" => uwreal([0.03658,0.00021], "fPi")),
    "J306" => Dict("mPi" => uwreal([0.08690,0.00019], "mPi"), "mK" => uwreal([0.11335,0.00019], "mK"), "fPi" => uwreal([0.03653,0.00013], "fPi")),
    "J303" => Dict("mPi" => uwreal([0.06467,0.00022], "mPi"), "mK" => uwreal([0.11963,0.00019], "mK"), "fPi" => uwreal([0.03439,0.00015], "fPi")),
    "J304" => Dict("mPi" => uwreal([0.06561,0.00020], "mPi"), "mK" => uwreal([0.13187,0.00017], "mK"), "fPi" => uwreal([0.03418,0.00012], "fPi")),
    "E300" => Dict("mPi" => uwreal([0.04399,0.00012], "mPi"), "mK" => uwreal([0.12402,0.00009], "mK"), "fPi" => uwreal([0.03264,0.00012], "fPi")),
    "F300" => Dict("mPi" => uwreal([0.03381,0.00023], "mPi"), "mK" => uwreal([0.12358,0.00017], "mK"), "fPi" => uwreal([0.03168,0.00023], "fPi")),
    
    "J500" => Dict("mPi" => uwreal([0.08157,0.00017], "mPi"), "mK" => uwreal([0.08157,0.00017], "mK"), "fPi" => uwreal([0.02983,0.00010], "fPi")),
    "J501" => Dict("mPi" => uwreal([0.06590,0.00023], "mPi"), "mK" => uwreal([0.08796,0.00024], "mK"), "fPi" => uwreal([0.02855,0.00015], "fPi"))
)




# const E0_ens = Dict(
#     "A653" => Dict("fit_range" => [13,23], "mRho" => uwreal([0.42467,0.00154],"mRho"), "E2Pi" => uwreal([0.6735,0.0013],  "E2Pi"), "E0" => uwreal([0.4288,0.0015],  "E0")),
#     "A654" => Dict("fit_range" => [12,23], "mRho" => uwreal([0.40084,0.00259],"mRho"), "E2Pi" => uwreal([0.6203,0.0014],  "E2Pi"), "E0" => uwreal([0.4092,0.0024],  "E0")),
#     "H101" => Dict("fit_range" => [14,31], "mRho" => uwreal([0.37086,0.00219],"mRho"), "E2Pi" => uwreal([0.53613,0.00097],"E2Pi"), "E0" => uwreal([0.3735,0.0023],  "E0")),
#     "H102" => Dict("fit_range" => [17,28], "mRho" => uwreal([0.35029,0.00578],"mRho"), "E2Pi" => uwreal([0.49887,0.00099],"E2Pi"), "E0" => uwreal([0.3549,0.0055],  "E0")),
#     "U101" => Dict("fit_range" => [11,19], "mRho" => uwreal([0.36025,0.00777],"mRho"), "E2Pi" => uwreal([0.5711,0.0020],  "E2Pi"), "E0" => uwreal([0.3692,0.0068],  "E0")),
#     "H105" => Dict("fit_range" => [12,20], "mRho" => uwreal([0.34482,0.00805],"mRho"), "E2Pi" => uwreal([0.4619,0.0012],  "E2Pi"), "E0" => uwreal([0.3473,0.0071],  "E0")),
#     "N101" => Dict("fit_range" => [15,24], "mRho" => uwreal([0.35359,0.00570],"mRho"), "E2Pi" => uwreal([0.35679,0.00076],"E2Pi"), "E0" => uwreal([0.3389,0.0029],  "E0")),
#     "C101" => Dict("fit_range" => [15,31], "mRho" => uwreal([0.32556,0.00323],"mRho"), "E2Pi" => uwreal([0.32430,0.00092],"E2Pi"), "E0" => uwreal([0.3049,0.0018],  "E0")),
#     "D150" => Dict("fit_range" => [15,31], "mRho" => uwreal([0.31976,0.00330],"mRho"), "E2Pi" => uwreal([0.22659,0.00094],"E2Pi"), "E0" => uwreal([0.2237,0.0010],  "E0")),
#     "B450" => Dict("fit_range" => [14,31], "mRho" => uwreal([0.33387,0.00181],"mRho"), "E2Pi" => uwreal([0.50760,0.00064],"E2Pi"), "E0" => uwreal([0.3374,0.0018],  "E0")),
#     "H400" => Dict("fit_range" => [14,30], "mRho" => uwreal([0.33451,0.00270],"mRho"), "E2Pi" => uwreal([0.51057,0.00088],"E2Pi"), "E0" => uwreal([0.3378,0.0028],  "E0")),
#     "S400" => Dict("fit_range" => [15,32], "mRho" => uwreal([0.30843,0.00296],"mRho"), "E2Pi" => uwreal([0.47660,0.00052],"E2Pi"), "E0" => uwreal([0.3142,0.0029],  "E0")),
#     "N452" => Dict("fit_range" => [18,35], "mRho" => uwreal([0.31639,0.00252],"mRho"), "E2Pi" => uwreal([0.37674,0.00043],"E2Pi"), "E0" => uwreal([0.3173,0.0023],  "E0")),
#     "N451" => Dict("fit_range" => [14,31], "mRho" => uwreal([0.30698,0.00200],"mRho"), "E2Pi" => uwreal([0.34279,0.00059],"E2Pi"), "E0" => uwreal([0.3030,0.0015],  "E0")),
#     "D450" => Dict("fit_range" => [19,38], "mRho" => uwreal([0.29343,0.00201],"mRho"), "E2Pi" => uwreal([0.25771,0.00066],"E2Pi"), "E0" => uwreal([0.25264,0.00063],"E0")),
#     "D452" => Dict("fit_range" => [17,35], "mRho" => uwreal([0.28541,0.00204],"mRho"), "E2Pi" => uwreal([0.22941,0.00061],"E2Pi"), "E0" => uwreal([0.22468,0.00073],"E0")),
#     "H200" => Dict("fit_range" => [14,34], "mRho" => uwreal([0.28744,0.00249],"mRho"), "E2Pi" => uwreal([0.47799,0.00073],"E2Pi"), "E0" => uwreal([0.2921,0.0026],  "E0")),
#     "N202" => Dict("fit_range" => [17,37], "mRho" => uwreal([0.28073,0.00175],"mRho"), "E2Pi" => uwreal([0.37516,0.00045],"E2Pi"), "E0" => uwreal([0.2824,0.0018],  "E0")),
#     "N203" => Dict("fit_range" => [17,38], "mRho" => uwreal([0.26810,0.00215],"mRho"), "E2Pi" => uwreal([0.34519,0.00036],"E2Pi"), "E0" => uwreal([0.2700,0.0020],  "E0")),
#     "N200" => Dict("fit_range" => [18,36], "mRho" => uwreal([0.25929,0.00286],"mRho"), "E2Pi" => uwreal([0.32023,0.00033],"E2Pi"), "E0" => uwreal([0.2591,0.0024],  "E0")),
#     "D251" => Dict("fit_range" => [22,42], "mRho" => uwreal([0.25895,0.00181],"mRho"), "E2Pi" => uwreal([0.26913,0.00022],"E2Pi"), "E0" => uwreal([0.2514,0.0011],  "E0")),
#     "D200" => Dict("fit_range" => [28,41], "mRho" => uwreal([0.23920,0.00526],"mRho"), "E2Pi" => uwreal([0.23551,0.00031],"E2Pi"), "E0" => uwreal([0.2206,0.0025],  "E0")),
#     "E250" => Dict("fit_range" => [23,45], "mRho" => uwreal([0.23625,0.00211],"mRho"), "E2Pi" => uwreal([0.15588,0.00025],"E2Pi"), "E0" => uwreal([0.15448,0.00028],"E0")),
#     "N300" => Dict("fit_range" => [20,43], "mRho" => uwreal([0.22665,0.00158],"mRho"), "E2Pi" => uwreal([0.33655,0.00038],"E2Pi"), "E0" => uwreal([0.2292,0.0016],  "E0")),
#     "J307" => Dict("fit_range" => [25,45], "mRho" => uwreal([0.21558,0.00265],"mRho"), "E2Pi" => uwreal([0.28818,0.00061],"E2Pi"), "E0" => uwreal([0.2166,0.0027],  "E0")),
#     "N302" => Dict("fit_range" => [18,43], "mRho" => uwreal([0.21657,0.00212],"mRho"), "E2Pi" => uwreal([0.31442,0.00059],"E2Pi"), "E0" => uwreal([0.2197,0.0020],  "E0")),
#     "J306" => Dict("fit_range" => [22,41], "mRho" => uwreal([0.21056,0.00268],"mRho"), "E2Pi" => uwreal([0.26222,0.00025],"E2Pi"), "E0" => uwreal([0.2115,0.0024],  "E0")),
#     "J303" => Dict("fit_range" => [24,47], "mRho" => uwreal([0.20181,0.00183],"mRho"), "E2Pi" => uwreal([0.23512,0.00025],"E2Pi"), "E0" => uwreal([0.1991,0.0014],  "E0")),
#     "E300" => Dict("fit_range" => [27,51], "mRho" => uwreal([0.19067,0.00184],"mRho"), "E2Pi" => uwreal([0.15772,0.00014],"E2Pi"), "E0" => uwreal([0.15450,0.00034],"E0")),
#     "F300" => Dict("fit_range" => [33,51], "mRho" => uwreal([0.18503,0.00392],"mRho"), "E2Pi" => uwreal([0.11921,0.00026],"E2Pi"), "E0" => uwreal([0.11833,0.00068],"E0")),
#     "J500" => Dict("fit_range" => [32,58], "mRho" => uwreal([0.17019,0.00176],"mRho"), "E2Pi" => uwreal([0.25528,0.00022],"E2Pi"), "E0" => uwreal([0.1720,0.0018],  "E0")),
#     "J501" => Dict("fit_range" => [25,55], "mRho" => uwreal([0.16598,0.00186],"mRho"), "E2Pi" => uwreal([0.23649,0.00026],"E2Pi"), "E0" => uwreal([0.1682,0.0017],  "E0")),
#     "C102" => Dict("fit_range" => [13,25], "mRho" => uwreal([0.33063,0.00408],"mRho"), "E2Pi" => uwreal([0.3251,0.0010],  "E2Pi"), "E0" => uwreal([0.3077,0.0017],  "E0")),
#     "D451" => Dict("fit_range" => [19,35], "mRho" => uwreal([0.29622,0.00266],"mRho"), "E2Pi" => uwreal([0.25761,0.00046],"E2Pi"), "E0" => uwreal([0.25289,0.00076],"E0")),
#     "D201" => Dict("fit_range" => [24,38], "mRho" => uwreal([0.24816,0.00391],"mRho"), "E2Pi" => uwreal([0.23546,0.00048],"E2Pi"), "E0" => uwreal([0.2244,0.0015],  "E0")),
#     "J304" => Dict("fit_range" => [23,46], "mRho" => uwreal([0.20184,0.00173],"mRho"), "E2Pi" => uwreal([0.23616,0.00022],"E2Pi"), "E0" => uwreal([0.1994,0.0014],  "E0")),
#     "A650" => Dict("fit_range" => [16,23], "mRho" => uwreal([0.39003,0.00507],"mRho"), "E2Pi" => uwreal([0.6380,0.0016],  "E2Pi"), "E0" => uwreal([0.3976,0.0051],  "E0"))
# )