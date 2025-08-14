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
"N452" => uwreal([1.11412, 0.00058], "Zvc"), # TO BE UPDATED !
"S400" => uwreal([1.11159, 0.00089], "Zvc"),
"N451" => uwreal([1.11412, 0.00058], "Zvc"),
"D450" => uwreal([1.10790, 0.00026], "Zvc"),
"D451" => uwreal([1.10790, 0.00026], "Zvc"), # TO BE UPDATED!
"D452" => uwreal([1.10790, 0.00023], "Zvc"),

"N202" => uwreal([1.04843, 0.00085], "Zvc"), 
"N203" => uwreal([1.04534, 0.00039], "Zvc"),
"N200" => uwreal([1.04012, 0.00025], "Zvc"),
"D251" => uwreal([1.04012, 0.00045], "Zvc"), # TO BE UPDATED! 
"D200" => uwreal([1.03587, 0.00022], "Zvc"),
"D201" => uwreal([1.03587, 0.00022], "Zvc"), # TO BE UPDATED!
"E250" => uwreal([1.03310, 0.00011], "Zvc"),

"J307" => uwreal([0.96037, 0.00039], "Zvc"), # TO BE UPDATED!
"N300" => uwreal([0.97722, 0.00060], "Zvc"),
"N302" => uwreal([0.97241, 0.00030], "Zvc"),
"J306" => uwreal([0.96037, 0.00039], "Zvc"), # TO BE UPDATED!
"J303" => uwreal([0.96037, 0.00039], "Zvc"),
"J304" => uwreal([0.96037, 0.00039], "Zvc"), # TO BE UPDATED!
"E300" => uwreal([0.96639, 0.00026], "Zvc"),
"F300" => uwreal([0.96639, 0.00026], "Zvc"), # TO BE UPDATED!

"J500" => uwreal([0.93412, 0.00051], "Zvc"),
"J501" => uwreal([0.93412, 0.00051], "Zvc") # TO BE UPDATED!
)

m_ens = Dict(
    # pion and kaon masses extracted from 2206.06582v2
    # rho masses taken fron Dalibor (available also at 2203.08676v2)
    "A653" => Dict("m_pi" => uwreal([0.21193,0.00091], "mpi"), "m_K" => uwreal([0.21193,0.00091], "mK"), "m_rho" => uwreal([0.4240, 0.0088],    "mrho")), 
    "A654" => Dict("m_pi" => uwreal([0.16647,0.00121], "mpi"), "m_K" => uwreal([0.22712,0.00089], "mK"), "m_rho" => uwreal([0.3988, 0.0019],    "mrho")), 
    
    "H101" => Dict("m_pi" => uwreal([0.18217,0.00062], "mpi"), "m_K" => uwreal([0.18217,0.00062], "mK"), "m_rho" => uwreal([0.3709,0.0018], "mrho")),
    "H102" => Dict("m_pi" => uwreal([0.15395,0.00071], "mpi"), "m_K" => uwreal([0.19144,0.00057], "mK"), "m_rho" => uwreal([0.3559,0.0036], "mrho")),
    "H105" => Dict("m_pi" => uwreal([0.12136,0.00124], "mpi"), "m_K" => uwreal([0.20230,0.00061], "mK"), "m_rho" => uwreal([0.3468,0.0025], "mrho")),
    "N101" => Dict("m_pi" => uwreal([0.12150,0.00055], "mpi"), "m_K" => uwreal([0.20158,0.00031], "mK"), "m_rho" => uwreal([0.3368,0.0045], "mrho")),
    "C101" => Dict("m_pi" => uwreal([0.09569,0.00073], "mpi"), "m_K" => uwreal([0.20579,0.00034], "mK"), "m_rho" => uwreal([0.3262,0.0022], "mrho")),
    "C102" => Dict("m_pi" => uwreal([0.09671,0.00078], "mpi"), "m_K" => uwreal([0.21761,0.00047], "mK"), "m_rho" => uwreal([0.3308,0.0038], "mrho")),   
    "D150" => Dict("m_pi" => uwreal([0.05654,0.00094], "mpi"), "m_K" => uwreal([0.20835,0.00035], "mK"), "m_rho" => uwreal([0.3198,0.0026], "mrho")), 
    
    "B450" => Dict("m_pi" => uwreal([0.16063,0.00045], "mpi"), "m_K" => uwreal([0.16063,0.00045], "mK"), "m_rho" => uwreal([0.3336,0.0013], "mrho")),
    "S400" => Dict("m_pi" => uwreal([0.13506,0.00044], "mpi"), "m_K" => uwreal([0.17022,0.00039], "mK"), "m_rho" => uwreal([0.3094,0.0021], "mrho")),
    "N452" => Dict("m_pi" => uwreal([0.13546,0.00030], "mpi"), "m_K" => uwreal([0.17031,0.00026], "mK"), "m_rho" => uwreal([0.3163,0.0023], "mrho")),
    "N451" => Dict("m_pi" => uwreal([0.11072,0.00029], "mpi"), "m_K" => uwreal([0.17824,0.00018], "mK"), "m_rho" => uwreal([0.3071,0.0019], "mrho")),
    "D450" => Dict("m_pi" => uwreal([0.08329,0.00043], "mpi"), "m_K" => uwreal([0.18384,0.00018], "mK"), "m_rho" => uwreal([0.2934,0.0017], "mrho")),
    "D451" => Dict("m_pi" => uwreal([0.08359,0.00030], "mpi"), "m_K" => uwreal([0.19402,0.00014], "mK"), "m_rho" => uwreal([0.2942,0.0017], "mrho")),  
    "D452" => Dict("m_pi" => uwreal([0.05941,0.00055], "mpi"), "m_K" => uwreal([0.18651,0.00015], "mK"), "m_rho" => uwreal([0.2855,0.0016], "mrho")), 
    
    "H200" => Dict("m_pi" => uwreal([0.13535,0.00060], "mpi"), "m_K" => uwreal([0.13535,0.00060], "mK"), "m_rho" => uwreal([0.2878,0.0019], "mrho")),
    "N202" => Dict("m_pi" => uwreal([0.13424,0.00031], "mpi"), "m_K" => uwreal([0.13424,0.00031], "mK"), "m_rho" => uwreal([0.2808,0.0015], "mrho")),
    "N203" => Dict("m_pi" => uwreal([0.11254,0.00024], "mpi"), "m_K" => uwreal([0.14402,0.00020], "mK"), "m_rho" => uwreal([0.2684,0.0014], "mrho")),
    "N200" => Dict("m_pi" => uwreal([0.09234,0.00031], "mpi"), "m_K" => uwreal([0.15071,0.00023], "mK"), "m_rho" => uwreal([0.2591,0.0023], "mrho")),
    "D251" => Dict("m_pi" => uwreal([0.09293,0.00016], "mpi"), "m_K" => uwreal([0.15041,0.00012], "mK"), "m_rho" => uwreal([0.2584,0.0015], "mrho")), 
    "D200" => Dict("m_pi" => uwreal([0.06507,0.00028], "mpi"), "m_K" => uwreal([0.15630,0.00015], "mK"), "m_rho" => uwreal([0.2392,0.0051], "mrho")), 
    "D201" => Dict("m_pi" => uwreal([0.06499,0.00043], "mpi"), "m_K" => uwreal([0.16309,0.00024], "mK"), "m_rho" => uwreal([0.2477,0.0036], "mrho")), 
    "E250" => Dict("m_pi" => uwreal([0.04170,0.00041], "mpi"), "m_K" => uwreal([0.15924,0.00009], "mK"), "m_rho" => uwreal([0.2362,0.0019], "mrho")), 
    
    "N300" => Dict("m_pi" => uwreal([0.10569,0.00023], "mpi"), "m_K" => uwreal([0.10569,0.00023], "mK"), "m_rho" => uwreal([0.2302,0.0022], "mrho")),
    "J307" => Dict("m_pi" => uwreal([NaN, NaN], "mpi"), "m_K" => uwreal([NaN,NaN], "mK"), "m_rho" => uwreal([0.2156,0.0023], "mrho")),
    "N302" => Dict("m_pi" => uwreal([0.08690,0.00034], "mpi"), "m_K" => uwreal([0.11358,0.00028], "mK"), "m_rho" => uwreal([0.2172,0.0011], "mrho")),
    "J306" => Dict("m_pi" => uwreal([NaN,NaN], "mpi"), "m_K" => uwreal([NaN,NaN], "mK"), "m_rho" => uwreal([0.2107,0.0022], "mrho")),
    "J303" => Dict("m_pi" => uwreal([0.06475,0.00018], "mpi"), "m_K" => uwreal([0.11963,0.00016], "mK"), "m_rho" => uwreal([0.2016,0.0015], "mrho")),
    "J304" => Dict("m_pi" => uwreal([0.06550,0.00020], "mpi"), "m_K" => uwreal([0.13187,0.00017], "mK"), "m_rho" => uwreal([0.2020,0.0015], "mrho")), 
    "E300" => Dict("m_pi" => uwreal([0.04393,0.00016], "mpi"), "m_K" => uwreal([0.12372,0.00010], "mK"), "m_rho" => uwreal([0.1927,0.0010], "mrho")),
    "F300" => Dict("m_pi" => uwreal([NaN,NaN], "mpi"), "m_K" => uwreal([NaN,NaN], "mK"), "m_rho" => uwreal([0.1849,0.0038], "mrho")),
    
    "J500" => Dict("m_pi" => uwreal([0.08153,0.00019], "mpi"), "m_K" => uwreal([0.08153,0.00019], "mK"), "m_rho" => uwreal([0.1701,0.0011], "mrho")), 
    "J501" => Dict("m_pi" => uwreal([0.06582,0.00023], "mpi"), "m_K" => uwreal([0.08794,0.00022], "mK"), "m_rho" => uwreal([0.1646,0.0011], "mrho"))  
)

kcd_in = Dict( # values of kappa_charm. kappaC => interpolated value with error. kappaC_sim => 1st value used for HVP data. kappaC_sim_plus => second value used for HVP data
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

const b_values = [3.34, 3.40, 3.46, 3.55, 3.70, 3.85]
const hc = 197.3269804 #MeV fm

# scale setting
const t0sqrt_ph = uwreal([0.1449, 1e-8], "sqrtt0 [fm]")  # Regensburg with artificial error
const t0sqrt_ph_err = uwreal([0.1449, 0.0007], "sqrtt0 [fm] err")  # Regensburg with correct error

# const t0sqrt_ph_noerr = uwreal([0.1449, 0.0000], "sqrtt0 [fm]")  # Regensburg without error
# println("scale setting from BRUNO ET AL!!!")
const t0sqrt_ph_bruno = uwreal([0.1460, 0.0019], "sqrtt0 Bruno [fm]") 

#2211.03744 tab 5 first line
const t0_data = [2.204, 2.872, 3.682, 5.162, 8.613, 14.011]
const t0_error = [5, 10, 12, 16, 25, 39] .* 1e-3
# t0 sym from Bruno et al
const t0_data_bruno = [2.860, 3.659, 5.164, 8.595, 13.99]
const t0_data_bruno_err = [0.011, 0.016, 0.018, 0.025, 0.067]

const Ct0 = zeros(6, 6)
for i = 1:6
    Ct0[i,i] = t0_error[i] ^ 2    
end

const Ct0_bruno = zeros(5, 5)
for i = 1:5
    Ct0_bruno[i,i] = t0_data_bruno_err[i] ^ 2    
end

const t0_ = cobs(t0_data, Ct0, "t0sym/a2")
const t0_bruno_ = cobs(t0_data_bruno, Ct0_bruno, "t0sym/a2")

const a_ = t0sqrt_ph ./ sqrt.( t0_)

t0(beta::Float64) = t0_[b_values .== beta][1]
t0_bruno(beta::Float64) = t0_bruno_[b_values[2:end] .== beta][1]
a(beta::Float64)  = a_[b_values .== beta][1]

