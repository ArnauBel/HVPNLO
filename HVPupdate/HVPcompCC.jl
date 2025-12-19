# Import packages

using Revise

include("HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

using ProgressBars
using Suppressor
using TimerOutputs

# include uwreal constants

# include("HVPtool/uwConst.jl")

# BDIO path definition (set 'STD_DERIV = true' to use the standard sym. derivative in the impr.)

julia_script_directory = @__DIR__

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_tl   = joinpath(julia_script_directory,"..", "..", "HVPData", "tree_level_improv")
path_spec = joinpath(julia_script_directory,"..", "..", "HVPData", "spectroscopy")

pFVC_MLL  = joinpath(julia_script_directory, "..",  "..", "HVPData", "FSE_MLL")

# Blind analysis (Simon K.) safe ensembles: 
# SU(3) sym.  H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

# All considered ensembles are:
ensList = [
    "A653","A654","B450","C101","C102",
    "D150","D200","D201","D251","D450",
    "D451","D452","E250","E300","F300",
    "H101","H102","H200","H650","J303",
    "J304","J306","J307","J500","J501",
    "N101","N200","N202","N203","N300",
    "N302","N451","N452","S400"
    ]


ensInfo = EnsInfo.(ensList)

# We do not have charm or disconnected data for some of the ensembles
ensNOcharm = ["C102","D150","D201","D251","D451","F300","H200","H650","J304","J306","J307","J501","N451","N452"]
ensNOdisc  = ["F300","H650","J306"]

ensSPECdata = ["D200","E250"]  # J303

ensTAILfit  = ["A653","A654","B450"]

wpm = Dict{String, Vector{Float64}}()

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> 1D HVP computation (charm) [LO, NLOa, NLOb] <==========================##


diag = ""  # LO  NLOa  NLOb  NLOa&b
wind = ""  # NW  SD  ID  LD  ILD

IMPR_SET = ["1","2"] # ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

STD_DERIV = false
RESC      = false

BLIND = true

WRITE = true
# OVERWRITE = true  # Always true in here!! 


path_bdio = path_bdio_dict["local"]


BM = wind in ["NW","LD","ILD","LD2"]


@info("STARTING CHARM HVP COMPUTATION [diag. $diag; wind. $wind; resc. $RESC]")
STD_DERIV ? @info("STANDARD DERIVATIVE is being employed in the IMPROVEMENT") : nothing

# LO => (alpha/pi)^2 || NLO => (alpha/pi)^3
exp_diag = diag == "LO" ? 2 : 3

@time begin
    for ens in ensInfo

        @info("Computing HVP for ensemble $(ens.id)")
        ens.id ∈ ensNOdisc ? @info("  > NO DISC. DATA FOR $(ens.id)") : nothing
        if ens.id ∈ ensNOcharm
            @info("  > NO CHARM DATA FOR $(ens.id)\n    Skipping this ensemble")
            continue
        end

        println("   - Reading TMR...")
        TMRbeta = BDIOread_TMR(path_bdio,ens,diag,resc=RESC,beta=true)

        T = HVPobs.Data.get_T(ens.id)
        sym_points = Int64(T/2+1)
        t = collect(1:sym_points)

        if !RESC
            println("   - Reading t0...")
            # t0_sym = BDIOread_t0_SU3sym(path_bdio,ens)
            t0_sym = t0sym(ens.beta)
            a_ß = sqrtt0_ph / sqrt(t0_sym); tfm_ß = a_ß.*(collect(1:length(TMRbeta)).-1)
        else
            println("   - Reading fPi...")
            a_ß = hbarc * fPiph(ens.beta) / fPi_ph; tfm_ß = a_ß.*(collect(1:length(TMRbeta)).-1)
        end

        # define windowed kernel
        TMRwbeta = wind == "NW" ? TMRbeta : TMRbeta .* Window(wind)(tfm_ß)

        for impr_set in IMPR_SET

            println("   - Starting set "*impr_set)

            println("      - Reading HVP(!!)...")
            HVP, info = BDIOread_HVPens(path_bdio,diag,wind,ens.id,impr_set,info=true,resc=RESC,STD=STD_DERIV,BLIND=BLIND)

            println("      - Reading corr...")
            corr = BDIOread_corr(path_bdio,ens,impr_set,STD=STD_DERIV)

            println("      - Computing HVPs...")

            for key in ["gCCconn_ll_sim","gCCconn_lc_sim","gCCconn_ll_sim+","gCCconn_lc_sim+"]
                println("         - $key  &  $(key[1:8]*"SU3"*key[8:end])")

                obs = corr[key][t]

                # int      = obs .* TMRw[t]
                int_SU3  = obs .* TMRwbeta[t]
                # amu      = (alpha/pi)^exp_diag * sum(int) * 1e10
                amu_SU3  = (alpha/pi)^exp_diag * sum(int_SU3) * 1e10

                # HVP[key] = amu
                HVP[key[1:8]*"SU3"*key[8:end]] = amu_SU3
            end

            println("         - Reading kappaC target...")

            # kappaC_tar = BDIOread_KappaC_tar(path_bdio,ens)
            mDs_ph_prime, mDs_beta, Ds_dict = BDIOread_mDs_kappaC(path_bdio,ens)
            kappaC_tar = Ds_dict["kappaC"]

            println("         - Interpolating charm contributions...")

            # kappaC_tar = uwreal([kcd_in[ens.id]["kappaC"],kcd_in[ens.id]["kappaC_err"]], "kappaC target")   # extract KappaC from tables
            kappaC = [kcd_in[ens.id]["kappaC_sim"],kcd_in[ens.id]["kappaC_sim_plus"]]

            for key in ["gCCconn_SU3_ll","gCCconn_SU3_lc"] # ["gCCconn_ll","gCCconn_lc","gCCconn_SU3_ll","gCCconn_SU3_lc"]
                obs  = [HVP["$(key)_sim"],HVP["$(key)_sim+"]]
                
                par, _  = lin_fit(kappaC,obs,wpm=wpm,lineprint=false)
                obs_tar = y_lin_fit(par,kappaC_tar)

                if kappaC[1] <= value(kappaC_tar) <= kappaC[2]
                    HVP[key] = obs_tar
                    info["HVPsyst"][key] = 0.0
                else
                    HVP[key] = obs_tar

                    DkappaC = abs(kappaC[1]-kappaC[2])
                    DkappaC_tar = minimum([abs(value(kappaC_tar)-kappaC[1]),abs(value(kappaC_tar)-kappaC[2])])
                    uwerr(obs_tar); info["HVPsyst"][key] = 1.5 * (abs(DkappaC_tar) / DkappaC)^2 * err(obs_tar)
                end
            end

            # if wind == "SD" && ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc
            #     for (k,key) in enumerate(["gCCdisc_cc","gC8disc_cc"])
            #         println("         - $key")

            #         obs = corr[key][t]

            #         int = k==1 ? (obs .* TMRwbeta[t]) : (obs .* TMRw[t])
            #         amu = (alpha/pi)^exp_diag * sum(int) * 1e10

            #         HVP[key] = amu
            #     end
            # end

            if WRITE
                println("      - Writing BDIO & JDL2...")

                DERstr = STD_DERIV ? "_std" : "" 
                BLstr  = BLIND ? "Blind" : ""
                RESstr = RESC ? "_resc" : ""
                pBDIO  = create_path(path_bdio,["HVP&FVC",wind,ens.id,"$(ens.id)_$(BLstr)HVP$(diag)_set$(impr_set)$(RESstr)$(DERstr)"],OVERWRITE=true)
                pjdl2  = create_path(path_bdio,["HVP&FVC",wind,ens.id,"$(ens.id)_$(BLstr)HVP$(diag)_info_set$(impr_set)$(RESstr)$(DERstr).jld2"],OVERWRITE=true)

                io = IOBuffer()
                write(io, "$(ens.id) HVP")
                fb = ALPHAdobs_create(pBDIO, io)

                extra = Dict{String, Any}("ens" => ens.id, "impr_set" => impr_set, "diag" => diag, "wind" => wind)
                ALPHAdobs_write(fb, HVP, extra=extra)

                ALPHAdobs_close(fb)

                save(pjdl2,"HVPinfo",info)
            end
        end # end impr_set loop
    end # end ens loop
end # end timer


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> 1D HVP SDsub computation (charm) [LO, NLOa, NLOb] <==========================##

# diag = "NLOa"  # LO  NLOa  NLOb  NLOa&b

# IMPR_SET = ["1"] # ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

# STD_DERIV = false
# RESC      = false

# BLIND = false

# WRITE = false
# # OVERWRITE = true  # Always true in here!! 


# path_bdio = path_bdio_dict["local"]


# @info("STARTING CHARM HVP COMPUTATION [diag. $diag; wind. subtracted SD]")
# STD_DERIV ? @info("STANDARD DERIVATIVE is being employed in the IMPROVEMENT") : nothing


# # LO => (alpha/pi)^2 || NLO => (alpha/pi)^3
# exp_diag = diag == "LO" ? 2 : 3

# @time begin
#     for ens in ensInfo

#         @info("Computing HVP for ensemble $(ens.id)")
#         ens.id ∈ ensNOdisc ? @info("  > NO DISC. DATA FOR $(ens.id)") : nothing
#         if ens.id ∈ ensNOcharm
#             @info("  > NO CHARM DATA FOR $(ens.id)\n    Skipping this ensemble")
#             continue
#         end

#         println("   - Reading TMR...")
#         TMRbeta = BDIOread_TMR(path_bdio,ens,diag,resc=RESC,beta=true)

#         T = HVPobs.Data.get_T(ens.id)
#         sym_points = Int64(T/2+1)
#         t = collect(1:sym_points)

#         if !RESC
#             println("   - Reading t0...")
#             # t0_sym = BDIOread_t0_SU3sym(path_bdio,ens)
#             t0_sym = t0sym(ens.beta)
#             a_ß = sqrtt0_ph / sqrt(t0_sym); tfm_ß = a_ß.*(collect(1:length(TMRbeta)).-1)
#         else
#             println("   - Reading fPi...")
#             a_ß = hbarc * fPiph(ens.beta) / fPi_ph; tfm_ß = a_ß.*(collect(1:length(TMRbeta)).-1)
#         end

#         TMRwbeta = TMRbeta .* Window("SD")(tfm_ß)

#         TMRb_SU3(Q::Float64)   = ((16/(Q/factor_ß)^2)^2 * π^2 * (massmu/factor_ß)^2) .* C4[diag].((massmu/factor_ß).*(t .- 1)) .* sin.((Q/factor_ß/4) .* (t .- 1)).^4
#         TMRsub_SU3(Q::Float64) = TMRwbeta .- (Window("SD")(0) .* TMRb_SU3(Q))
        

#         for impr_set in IMPR_SET

#             println("   - Starting set "*impr_set)

#             println("      - Reading HVP(!!)...")
#             HVP, info = BDIOread_HVPens(path_bdio,diag,"SDsub",ens.id,impr_set,info=true,resc=RESC,STD=STD_DERIV,BLIND=BLIND)

#             println("      - Reading corr...")
#             corr = BDIOread_corr(path_bdio,ens,impr_set,STD=STD_DERIV)

#             println("      - Computing HVPs...\n        For Q = $Qlist")

#             HVPQ = Dict{String, Array{uwreal}}()
#             # HVP = Dict{String, Array{uwreal}}()
#             HVPsyst = Dict{String, Array{Float64}}()


#             for key in ["gCCconn_ll_sim","gCCconn_lc_sim","gCCconn_ll_sim+","gCCconn_lc_sim+"]
#                 println("         - $key  &  $(key[1:8])SU3$(key[8:end])  &  ∆lc_b$(key[8:end])")

#                 obs  = corr[key][t]
#                 # ∆obs = obs .- 2 .* corr["g33_"*key[9:10]][t]

#                 # HVPQ[key] = []
#                 # HVPQ[key[1:8]*"SU3"*key[8:end]] = []
#                 # HVPQ["∆lc_b"*key[8:end]] = []
#                 # HVPQ["∆lc_b_"*key[5:6]*"_beta_"*key[13:end]] = []
#                 for Q in Qlist
#                     # int    = obs .* TMRsub(Q)
#                     intSU3 = obs .* TMRsub_SU3(Q)
#                     # intb   = ∆obs .* TMRb(Q)
#                     # intbSU3 = ∆obs .* TMRb_SU3(Q)

#                     # amu    = (alpha/pi)^exp_diag * sum(int) * 1e10
#                     amuSU3 = (alpha/pi)^exp_diag * sum(intSU3) * 1e10
#                     # amub   = (alpha/pi)^exp_diag * sum(intb) * 1e10
#                     # amubSU3 = (alpha/pi)^exp_diag * sum(intbSU3) * 1e10

#                     # push!(HVPQ[key], amu)
#                     push!(HVPQ[key[1:8]*"SU3"*key[8:end]], amuSU3)
#                     # push!(HVPQ["∆lc_b"*key[8:end]], amub)
#                     # push!(HVPQ["∆lc_b_"*key[5:6]*"_beta_"*key[13:end]], amubSU3)
#                 end
#             end

#             println("      - Reading kappaC target...")

#             mDs_ph_prime, mDs_beta, Ds_dict = BDIOread_mDs_kappaC(path_bdio,ens)
#             kappaC_tar = Ds_dict["kappaC"]
#             println("         - Interpolating charm contribution...")

#             kappaC = [kcd_in[ens.id]["kappaC_sim"],kcd_in[ens.id]["kappaC_sim_plus"]]

#             for key in ["gCCconn_SU3_ll","gCCconn_SU3_lc"] # ["gCCconn_ll","gCCconn_lc","gCCconn_SU3_ll","gCCconn_SU3_lc","∆lc_b_ll","∆lc_b_lc"]
#                 HVPQ[key] = []
#                 HVPsyst[key] = []
#                 for i=1:length(Qlist)
#                     obs  = [HVPQ["$(key)_sim"][i],HVPQ["$(key)_sim+"][i]]
                
#                     par, _  = lin_fit(kappaC,obs,lineprint=false)
#                     obs_tar = y_lin_fit(par,kappaC_tar)

#                     if kappaC[1] <= value(kappaC_tar) <= kappaC[2]
#                         push!(HVPQ[key], obs_tar)
#                         push!(HVPsyst[key], 0.0)
#                     else
#                         push!(HVPQ[key], obs_tar)

#                         DkappaC = abs(kappaC[1]-kappaC[2])
#                         DkappaC_tar = minimum([abs(value(kappaC_tar)-kappaC[1]),abs(value(kappaC_tar)-kappaC[2])])
#                         uwerr(obs_tar); push!(HVPsyst[key], 1.5 * (abs(DkappaC_tar) / DkappaC)^2 * err(obs_tar))
#                     end
#                 end
#             end

#             # if (ens.id ∉ ensNOdisc && ens.kappa_l != ens.kappa_s)
#             #     println("      - Computing ∆ls(amu)...")
#             #     for discr in ["ll","lc"]
#             #         key = "∆ls_amu_$discr"
#             #         println("         - $key")
#             #         obs = corr["g88_$discr"][t] .- corr["g33_$discr"][t]

#             #         int = obs .* TMRw
#             #         amu = (alpha/pi)^exp_diag * sum(int) * 1e10

#             #         HVP[key] = [amu]
#             #     end
#             #     println("      - Computing gCdisc...")
#             #     for (k,key) in enumerate(["gCCdisc_cc","gC8disc_cc"])
#             #         println("         - $key")

#             #         obs = corr[key][t]

#             #         int = k==1 ? (obs .* TMRwbeta) : (obs .* TMRw)
#             #         amu = (alpha/pi)^exp_diag * sum(int) * 1e10

#             #         HVP[key] = [amu]
#             #     end
#             # else
#             #     println("         - 'G33 = G88 & GCCdisc = 0' for $(ens.id)")
#             # end

#             for key in ["gCCconn_SU3_ll","gCCconn_SU3_lc"]
#                 HVP[]
#                 info["HVPsyst"][key] = HVPsyst[key]
#             end

#             if WRITE
#                 println("      - Writing BDIO & JDL2...")

#                 DERstr = STD_DERIV ? "_std" : "" 
#                 RESstr = RESC ? "_resc" : ""
#                 pBDIO  = create_path(path_bdio,["HVP&FVC","SDsub",ens.id,"$(ens.id)_HVP$(diag)_set$(impr_set)$(RESstr)$(DERstr)"],OVERWRITE=true)
#                 pjdl2  = create_path(path_bdio,["HVP&FVC","SDsub",ens.id,"$(ens.id)_HVP$(diag)_info_set$(impr_set)$(RESstr)$(DERstr).jld2"],OVERWRITE=true)

#                 io = IOBuffer()
#                 write(io, "$(ens.id) HVP")
#                 fb = ALPHAdobs_create(pBDIO, io)

#                 extra = Dict{String, Any}("ens" => ens.id, "impr_set" => impr_set, "diag" => diag, "wind" => "SDsub", "Qlist" => Qlist)
#                 ALPHAdobs_write(fb, HVPQ, extra=extra)
#                 (ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc) ? ALPHAdobs_write(fb, HVP, extra=extra) : nothing

#                 ALPHAdobs_close(fb)

#                 HVPinfo = Dict{String,Dict}("HVPsyst" => HVPsyst) 
#                 save(pjdl2,"HVPinfo",HVPinfo)
#             end
#         end # end impr_set loop
#     end # end ens loop
# end # end timer



## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> READING TEST <==========================##

diag     = "NLOa"  #  "LO"  "NLOa"  "NLOb"  "NLOa&b"  "NLOc"
wind     = "ID"  #  "NW"  "SD"  "ID"  "LD"  "ILD"
ensid    = "E250"
impr_set = "1"

STD   = false
VREF  = false
RESC  = false

BLIND = false

path_bdio = path_bdio_dict["local"]

# ensid = "D200"; HVP, info = BDIOread_HVPens(path_bdio,diag,wind,ensid,impr_set,info=true,resc=RESC,STD=STD,BLIND=BLIND); println("$(ensid) -> $(print_uwreal(HVP["g33_ll"]*10))")
HVP, info = BDIOread_HVPens(path_bdio,diag,wind,ensid,impr_set,info=true,resc=RESC,STD=STD,BLIND=BLIND)

@info("Ready")

##

print_uwreal(HVP["gCCconn_SU3_lc"])
details(HVP["gCCconn_SU3_lc"])