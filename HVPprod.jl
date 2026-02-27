# Import packages

using Revise

include("HVPtool/HVPtool.jl")
using .HVPtool
using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
# using OrderedCollections

using TimerOutputs
using Suppressor

# include uwreal constants

# include("HVPtool/uwConst.jl")

# Path definition

julia_script_directory = @__DIR__

path_hvp_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "HVPData"),
    "SSD"   => joinpath(julia_script_directory, "..", "ObsExternal","PortableSSD","HVPData"),
    "clust" => joinpath(julia_script_directory, "..", "ObsExternal", "mogon_mount", "HVPdata")
)

path_bdio_dict = Dict{String,String}(
    "local" => joinpath(julia_script_directory, "..", "ObsBDIO"),
    "SSD"   => joinpath(julia_script_directory, "..", "ObsExternal","PortableSSD","ObsBDIO"),
    "clust" => joinpath(julia_script_directory, "..", "ObsExternal", "mogon_mount", "ObsBDIO")
)

path_coef  = joinpath(julia_script_directory, "..", "KernelCoeff")

path_HVP   = joinpath(path_hvp_dict["local"], "HVP_data")
path_rw_   = joinpath(path_hvp_dict["local"], "reweight")
path_ms    = joinpath(path_hvp_dict["local"], "ms_t0_dat")

path_fvcPI    = joinpath(path_hvp_dict["local"], "FSE_HP", "inf", "JKMPI_Mvmd")  # _Mvmd
path_fvcK     = joinpath(path_hvp_dict["local"], "FSE_HP", "inf", "JKMK")
path_fvcPIref = joinpath(path_hvp_dict["local"], "FSE_HP", "ref", "JKMPI_Mvmd")
path_fvcKref  = joinpath(path_hvp_dict["local"], "FSE_HP", "ref", "JKMK")


# Total list of fully functional ensembles:
# "A653","A654","B450","C101","C102","D150","D200","D201","D251","D450","D451","D452","E250","E300","F300","H101","H102","H200","J303","J304","J306","J307","J500","J501","N101","N200","N202","N203","N300","N302","N451","N452","S400"

# Problematic ensembles: 
# "H105"

ensList = ["A653","A654","B450","C101","C102","D150","D200","D201","D251","D450","D451","D452","E250","E300","F300","H101","H102","H200","H650","J303","J304","J306","J307","J500","J501","N101","N200","N202","N203","N300","N302","N451","N452","S400"]

# We do not have charm or disconnected data for some of the ensembles
ensNOcharm = ["C102","D150","D201","D251","D451","F300","H200","H650","J304","J306","J307","J501","N451","N452"]
ensNOdisc  = ["F300","J306"]

# Ensemble check
# ensInfo, bad_ensInfo = ensCheck(EnsInfo.(ensList), ensNOcharm, ensNOdisc, path_HVP, path_rw_, path_ms, path_fvcPI, data_status=true)
ensInfo =  EnsInfo.(ensList)

# isempty(bad_ensInfo) ? @info("Enough information has been found for all ensembles") : @info("Not enough information has been found concerning ensables $(join(bad_ensInfo, ", "))\n")

@info("Ready")

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> Corr. extraction and treatment <==========================##

IMPR      = true
RENORM    = true
STD_DERIV = false

IMPR_SET  = ["1","2"]  # ["1"] ["2"] ["1","2"] ["1old","2"] ["1","1old","2"]

OVERWRITE = false


path_bdio_w =  path_bdio_dict["local"]

@time begin
    for ens in ensInfo

        @info("Reading corr ensemble: $(ens.id)")
        ens.id ∈ ensNOcharm ? @info("  > NO CHARM DATA FOR $(ens.id)") : nothing
        ens.id ∈ ensNOdisc  ? @info("  > NO DISC. DATA FOR $(ens.id)") : nothing

        # choose defleated rw if possible
        path_rw = !isempty(filter(x->occursin(ens.id, x), readdir(joinpath(path_rw_,"reweight_deflated"), join=true))) ? joinpath(path_rw_,"reweight_deflated") : path_rw_

        for impr_set in IMPR_SET
            println("   - Impr Set: ", impr_set)
            println("      - G33 ll and lc correlator...")

            if ens.id ∉ ["C101","E300","J303"]
                g33_ll, g33_lc = corr33(path_HVP, ens, path_rw=path_rw, impr=IMPR, impr_set=impr_set, cons=true, frw_bcwd=true, std=STD_DERIV)
            else
                g33_ll = corr33(path_HVP, ens, path_rw=path_rw, impr=IMPR, impr_set=impr_set, cons=false, frw_bcwd=true, std=STD_DERIV)
                println("        (No LMA for 'lc' -> using ratio)")
                g33stoc_ll, g33stoc_lc = corr33(path_HVP, ens, path_rw=path_rw, impr=IMPR, impr_set=impr_set, cons=true, frw_bcwd=true, std=STD_DERIV, lma=false)
                g33_lc_obs = g33_ll.obs .* (g33stoc_lc.obs./g33stoc_ll.obs)
                
                g33_lc = Corr(g33_lc_obs, ens.id, "G33lc")
            end

            println("      - G88 connected ll and lc correlator...")        
            g88_ll_conn, g88_lc_conn = corr88_conn(path_HVP, ens, g33_ll, g33_lc=g33_lc, path_rw=path_rw, impr=IMPR, impr_set=impr_set, cons=true, frw_bcwd=true, std=STD_DERIV)
            
            println("      - G08 connected ll and lc correlator...")
            g08_ll_conn, g08_lc_conn = corr08_conn(g33_ll, g88_ll_conn, g33_lc=g33_lc, g88_lc=g88_lc_conn)

            if ens.id ∉ ensNOcharm
                println("      - Gcc connected ll and lc correlator...")
                gcc_ll_conn, gcc_lc_conn = corrC_conn(path_HVP, ens, path_rw=path_rw, impr=IMPR, impr_set=impr_set, cons=true, frw_bcwd=true, std=STD_DERIV, plus=false)
                gcc_ll_conn_p, gcc_lc_conn_p = corrC_conn(path_HVP, ens, path_rw=path_rw, impr=IMPR, impr_set=impr_set, cons=true, frw_bcwd=true, std=STD_DERIV, plus=true)

            end

            if ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc
                println("      - G88, Gcc, G08, Gc8 disconnected...")
                g88_ll_disc, g88_lc_disc, g88_cc_disc = corrDisconnected(path_HVP, ens, "88", path_rw=path_rw, impr=IMPR, impr_set=impr_set, std=STD_DERIV)
                g08_ll_disc, g08_lc_disc, g08_cc_disc = corrDisconnected(path_HVP, ens, "08", path_rw=path_rw, impr=IMPR, impr_set=impr_set, std=STD_DERIV)
                g80_ll_disc, g80_lc_disc, g80_cc_disc = corrDisconnected(path_HVP, ens, "80", path_rw=path_rw, impr=IMPR, impr_set=impr_set, std=STD_DERIV)

                gcc_ll_disc, gcc_lc_disc, gcc_cc_disc = corrDisconnected(path_HVP, ens, "cc", path_rw=path_rw, impr=false, std=STD_DERIV)
                gc8_ll_disc, gc8_lc_disc, gc8_cc_disc = corrDisconnected(path_HVP, ens, "c8", path_rw=path_rw, impr=false, std=STD_DERIV)
                g8c_ll_disc, g8c_lc_disc, g8c_cc_disc = corrDisconnected(path_HVP, ens, "8c", path_rw=path_rw, impr=false, std=STD_DERIV)
            end

            if RENORM
                println("      - Renormalization...")

                Z3 = get_Z3(ens, impr_set=impr_set)
                renormalize!(g33_ll, Z3^2)
                renormalize!(g33_lc, Z3)
                
                Z8 = get_Z8(ens, impr_set=impr_set)
                renormalize!(g88_ll_conn, Z8^2)
                renormalize!(g88_lc_conn, Z8)
                if ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc
                    renormalize!(g88_ll_disc, Z8^2)
                    renormalize!(g88_lc_disc, Z8)
                end
    
                if ens.id ∉ ensNOcharm
                    Zcc = Zvc_l[ens.id]
                    renormalize!(gcc_ll_conn, Zcc*Zcc)
                    renormalize!(gcc_lc_conn, Zcc)
                    renormalize!(gcc_ll_conn_p, Zcc*Zcc)
                    renormalize!(gcc_lc_conn_p, Zcc)
                end
                if ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc && ens.id ∉ ensNOcharm
                    Zcc = Zvc_l[ens.id]
                    renormalize!(gcc_ll_disc, Zcc*Zcc)
                    renormalize!(gcc_lc_disc, Zcc)
                end
    
                Z08 = get_Z08(ens, impr_set=impr_set)
                renormalize!(g08_ll_conn, Z8*Z08)
                renormalize!(g08_lc_conn, Z08)
                if ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc
                    renormalize!(g08_ll_disc, Z8*Z08)
                    renormalize!(g08_lc_disc, Z08)
                    renormalize!(g80_ll_disc, Z8*Z08)
                    renormalize!(g80_lc_disc, Z8)
                end
            end
            
            println("      - Writing BDIO...")

            data_corr = Dict{String, Array{uwreal}}()

            data_corr["g33_ll"] = g33_ll.obs
            data_corr["g33_lc"] = g33_lc.obs
            
            data_corr["g88conn_ll"] = g88_ll_conn.obs
            data_corr["g88conn_lc"] = g88_lc_conn.obs

            if ens.id ∉ ensNOcharm
                data_corr["gCCconn_ll_sim"] = gcc_ll_conn.obs
                data_corr["gCCconn_lc_sim"] = gcc_lc_conn.obs
                data_corr["gCCconn_ll_sim+"] = gcc_ll_conn_p.obs
                data_corr["gCCconn_lc_sim+"] = gcc_lc_conn_p.obs
            end

            if ens.kappa_l != ens.kappa_s && ens.id ∉ ensNOdisc
                data_corr["g08conn_ll"] = g08_ll_conn.obs
                data_corr["g08conn_lc"] = g08_lc_conn.obs

                data_corr["g88disc_ll"] = g88_ll_disc.obs
                data_corr["g88disc_lc"] = g88_lc_disc.obs

                data_corr["gCCdisc_ll"] = gcc_ll_disc.obs
                data_corr["gCCdisc_lc"] = gcc_lc_disc.obs
                data_corr["gCCdisc_cc"] = gcc_cc_disc.obs

                data_corr["g08disc_ll"] = g08_ll_disc.obs
                data_corr["g08disc_lc"] = g08_lc_disc.obs

                data_corr["g80disc_ll"] = g80_ll_disc.obs
                data_corr["g80disc_lc"] = g80_lc_disc.obs
    
                data_corr["gC8disc_ll"] = gc8_ll_disc.obs
                data_corr["gC8disc_lc"] = gc8_lc_disc.obs
                data_corr["gC8disc_cc"] = gc8_cc_disc.obs

                data_corr["g8Cdisc_ll"] = g8c_ll_disc.obs
                data_corr["g8Cdisc_lc"] = g8c_lc_disc.obs
                data_corr["g8Cdisc_cc"] = g8c_cc_disc.obs
            end

            DERstr = STD_DERIV ? "_std" : "" 
            pBDIO  = create_path(path_bdio_w,["Corr&Kernel&t0",ens.id,"$(ens.id)_corr_set$(impr_set)$(DERstr)"],OVERWRITE=OVERWRITE)
            
            io = IOBuffer()
            write(io, "$(ens.id)  HVP correlators, improvement set $(impr_set)")
            fb = ALPHAdobs_create(pBDIO, io)

            extra = Dict{String, Any}("ens" => ens.id, "impr_set" => impr_set, "std" => STD_DERIV)
            ALPHAdobs_write(fb, data_corr, extra=extra)
            
            ALPHAdobs_close(fb)
        end # end impr_set loop
    end # end ens loop
end # end timer

##==========================> FVC extraction <==========================##

OVERWRITE = false  # Allows to erase data and overwrite it with new data, use carefully !!


path_bdio_w =  path_bdio_dict["local"]

@time begin
    for ens in ensInfo
        @info("Computing Hansen-Patella FVC for ensemble: $(ens.id)")

        println("   - Computing FVC...")

        fvc_rawPi = get_fvc(path_fvcPI, ens.id)
        rowsPi, _ = size(fvc_rawPi)
        FVCPi = [vcat(sum(fvc_rawPi[1:i, :], dims=1)...) for i in 1:rowsPi]

        fvc_rawK = get_fvc(path_fvcK, ens.id)
        rowsK, _ = size(fvc_rawK)
        FVCK = [vcat(sum(fvc_rawK[1:i, :], dims=1)...) for i in 1:rowsK]

        println("   - Computing FVC (at Lref)...")

        fvc_rawPiref = get_fvc(path_fvcPIref, ens.id)
        rowsPiref, _ = size(fvc_rawPiref)
        FVCPi_ref = [vcat(sum(fvc_rawPiref[1:i, :], dims=1)...) for i in 1:rowsPiref]

        fvc_rawKref = get_fvc(path_fvcKref, ens.id)
        rowsKref, _ = size(fvc_rawKref)
        FVCK_ref = [vcat(sum(fvc_rawKref[1:i, :], dims=1)...) for i in 1:rowsKref]

        println("   - Writing BDIO...")

        data_FVC = Dict{String, Array{uwreal}}(
            "FVCPi1" => FVCPi[1],
            "FVCPi2" => FVCPi[2],
            "FVCPi3" => FVCPi[3],
            "FVCPi4" => FVCPi[4],
            "FVCPi5" => FVCPi[5],
            "FVCPi6" => FVCPi[6],
            "FVCK1"  => FVCK[1],
            "FVCK2"  => FVCK[2],
            "FVCK3"  => FVCK[3],
            "FVCK4"  => FVCK[4],
            "FVCK5"  => FVCK[5],
            "FVCK6"  => FVCK[6],
        )

        data_FVCref = Dict{String, Array{uwreal}}(
            "FVCPi1" => FVCPi[1] - FVCPi_ref[1],
            "FVCPi2" => FVCPi[2] - FVCPi_ref[2],
            "FVCPi3" => FVCPi[3] - FVCPi_ref[3],
            "FVCPi4" => FVCPi[4] - FVCPi_ref[4],
            "FVCPi5" => FVCPi[5] - FVCPi_ref[5],
            "FVCPi6" => FVCPi[6] - FVCPi_ref[6],
            "FVCK1"  => FVCK[1] - FVCK_ref[1],
            "FVCK2"  => FVCK[2] - FVCK_ref[2],
            "FVCK3"  => FVCK[3] - FVCK_ref[3],
            "FVCK4"  => FVCK[4] - FVCK_ref[4],
            "FVCK5"  => FVCK[5] - FVCK_ref[5],
            "FVCK6"  => FVCK[6] - FVCK_ref[6],
        )

        pBDIO  = create_path(path_bdio_w,["Corr&Kernel&t0",ens.id,"$(ens.id)_FVC"],OVERWRITE=OVERWRITE)

        io = IOBuffer()
        write(io, "$(ens.id) FVC")
        fb = ALPHAdobs_create(pBDIO, io)
        # fb = ALPHAdobs_create(joinpath("/Users/cesc/Desktop","$(ens.id)_FVC"), io)

        extra = Dict{String, Any}("ens" => ens.id, "Vref" => false)
        ALPHAdobs_write(fb, data_FVC, extra=extra)
        
        ALPHAdobs_close(fb)

        pBDIOref  = create_path(path_bdio_w,["Corr&Kernel&t0",ens.id,"$(ens.id)_FVC_Vref"],OVERWRITE=OVERWRITE)

        io = IOBuffer()
        write(io, "$(ens.id) FVC (Vref)")
        fb = ALPHAdobs_create(pBDIOref, io)
        # fb = ALPHAdobs_create(joinpath("/Users/cesc/Desktop","$(ens.id)_FVC"), io)

        extra = Dict{String, Any}("ens" => ens.id, "Vref" => true)
        ALPHAdobs_write(fb, data_FVCref, extra=extra)
        
        ALPHAdobs_close(fb)
    end # end ens loop
end # end timer

##==========================> t0 computation <==========================##

OVERWRITE = false

PLOT = true


path_bdio_w = path_bdio_dict["local"]

@time begin
    for ens in ensInfo
        @info("Computing t0 for ensemble: $(ens.id)")

        # choose defleated rw if possible       
        path_rw = !isempty(filter(x->occursin(ens.id, x), readdir(joinpath(path_rw_,"reweight_deflated"), join=true))) ? joinpath(path_rw_,"reweight_deflated") : path_rw_

        println("   - Computing t0...")
        t0 = uwreal[]
        @suppress begin
            try
                t0 = get_t0(path_ms, ens, path_rw = path_rw, pl = PLOT, wpm = wpmm)
            catch
                println("     (loading uwerr with label = $(ens.id))")
                _ = BDIOread_corr(path_bdio_w,ens.id,"1",STD=false)
                t0 = get_t0(path_ms, ens, path_rw = path_rw, pl = PLOT)
            end
        end
        
        data_t0 = Dict{String,Array{uwreal}}("t0" => [t0])

        println("   - Writing BDIO...")

        pBDIO  = create_path(path_bdio_w,["Corr&Kernel&t0",ens.id,"$(ens.id)_t0"],OVERWRITE=OVERWRITE)

        io = IOBuffer()
        write(io, "$(ens.id) t0")
        fb = ALPHAdobs_create(pBDIO, io)

        extra = Dict{String, Any}("ens" => ens.id)
        ALPHAdobs_write(fb, data_t0, extra=extra)
        
        ALPHAdobs_close(fb)
    end # end ens loop
end # end timer

##==========================> LO & NLO TMR computation <==========================##

scale = ""  # t0  t0su3  fPi  fPiph

COMPTMRc  = false
OVERWRITE = false

path_bdio_w = path_bdio_dict["local"]

@time begin
    for ens in ensInfo
        @info("Computing TMRs for ensemble: $(ens.id)\n With scale: $scale")
        
        if scale == "t0"
            t0 = BDIOread_t0(path_bdio_w,ens.id)
            factor = hbarc * sqrt(t0)/sqrtt0_ph
        elseif scale  == "t0su3"
            # t0 = BDIOread_t0_SU3sym(path_bdio_w,ens.beta)
            t0 = t0sym(ens.beta)
            factor = hbarc * sqrt(t0)/sqrtt0_ph
        elseif scale == "fPi"
            fPi = BDIOread_fPS(path_bdio_w,ens.id)["fPi"]
            factor = fPi_ph/fPi
        elseif scale == "fPiph"
            fPi = fPiph(ens.beta)
            factor = fPi_ph/fPi
        else
            error("Scale = $scale cannot be accepted; please choose between 't0',  't0su3', 'fPi' or 'fPiph'")
        end

        t_hat = (massmu/factor) .* collect(0:HVPobs.Data.get_T(ens.id))
        t_hat = t_hat[t_hat .< 4.]  # we compute up to T or ≈7.47fm

        println("   - TMR for diagram 'LO'...")
        TMR = factor^2 .* Tildef2(t_hat,path_coef)
        println("   - TMR for diagram 'NLOa'...")
        TMRa = factor^2 .* Tildef4a(t_hat,path_coef)
        println("   - TMR for diagram 'NLOb'...")
        TMRb = factor^2 .* Tildef4b(t_hat,path_coef)
        # COMPTMRc = ens.id ∉ ensNOcharm
        if  COMPTMRc
            println("   - TMR for diagram 'NLOc'...")
            TMRc = factor^4 .* Tildef4c(t_hat,path_coef)
        end


        data_TMRLO = Dict{String, Array{uwreal}}(
            "TMR" => TMR
        )

        data_TMRNLOab = Dict{String, Array{uwreal}}(
            "TMRa" => TMRa,
            "TMRb" => TMRb,
        )

        if  COMPTMRc
            data_TMRNLOc = Dict{String, Array{uwreal}}(
                "TMRc" => TMRc,
            )
        end

        println("   - Writing BDIO...")

        pBDIOLO  = create_path(path_bdio_w,["Corr&Kernel&t0",ens.id,"$(ens.id)_TMR$(scale)_LO"]     ,OVERWRITE=OVERWRITE)
        pBDIONLOab = create_path(path_bdio_w,["Corr&Kernel&t0",ens.id,"$(ens.id)_TMR$(scale)_NLOab"],OVERWRITE=OVERWRITE)
        if  COMPTMRc
            pBDIONLOc  = create_path(path_bdio_w,["Corr&Kernel&t0",ens.id,"$(ens.id)_TMR$(scale)_NLOc"] ,OVERWRITE=OVERWRITE)
        end

        io = IOBuffer()
        write(io, "$(ens.id) TMR LO Kernel")
        fb = ALPHAdobs_create(pBDIOLO, io)

        extra = Dict{String, Any}("ens" => ens.id, "scale" => scale)
        ALPHAdobs_write(fb, data_TMRLO, extra=extra)
        
        ALPHAdobs_close(fb)


        io = IOBuffer()
        write(io, "$(ens.id) TMR NLOa & NLOb Kernels")
        fb = ALPHAdobs_create(pBDIONLOab, io)

        extra = Dict{String, Any}("ens" => ens.id, "scale" => scale)
        ALPHAdobs_write(fb, data_TMRNLOab, extra=extra)
        
        ALPHAdobs_close(fb)

        if  COMPTMRc
            io = IOBuffer()
            write(io, "$(ens.id) TMR NLOc Kernels")
            fb = ALPHAdobs_create(pBDIONLOc, io)

            extra = Dict{String, Any}("ens" => ens.id, "scale" => scale)
            ALPHAdobs_write(fb, data_TMRNLOc, extra=extra)
            
            ALPHAdobs_close(fb)
        end
    end # end ens loop
end # end timer

## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> NLO TMR blinding <==========================##

# using Distributions

scale = ""  # t0  fPi

OVERWRITE = false  # Allows to erase data and overwrite it with new data, use carefully !!


path_bdio = path_bdio_dict["local"]

# BLIND = value(BDIOread_TMR(path_bdio,"A653","NLOb",BLIND=true)[2] / BDIOread_TMR(path_bdio,"A653","NLOb",BLIND=false)[2])  #  rand(Uniform(0.5, 2))

@time begin
    for ens in ensInfo
        @info("Computing blinded TMRs for ensemble: $(ens.id)")

        println("   - Reading t0...")

        if scale == "t0"
            t0 = BDIOread_t0(path_bdio,ens.id)
            factor = hbarc * sqrt(t0)/sqrtt0_ph
        elseif scale == "fPi"
            fPi = BDIOread_fPS(path_bdio,ens.id)["fPi"]
            factor = fPi_ph/fPi
        else
            error("Scale = $scale cannot be accepted; please choose between 't0',  't0su3', 'fPi' or 'fPiph'")
        end

        t_hat = (massmu/factor) .* collect(0:HVPobs.Data.get_T(ens.id))
        t_hat = t_hat[t_hat .< 4.]  # we compute up to T or ≈7.47fm

        println("   - TMR for diagram 'NLOa'...")
        TMRa = factor^2 .* BLIND .* Tildef4a(t_hat,path_coef)
        println("   - TMR for diagram 'NLOb'...")
        TMRb = factor^2 .* BLIND .* Tildef4b(t_hat,path_coef)

        data_TMRNLOab = Dict{String, Array{uwreal}}(
            "TMRa" => TMRa,
            "TMRb" => TMRb,
        )
        
        println("   - Writing BDIO...")

        pBDIONLOab = create_path(path_bdio,["Corr&Kernel&t0",ens.id,"$(ens.id)_BlindTMR$(scale)_NLOab"],OVERWRITE=OVERWRITE)

        io = IOBuffer()
        write(io, "$(ens.id) TMR NLOa & NLOb Kernels")
        fb = ALPHAdobs_create(pBDIONLOab, io)

        extra = Dict{String, Any}("ens" => ens.id, "scale" => scale)
        ALPHAdobs_write(fb, data_TMRNLOab, extra=extra)
        
        ALPHAdobs_close(fb)
    end # end ens loop
end # end timer

# println("- Writing TXT...")

# open(joinpath(path_bdio,"..","BlindVal.txt"), "w") do io
#     write(io, string(BLIND))  # Convert to string and write to file
# end

# @info("Blinding complete!")


## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##
## <<------------------------------------------------------------------------------------------------------------------------>> ##

##==========================> READING TEST <==========================##

diag     = ""
ensid    = ""
impr_set = ""

STD_DERIV = false
VREF      = false
RESC      = false

BLIND = false

path_bdio_r = path_bdio_dict["local"]

t0  = BDIOread_t0(path_bdio_r,ensid)
t0s = BDIOread_t0_SU3sym(path_bdio_r,ensid)

fPi = BDIOread_fPS(path_bdio_r,ensid)["fPi"]

TMR = BDIOread_TMR(path_bdio_r,ensid,diag,resc=RESC,beta=false,BLIND=BLIND)

if ensid ∉ ensNOcharm
    TMRb = BDIOread_TMR(path_bdio_r,ensid,diag,resc=RESC,beta=true,BLIND=BLIND)
end

corr = BDIOread_corr(path_bdio_r,ensid,impr_set,STD=STD_DERIV)

fvc_hp  = BDIOread_FVCcorr(path_bdio_r,ensid,Vref=VREF)

t_mll, fvc_mll  = TXTread_FVCcorr_GS(joinpath(julia_script_directory, "..", "HVPData", "FSE_MLL"),ensid,Vref=VREF)
