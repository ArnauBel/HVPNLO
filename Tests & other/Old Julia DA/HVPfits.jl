# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

using ProgressBars
using Suppressor

# Include Isovector Model and Model Average tools

include("tools/const.jl")

# BDIO path definition

julia_script_directory = @__DIR__

path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "S400", "N203", "N200", "D200", "N302"]
ensInfo = EnsInfo.(ensList)

DictComptoKey = Dict{String,Vector{String}}(
    "33"      => ["g33_ll","g33_lc"],
    "88"      => ["g88_ll","g88_lc"],
    "cc conn" => ["gcc_ll_conn","gcc_lc_conn"],
    "cc disc" => ["gcc_cc_disc"],
    "c8 disc" => ["gc8_cc_disc"],

    "3333"    => ["g3333_ll","g3333_lc"],
    "8888"    => ["g8888_ll","g8888_lc"],
    "CCCC"    => ["gCCCC_ll","gCCCC_lc"],
    "3388"    => ["g3388_ll","g3388_lc"],
    "33CC"    => ["g33CC_ll","g33CC_lc"],
    "88CC"    => ["g88CC_ll","g88CC_lc"]
)

ensSU3 = []
for ens in ensInfo
    ens.kappa_l == ens.kappa_s ? push!(ensSU3,ens) : continue
end

println("Ens in SU(3) flavour sym. point: $(getfield.(ensSU3,:id))")

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> Fits [a & b] <==========================##

COMP = "33" # 33, 88, cc conn, cc disc, c8 disc


mdof = COMP in ["cc disc","c8 disc"] ? 3 : 4                                              # impose a minimum to the dof for all fits 
nens = COMP in ["cc disc","c8 disc"] ? length(ensInfo)-length(ensSU3) : length(ensInfo)   # number of data ensembles
type_basemodel = "phi4"      #  phi4  simple
type_DA        = "All-(a4)"  # All  All-(a4)  All-(a4,a2phi4)  ||  Only(phi2inv,logphi2)  Only(a3cutoff,phi2inv,logphi2)
include("isovModel.jl")

baseDAdirect = joinpath(path_bdio,"DA","base[$type_basemodel]")
!ispath(baseDAdirect) ? mkdir(baseDAdirect) : nothing

DAdirect = joinpath(baseDAdirect,type_DA)
!ispath(DAdirect) ? mkdir(DAdirect) : nothing


mykeys = DictComptoKey[COMP]

phi2_ph = 8*(t0_ph*(mpi_ph*1e-3)/hbarc)^2
phi4_ph = 8*t0_ph^2*(((mK_ph*1e-3)/hbarc)^2 + 0.5*((mpi_ph*1e-3)/hbarc)^2)

@info(" Fitting for diagr. 'a' & 'b' and component $COMP\n         -> Basemodel: $type_basemodel & DA: $type_DA")

xdata = []
ydata = Dict{String, Dict}("1" => Dict{String,  Vector{Vector{uwreal}}}(), "2" => Dict{String, Vector{Vector{uwreal}}}())

i = 0
for ens in ensInfo
    println("- Reading data ensemble: $(ens.id)")
    if (COMP ∉ ["cc disc","c8 disc"]) || (ens.kappa_l != ens.kappa_s)
        i += 1

        println("   - Reading t0...")

        t0 = uwreal(0.0)
        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_t0"),"r")
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            t0 = ALPHAdobs_read_next(fb, size=sz, keys=ks)["t0"][1]
        end

        println("   - Reading HVP...")

        HVP = Dict{String, Dict}("1" => Dict{String, Vector{uwreal}}(), "2" => Dict{String, Vector{uwreal}}())

        for impr_set in ["1","2"]
            fb = BDIO_open(joinpath(path_bdio,"HVP&FVC",ens.id,"HVP","$(ens.id)_HVP_set$impr_set"),"r")
            val = Dict{String, Dict{String, Vector{uwreal}}}()
            while ALPHAdobs_next_p(fb)
                d = ALPHAdobs_read_parameters(fb)
                sz = tuple(d["size"]...)
                ks = collect(d["keys"])
                val["HVP"] = ALPHAdobs_read_next(fb, size=sz, keys=ks)
            end
            BDIO_close!(fb)
            info = load(joinpath(path_bdio,"HVP&FVC",ens.id,"HVP","$(ens.id)_HVPinfo_set$(impr_set).jld2"), "HVPinfo")
            HVPdict = merge(val,info)
            # apply syst.!
            if COMP != "cc conn"
                uwreal_syst = Dict{String, Vector{uwreal}}(); [uwreal_syst[key] = [uwreal([0.0, syst], "syst from BM/cut-off") for syst in HVPdict["HVPsyst"][key]] for key in mykeys]
            else
                uwreal_syst = Dict{String, Vector{uwreal}}(); uwreal_syst["gcc_ll_conn"] = uwreal_syst["gcc_lc_conn"] = [uwreal([0.0,0.0], "syst from BM/cut-off"),uwreal([0.0,0.0], "syst from BM/cut-off")]
            end
            [HVP[impr_set][key] = HVPdict["HVP"][key] .+ uwreal_syst[key] for key in mykeys]
        end

        println("   - Reading FVC...")

        FVC = Dict{String, Any}()

        fb_ab = BDIO_open(joinpath(path_bdio,"HVP&FVC",ens.id,"FVC","$(ens.id)_FVC_ab"),"r")
        val_ab = Dict{String, Matrix{uwreal}}()
        while ALPHAdobs_next_p(fb_ab)
            d = ALPHAdobs_read_parameters(fb_ab)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            val_ab = ALPHAdobs_read_next(fb_ab, size=sz, keys=ks)
        end
        BDIO_close!(fb_ab)
        for FVCtype in ["FVCPi","FVCK"]
            FVC[FVCtype] = [val_ab[FVCtype][end,1],val_ab[FVCtype][end,2]]
        end

        if ens.kappa_l == ens.kappa_s
            if COMP in ["33","88"]; myFVC = 1.5.*FVC["FVCPi"]; else myFVC = [0.0,0.0]; end
        else
            if COMP == "33"
                myFVC = FVC["FVCPi"] .+ FVC["FVCK"]
            elseif COMP == "88"
                myFVC = (2/9).*FVC["FVCK"] # (2/3).*FVC["FVCK"]
            else
                myFVC = [0.0,0.0]
            end
        end

        println("   - Creating 'x' data points...")

        push!(xdata, [1 / (8*t0), 8*t0*m_ens[ens.id]["m_pi"]^2, 8*t0*(m_ens[ens.id]["m_K"]^2+0.5*m_ens[ens.id]["m_pi"]^2)])

        println("   - Creating 'y' data points...")

        for impr_set in ["1","2"]
            for key in mykeys
                i == 1 ? ydata[impr_set][key] = Vector{Vector{uwreal}}() : nothing
                push!(ydata[impr_set][key], HVP[impr_set][key] .+ myFVC)
            end
        end

    else

        println("   - The contributions from $COMP are 0 for this ensemble !!")

    end
end

@info(" Fitting points")

# xdata = hcat([[1 / (8*t0[ens.id]), 8*t0[ens.id]*m_ens[ens.id]["m_pi"]^2, 8*t0[ens.id]*(m_ens[ens.id]["m_K"]^2+0.5*m_ens[ens.id]["m_pi"]^2)] for ens in ensInfo]...)
xdata = hcat(xdata...)
xdata = [xdata[1,:] xdata[2,:] xdata[3,:]]

fit = Dict{String, Dict}(); par = Dict{String,Dict}()
for (j,diag) in enumerate(["a","b"])
    println("   - Starting diag. "*diag)
    fit_diag = Dict{String, Dict}("1" => Dict{String, Vector{FitRes}}(), "2" => Dict{String, Vector{FitRes}}())
    par_diag = Dict{String, Dict}("1" => Dict{String, Vector{Vector{uwreal}}}(), "2" => Dict{String, Vector{Vector{uwreal}}}())
    for impr_set in ["1","2"]
        println("      - Starting set "*impr_set)
        for key in mykeys
            println("         - Fitting for comp. $key ...")
            fit_diag[impr_set][key] = Vector{FitRes}(); par_diag[impr_set][key] = Vector{uwreal}()
            for i in ProgressBar(collect(1:length(f_tot_isov)))
                @suppress begin
                    myfit = fit_routine(f_tot_isov[i], value.(xdata), hcat(ydata[impr_set][key]...)[j,:], n_par_tot_isov[i], pval=true)
                    push!(fit_diag[impr_set][key], myfit)
                    push!(par_diag[impr_set][key], myfit.param)
                end
            end
        end
    end
    fit[diag] = fit_diag
    par[diag] = par_diag
end

println("- Printing BDIO & JDL2...")

joinpath(DAdirect,"Fit")

diagpath = joinpath(DAdirect,"Fit","a&b")
!ispath(diagpath) ? mkdir(diagpath) : nothing

pcomp = joinpath(diagpath,COMP)
!ispath(pcomp) ? mkdir(pcomp) : nothing

pFitRes = joinpath(pcomp,"FitRes.jld2")
save(pFitRes,"FitRes",fit)

io = IOBuffer()
write(io, "parameters")

fb = ALPHAdobs_create(joinpath(pcomp,"param"), io) 

for i in collect(1:length(par["a"]["1"][mykeys[1]]))
    for diag in ["a","b"]
        for impr_set in ["1","2"]
            parDict = Dict{String,Array{uwreal}}()
            for key in mykeys
                parDict["diag$(diag)_$(key)_set$(impr_set):[$i]"] = par[diag][impr_set][key][i]

            end
            extra = Dict{String, Any}("Diag" => diag, "Set" => impr_set)
            ALPHAdobs_write(fb, parDict, extra=extra)
        end
    end
end
ALPHAdobs_close(fb)

io = IOBuffer()
write(io, "xydata")

fb = ALPHAdobs_create(joinpath(pcomp,"xydata"), io)

xDict = Dict{String,Array{uwreal}}("xdata" => xdata)
ALPHAdobs_write(fb, xDict)
for impr_set in ["1","2"]
    yDict = Dict{String,Array{uwreal}}()
    for key in mykeys
        yDict["$(key)_set$(impr_set)"] = hcat(ydata[impr_set][key]...)
    end
    extra = Dict{String, Any}("Set" => impr_set)
    ALPHAdobs_write(fb, yDict, extra=extra)
end
ALPHAdobs_close(fb)

println("- Printing Model information...")

infoDict = Dict{String,Any}(
    "length" => length(f_tot_isov),
    "nens" => nens,
    "n_par_tot_isov" => n_par_tot_isov,
    "label_tot_isov" => label_tot_isov,
)

pinfo = joinpath(pcomp,"ModelInfo.jld2")
save(pinfo,"info",infoDict)


##==========================> Fits [c] <==========================##

COMP = "88CC" # 3333, 8888, CCCC, 3388, 33CC, 88CC

mdof = 4                 # impose a minimum to the dof for all fits 
nens = length(ensInfo)   # number of data ensembles
type_basemodel = "phi4"
type_DA        = "All-(a4)"
include("isovModel.jl")

baseDAdirect = joinpath(path_bdio,"DA","base[$type_basemodel]")
!ispath(baseDAdirect) ? mkdir(baseDAdirect) : nothing

DAdirect = joinpath(baseDAdirect,type_DA)
!ispath(DAdirect) ? mkdir(DAdirect) : nothing

mykeys = DictComptoKey[COMP]

phi2_ph = 8*(t0_ph*(mpi_ph*1e-3)/hbarc)^2
phi4_ph = 8*t0_ph^2*(((mK_ph*1e-3)/hbarc)^2 + 0.5*((mpi_ph*1e-3)/hbarc)^2)

@info(" Fitting for diagr. 'c' and component $COMP")

xdata = []
ydata = Dict{String, Dict}("1" => Dict{String,  Vector{uwreal}}(), "2" => Dict{String, Vector{uwreal}}())

for (i,ens) in enumerate(ensInfo)
    println("- Reading data ensemble: $(ens.id)")

    println("   - Reading t0...")

    t0 = uwreal(0.0)
    fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_t0"),"r")
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        t0 = ALPHAdobs_read_next(fb, size=sz, keys=ks)["t0"][1]
    end

    for impr_set in ["1","2"]
        println("   - Impr. set $impr_set")

        println("      - Reading HVP...")

        HVP = Dict{String, uwreal}()

        fb = BDIO_open(joinpath(path_bdio,"HVP&FVC",ens.id,"HVP","$(ens.id)_HVPC_set$(impr_set)"),"r")
        val = Dict{String, Dict{String, uwreal}}()
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            nobs = d["nobs"]
            dims = d["dimensions"]
            ks = collect(d["keys"])

            val["HVP"] =  ALPHAdobs_read_next(fb, keys=ks)
        end
        BDIO_close!(fb)
        info = load(joinpath(path_bdio,"HVP&FVC",ens.id,"HVP","$(ens.id)_HVPinfoC_set$(impr_set).jld2"), "HVPinfo")
        HVPdict = merge(val,info)
        # apply syst.!
        uwreal_syst = Dict{String, uwreal}()
        if COMP in ["3333","8888","3388"]
            [uwreal_syst[key] = uwreal([0.0, HVPdict["HVPsyst"][key]], "syst from BM/cut-off") for key in mykeys]
        else
            [uwreal_syst[key] = uwreal([0.0, 0.0] , "syst from BM/cut-off") for key in mykeys]
        end
        [HVP[key] = HVPdict["HVP"][key] + uwreal_syst[key] for key in mykeys]

        if COMP in ["3333","8888","3388"]
            println("      - Reading FVC...")

            FVC = Dict{String, Vector{uwreal}}()

            fb_c = BDIO_open(joinpath(path_bdio,"HVP&FVC",ens.id,"FVC","$(ens.id)_FVC_c_set$(impr_set)"),"r")
            while ALPHAdobs_next_p(fb_c)
                d = ALPHAdobs_read_parameters(fb_c)
                sz = tuple(d["size"]...)
                ks = collect(d["keys"])
                FVC = ALPHAdobs_read_next(fb_c, size=sz, keys=ks)
            end
            BDIO_close!(fb_c)
        end

        println("      - Creating 'y' data points...")

        for key in mykeys
            i == 1 ? ydata[impr_set][key] = Vector{uwreal}() : nothing
            if COMP in ["3333","8888","3388"]
                push!(ydata[impr_set][key], HVP[key] + FVC[key][end])
            else
                push!(ydata[impr_set][key], HVP[key])
            end
        end
    end

    println("   - Creating 'x' data points...")

    push!(xdata, [1 / (8*t0), 8*t0*m_ens[ens.id]["m_pi"]^2, 8*t0*(m_ens[ens.id]["m_K"]^2+0.5*m_ens[ens.id]["m_pi"]^2)])
end

@info(" Fitting points")

# xdata = hcat([[1 / (8*t0[ens.id]), 8*t0[ens.id]*m_ens[ens.id]["m_pi"]^2, 8*t0[ens.id]*(m_ens[ens.id]["m_K"]^2+0.5*m_ens[ens.id]["m_pi"]^2)] for ens in ensInfo]...)
xdata = hcat(xdata...)
xdata = [xdata[1,:] xdata[2,:] xdata[3,:]]

fit = Dict{String, Dict}("1" => Dict{String, Vector{FitRes}}(), "2" => Dict{String, Vector{FitRes}}())
par = Dict{String, Dict}("1" => Dict{String, Vector{Vector{uwreal}}}(), "2" => Dict{String, Vector{Vector{uwreal}}}())
for impr_set in ["1","2"]
    println("      - Starting set "*impr_set)
    for key in mykeys
        println("         - Fitting for comp. $key ...")
        fit[impr_set][key] = Vector{FitRes}(); par[impr_set][key] = Vector{uwreal}()
        for i in ProgressBar(collect(1:length(f_tot_isov)))
            @suppress begin
                myfit = fit_routine(f_tot_isov[i], value.(xdata), ydata[impr_set][key], n_par_tot_isov[i], pval=true)
                push!(fit[impr_set][key], myfit)
                push!(par[impr_set][key], myfit.param)
            end
        end
    end
end


println("- Printing BDIO & JDL2...")

joinpath(DAdirect,"Fit")

diagpath = joinpath(DAdirect,"Fit","c")
!ispath(diagpath) ? mkdir(diagpath) : nothing

pcomp = joinpath(diagpath,COMP)
!ispath(pcomp) ? mkdir(pcomp) : nothing

pFitRes = joinpath(pcomp,"FitRes.jld2")
save(pFitRes,"FitRes",fit)

io = IOBuffer()
write(io, "parameters")

fb = ALPHAdobs_create(joinpath(pcomp,"param"), io) 

for i in collect(1:length(par["1"][mykeys[1]]))
    for impr_set in ["1","2"]
        parDict = Dict{String,Array{uwreal}}()
        for key in mykeys
            parDict["diagc_$(key)_set$(impr_set):[$i]"] = par[impr_set][key][i]
        end
        extra = Dict{String, Any}("Diag" => "c", "Set" => impr_set)
        ALPHAdobs_write(fb, parDict, extra=extra)
    end
end
ALPHAdobs_close(fb)

io = IOBuffer()
write(io, "xydata")

fb = ALPHAdobs_create(joinpath(pcomp,"xydata"), io)

xDict = Dict{String,Array{uwreal}}("xdata" => xdata)
ALPHAdobs_write(fb, xDict)
for impr_set in ["1","2"]
    yDict = Dict{String,Array{uwreal}}()
    for key in mykeys
        yDict["$(key)_set$(impr_set)"] = hcat(ydata[impr_set][key]...)
    end
    extra = Dict{String, Any}("Set" => impr_set)
    ALPHAdobs_write(fb, yDict, extra=extra)
end
ALPHAdobs_close(fb)

println("- Printing Model information...")

infoDict = Dict{String,Any}(
    "length" => length(f_tot_isov),
    "nens" => nens,
    "n_par_tot_isov" => n_par_tot_isov,
    "label_tot_isov" => label_tot_isov,
)

pinfo = joinpath(pcomp,"ModelInfo.jld2")
save(pinfo,"info",infoDict)

##==========================> Fits [LO] <==========================##

COMP = "33" # 33, 88, cc conn, cc disc, c8 disc


mdof = COMP in ["cc disc","c8 disc"] ? 3 : 4                                              # impose a minimum to the dof for all fits 
nens = COMP in ["cc disc","c8 disc"] ? length(ensInfo)-length(ensSU3) : length(ensInfo)   # number of data ensembles
type_basemodel = "simple"                 #  phi4  simple
type_DA        = "Only(phi2inv,logphi2)"  # All  All-(a4)  All-(a4,a2phi4)  ||  Only(phi2inv,logphi2)  Only(a3cutoff,phi2inv,logphi2)
include("isovModel.jl")

baseDAdirect = joinpath(path_bdio,"DA","base[$type_basemodel]")
!ispath(baseDAdirect) ? mkdir(baseDAdirect) : nothing

DAdirect = joinpath(baseDAdirect,type_DA)
!ispath(DAdirect) ? mkdir(DAdirect) : nothing


mykeys = DictComptoKey[COMP]

phi2_ph = 8*(t0_ph*(mpi_ph*1e-3)/hbarc)^2
phi4_ph = 8*t0_ph^2*(((mK_ph*1e-3)/hbarc)^2 + 0.5*((mpi_ph*1e-3)/hbarc)^2)

@info(" Fitting for LO and component $COMP")

xdata = []
ydata = Dict{String, Dict}("1" => Dict{String,  Vector{uwreal}}(), "2" => Dict{String, Vector{uwreal}}())

i = 0
for ens in ensInfo
    println("- Reading data ensemble: $(ens.id)")
    if (COMP ∉ ["cc disc","c8 disc"]) || (ens.kappa_l != ens.kappa_s)
        i += 1

        println("   - Reading t0...")

        t0 = uwreal(0.0)
        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_t0"),"r")
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            t0 = ALPHAdobs_read_next(fb, size=sz, keys=ks)["t0"][1]
        end

        println("   - Reading HVP...")

        HVP = Dict{String, Dict}("1" => Dict{String, uwreal}(), "2" => Dict{String, uwreal}())

        for impr_set in ["1","2"]
            fb = BDIO_open(joinpath(path_bdio,"HVP&FVC",ens.id,"HVP","$(ens.id)_HVPLO_set$(impr_set)"),"r")
            val = Dict{String, Dict{String, uwreal}}()
            while ALPHAdobs_next_p(fb)
                d = ALPHAdobs_read_parameters(fb)
                nobs = d["nobs"]
                dims = d["dimensions"]
                ks = collect(d["keys"])
                val["HVP"] =  ALPHAdobs_read_next(fb, keys=ks)
            end
            BDIO_close!(fb)
            info = load(joinpath(path_bdio,"HVP&FVC",ens.id,"HVP","$(ens.id)_HVPLOinfo_set$(impr_set).jld2"), "HVPinfo")
            HVPdict = merge(val,info)
            # apply syst.!
            if COMP != "cc conn"
                uwreal_syst = Dict{String, uwreal}(); [uwreal_syst[key] = uwreal([0.0, HVPdict["HVPsyst"][key]], "syst from BM/cut-off") for key in mykeys]
            else
                uwreal_syst = Dict{String, uwreal}(); uwreal_syst["gcc_ll_conn"] = uwreal_syst["gcc_lc_conn"] = uwreal([0.0,0.0], "syst from BM/cut-off")
            end
            [HVP[impr_set][key] = HVPdict["HVP"][key] + uwreal_syst[key] for key in mykeys]
        end

        println("   - Reading FVC...")

        FVC = Dict{String, Any}()

        fb_LO = BDIO_open(joinpath(path_bdio,"HVP&FVC",ens.id,"FVC","$(ens.id)_FVC_LO"),"r")
        val_LO = Dict{String, uwreal}()
        while ALPHAdobs_next_p(fb_LO)
            d = ALPHAdobs_read_parameters(fb_LO)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            val_LO = ALPHAdobs_read_next(fb_LO, size=sz, keys=ks)
        end
        BDIO_close!(fb_LO)
        [FVC[FVCtype] = val_LO[FVCtype][end] for FVCtype in ["FVCPi","FVCK"]]

        if ens.kappa_l == ens.kappa_s
            if COMP in ["33","88"]; myFVC = 1.5*FVC["FVCPi"]; else myFVC = 0.0; end
        else
            if COMP == "33"
                myFVC = FVC["FVCPi"] + FVC["FVCK"]
            elseif COMP == "88"
                myFVC = (2/9)*FVC["FVCK"] # (2/3).*FVC["FVCK"]
            else
                myFVC = 0.0
            end
        end

        println("   - Creating 'x' data points...")

        push!(xdata, [1 / (8*t0), 8*t0*m_ens[ens.id]["m_pi"]^2, 8*t0*(m_ens[ens.id]["m_K"]^2+0.5*m_ens[ens.id]["m_pi"]^2)])

        println("   - Creating 'y' data points...")

        for impr_set in ["1","2"]
            for key in mykeys
                i == 1 ? ydata[impr_set][key] = Vector{uwreal}() : nothing
                push!(ydata[impr_set][key], HVP[impr_set][key] + myFVC)
            end
        end

    else

        println("   - The contributions from $COMP are 0 for this ensemble !!")

    end
end

@info(" Fitting points for LO")

xdata = hcat(xdata...)
xdata = [xdata[1,:] xdata[2,:] xdata[3,:]]

fit = Dict{String, Dict}(); par = Dict{String,Dict}()

fit = Dict{String, Dict}("1" => Dict{String, Vector{FitRes}}(), "2" => Dict{String, Vector{FitRes}}())
par = Dict{String, Dict}("1" => Dict{String, Vector{Vector{uwreal}}}(), "2" => Dict{String, Vector{Vector{uwreal}}}())
for impr_set in ["1","2"]
    println("      - Starting set "*impr_set)
    for key in mykeys
        println("         - Fitting for comp. $key ...")
        fit[impr_set][key] = Vector{FitRes}(); par[impr_set][key] = Vector{uwreal}()
        for i in ProgressBar(collect(1:length(f_tot_isov)))
            @suppress begin
                myfit = fit_routine(f_tot_isov[i], value.(xdata), ydata[impr_set][key], n_par_tot_isov[i], pval=true)
                push!(fit[impr_set][key], myfit)
                push!(par[impr_set][key], myfit.param)
            end
        end
    end
end

println("- Printing BDIO & JDL2...")

!ispath(joinpath(DAdirect,"Fit")) ? mkdir(joinpath(DAdirect,"Fit")) : nothing

diagpath = joinpath(DAdirect,"Fit","LO")
!ispath(diagpath) ? mkdir(diagpath) : nothing

pcomp = joinpath(diagpath,COMP)
!ispath(pcomp) ? mkdir(pcomp) : nothing

pFitRes = joinpath(pcomp,"FitRes.jld2")
save(pFitRes,"FitRes",fit)

io = IOBuffer()
write(io, "parameters")

fb = ALPHAdobs_create(joinpath(pcomp,"param"), io) 

for i in collect(1:length(par["1"][mykeys[1]]))
    for impr_set in ["1","2"]
        parDict = Dict{String,Array{uwreal}}()
        for key in mykeys
            parDict["diagLO_$(key)_set$(impr_set):[$i]"] = par[impr_set][key][i]
        end
        extra = Dict{String, Any}("Diag" => "LO", "Set" => impr_set)
        ALPHAdobs_write(fb, parDict, extra=extra)
    end
end
ALPHAdobs_close(fb)

io = IOBuffer()
write(io, "xydata")

fb = ALPHAdobs_create(joinpath(pcomp,"xydata"), io)

xDict = Dict{String,Array{uwreal}}("xdata" => xdata)
ALPHAdobs_write(fb, xDict)
for impr_set in ["1","2"]
    yDict = Dict{String,Array{uwreal}}()
    for key in mykeys
        yDict["$(key)_set$(impr_set)"] = ydata[impr_set][key]
    end
    extra = Dict{String, Any}("Set" => impr_set)
    ALPHAdobs_write(fb, yDict, extra=extra)
end
ALPHAdobs_close(fb)

println("- Printing Model information...")

infoDict = Dict{String,Any}(
    "length" => length(f_tot_isov),
    "nens" => nens,
    "n_par_tot_isov" => n_par_tot_isov,
    "label_tot_isov" => label_tot_isov,
)

pinfo = joinpath(pcomp,"ModelInfo.jld2")
save(pinfo,"info",infoDict)

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> 4-fit method [LO] <==========================##

COMP = "33" # 33, 88, cc conn, cc disc, c8 disc

nens = COMP in ["cc disc","c8 disc"] ? length(ensInfo)-length(ensSU3) : length(ensInfo)   # number of data ensembles
mpi_cut = "<360"  #  all  <360  <300

if mpi_cut != "all"
    cut = parse(Int64,mpi_cut[2:4])
end


baseDAdirect = joinpath(path_bdio,"DA","4-fit")
!ispath(baseDAdirect) ? mkdir(baseDAdirect) : nothing

DAdirect = joinpath(baseDAdirect,"mpi[$mpi_cut]")
!ispath(DAdirect) ? mkdir(DAdirect) : nothing


mykeys = DictComptoKey[COMP]

phi2_ph = 8*(t0_ph*(mpi_ph*1e-3)/hbarc)^2
phi4_ph = 8*t0_ph^2*(((mK_ph*1e-3)/hbarc)^2 + 0.5*((mpi_ph*1e-3)/hbarc)^2)

@info(" Fitting for LO and component $COMP")

xdata = []
ydata = Dict{String, Dict}("1" => Dict{String,  Vector{uwreal}}(), "2" => Dict{String, Vector{uwreal}}())

i = 0; j = 0
for ens in ensInfo
    println("- Reading data ensemble: $(ens.id)")
    if (COMP ∉ ["cc disc","c8 disc"]) || (ens.kappa_l != ens.kappa_s)

        println("   - Reading t0...")

        t0 = uwreal(0.0)
        fb = BDIO_open(joinpath(path_bdio,"Corr&Kernel&t0",ens.id,ens.id*"_t0"),"r")
        while ALPHAdobs_next_p(fb)
            d = ALPHAdobs_read_parameters(fb)
            sz = tuple(d["size"]...)
            ks = collect(d["keys"])
            t0 = ALPHAdobs_read_next(fb, size=sz, keys=ks)["t0"][1]
        end

        factor = hbarc * sqrt(t0)/t0_ph  
        mpi = m_ens[ens.id]["m_pi"]*factor * 1e3

        if mpi_cut == "all" || mpi < cut
            i += 1

            println("   - Reading HVP...")

            HVP = Dict{String, Dict}("1" => Dict{String, uwreal}(), "2" => Dict{String, uwreal}())

            for impr_set in ["1","2"]
                fb = BDIO_open(joinpath(path_bdio,"HVP&FVC",ens.id,"HVP","$(ens.id)_HVPLO_set$(impr_set)"),"r")
                val = Dict{String, Dict{String, uwreal}}()
                while ALPHAdobs_next_p(fb)
                    d = ALPHAdobs_read_parameters(fb)
                    nobs = d["nobs"]
                    dims = d["dimensions"]
                    ks = collect(d["keys"])
                    val["HVP"] =  ALPHAdobs_read_next(fb, keys=ks)
                end
                BDIO_close!(fb)
                info = load(joinpath(path_bdio,"HVP&FVC",ens.id,"HVP","$(ens.id)_HVPLOinfo_set$(impr_set).jld2"), "HVPinfo")
                HVPdict = merge(val,info)
                # apply syst.!
                if COMP != "cc conn"
                    uwreal_syst = Dict{String, uwreal}(); [uwreal_syst[key] = uwreal([0.0, HVPdict["HVPsyst"][key]], "syst from BM/cut-off") for key in mykeys]
                else
                    uwreal_syst = Dict{String, uwreal}(); uwreal_syst["gcc_ll_conn"] = uwreal_syst["gcc_lc_conn"] = uwreal([0.0,0.0], "syst from BM/cut-off")
                end
                [HVP[impr_set][key] = HVPdict["HVP"][key] + uwreal_syst[key] for key in mykeys]
            end

            println("   - Reading FVC...")

            FVC = Dict{String, Any}()

            fb_LO = BDIO_open(joinpath(path_bdio,"HVP&FVC",ens.id,"FVC","$(ens.id)_FVC_LO"),"r")
            val_LO = Dict{String, uwreal}()
            while ALPHAdobs_next_p(fb_LO)
                d = ALPHAdobs_read_parameters(fb_LO)
                sz = tuple(d["size"]...)
                ks = collect(d["keys"])
                val_LO = ALPHAdobs_read_next(fb_LO, size=sz, keys=ks)
            end
            BDIO_close!(fb_LO)
            [FVC[FVCtype] = val_LO[FVCtype][end] for FVCtype in ["FVCPi","FVCK"]]

            if ens.kappa_l == ens.kappa_s
                if COMP in ["33","88"]; myFVC = 1.5*FVC["FVCPi"]; else myFVC = 0.0; end
            else
                if COMP == "33"
                    myFVC = FVC["FVCPi"] + FVC["FVCK"]
                elseif COMP == "88"
                    myFVC = (2/9)*FVC["FVCK"] # (2/3).*FVC["FVCK"]
                else
                    myFVC = 0.0
                end
            end

            println("   - Creating 'x' data points...")

            push!(xdata, [1 / (8*t0), 8*t0*m_ens[ens.id]["m_pi"]^2, 8*t0*(m_ens[ens.id]["m_K"]^2+0.5*m_ens[ens.id]["m_pi"]^2)])

            println("   - Creating 'y' data points...")

            for impr_set in ["1","2"]
                for key in mykeys
                    i == 1 ? ydata[impr_set][key] = Vector{uwreal}() : nothing
                    push!(ydata[impr_set][key], HVP[impr_set][key] + myFVC)
                end
            end
        else
            println("   - The pion mass cut ($mpi_cut) has excluded this ensemble")
            j += 1
        end
    else
        println("   - The contributions from $COMP are 0 for this ensemble !!")
    end
end

nens -= j

@info(" Fitting points for LO (data ens = $nens)")

# isov_model1(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+  p[4] .* (log.(x[:,2]) .- log.(value.(phi2_ph)))
# isov_model2(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+  p[4] .* (x[:,2].^2 .- value.(phi2_ph).^2)
# isov_model3(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+  p[4] .* (x[:,2].^(-1) .- value.(phi2_ph).^(-1))
# isov_model4(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+  p[4] .* (x[:,2] .* log.(x[:,2]) .- value.(phi2_ph) .* log.(value.(phi2_ph)))

# f_tot_isov     = [isov_model1,isov_model2,isov_model3,isov_model4]
# label_tot_isov = ["logphi2","phi2sqr","phi2inv","phi2log"]
# n_par_tot_isov = [4,4,4,4]

isov_model1(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+ p[4] .* x[:,1].^(3/2) .+  p[5] .* (log.(x[:,2]) .- log.(value.(phi2_ph)))
isov_model2(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+ p[4] .* x[:,1].^(3/2) .+  p[5] .* (x[:,2].^2 .- value.(phi2_ph).^2)
isov_model3(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+ p[4] .* x[:,1].^(3/2) .+  p[5] .* (x[:,2].^(-1) .- value.(phi2_ph).^(-1))
isov_model4(x,p) = p[1] .+ p[2] .* x[:,1] .+ p[3] .* (x[:,2] .- value.(phi2_ph)) .+ p[4] .* x[:,1].^(3/2) .+  p[5] .* (x[:,2] .* log.(x[:,2]) .- value.(phi2_ph) .* log.(value.(phi2_ph)))

f_tot_isov     = [isov_model1,isov_model2,isov_model3,isov_model4]
label_tot_isov = ["logphi2","phi2sqr","phi2inv","phi2log"]
n_par_tot_isov = [5,5,5,5]

xdata = hcat(xdata...)
xdata = [xdata[1,:] xdata[2,:] xdata[3,:]]

fit = Dict{String, Dict}(); par = Dict{String,Dict}()

fit = Dict{String, Dict}("1" => Dict{String, Vector{FitRes}}(), "2" => Dict{String, Vector{FitRes}}())
par = Dict{String, Dict}("1" => Dict{String, Vector{Vector{uwreal}}}(), "2" => Dict{String, Vector{Vector{uwreal}}}())
for impr_set in ["1","2"]
    println("      - Starting set "*impr_set)
    for key in mykeys
        println("         - Fitting for comp. $key ...")
        fit[impr_set][key] = Vector{FitRes}(); par[impr_set][key] = Vector{uwreal}()
        for i in ProgressBar(collect(1:length(f_tot_isov)))
            @suppress begin
                myfit = fit_routine(f_tot_isov[i], value.(xdata), ydata[impr_set][key], n_par_tot_isov[i], pval=true)
                push!(fit[impr_set][key], myfit)
                push!(par[impr_set][key], myfit.param)
            end
        end
    end
end


println("- Printing BDIO & JDL2...")

!ispath(joinpath(DAdirect,"Fit")) ? mkdir(joinpath(DAdirect,"Fit")) : nothing

diagpath = joinpath(DAdirect,"Fit","LO")
!ispath(diagpath) ? mkdir(diagpath) : nothing

pcomp = joinpath(diagpath,COMP)
!ispath(pcomp) ? mkdir(pcomp) : nothing

pFitRes = joinpath(pcomp,"FitRes.jld2")
save(pFitRes,"FitRes",fit)

io = IOBuffer()
write(io, "parameters")

fb = ALPHAdobs_create(joinpath(pcomp,"param"), io) 

for i in collect(1:length(par["1"][mykeys[1]]))
    for impr_set in ["1","2"]
        parDict = Dict{String,Array{uwreal}}()
        for key in mykeys
            parDict["diagLO_$(key)_set$(impr_set):[$i]"] = par[impr_set][key][i]
        end
        extra = Dict{String, Any}("Diag" => "LO", "Set" => impr_set)
        ALPHAdobs_write(fb, parDict, extra=extra)
    end
end
ALPHAdobs_close(fb)

io = IOBuffer()
write(io, "xydata")

fb = ALPHAdobs_create(joinpath(pcomp,"xydata"), io)

xDict = Dict{String,Array{uwreal}}("xdata" => xdata)
ALPHAdobs_write(fb, xDict)
for impr_set in ["1","2"]
    yDict = Dict{String,Array{uwreal}}()
    for key in mykeys
        yDict["$(key)_set$(impr_set)"] = ydata[impr_set][key]
    end
    extra = Dict{String, Any}("Set" => impr_set)
    ALPHAdobs_write(fb, yDict, extra=extra)
end
ALPHAdobs_close(fb)


println("- Printing Model information...")

infoDict = Dict{String,Any}(
    "length" => length(f_tot_isov),
    "nens" => nens,
    "n_par_tot_isov" => n_par_tot_isov,
    "label_tot_isov" => label_tot_isov,
)

pinfo = joinpath(pcomp,"ModelInfo.jld2")
save(pinfo,"info",infoDict)

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> READING TEST <==========================##

type_basemodel = "4-fit"      #  phi4  simple  4-fit
type_DA        = "All-(a4)"   #  All  All-(a4)  All-(a4,a2phi4)  ||  Only(phi2inv,logphi2)  Only(a3cutoff,phi2inv,logphi2)
mass_cut       = "<360"        #  all  <360  <300

impr_set = "1"
comp = "33"
diag = "LO"   # "a$b"  "c"  "LO"
extract_data = "param"   # "FitRes"  "param"  "xydata"  "MA"  "info"

if diag == "LO"
    DIAG = ["LO"]
elseif diag == "a&b"
    DIAG = ["a","b"]
elseif diag == "c"
    DIAG = ["c"]
end

if type_basemodel == "4-fit"
    fittype_dict = joinpath(path_bdio,"DA","4-fit","mpi[$mass_cut]")
else
    fittype_dict = joinpath(path_bdio,"DA","base[$type_basemodel]",type_DA)
end

mykeys = DictComptoKey[comp]


if extract_data == "FitRes"
    res = load(joinpath(fittype_dict,"Fit",diag,comp,"FitRes.jld2"), "FitRes")
elseif extract_data == "param"
    modelinfo = load(joinpath(fittype_dict,"Fit",diag,comp,"ModelInfo.jld2"), "info")

    fb = BDIO_open(joinpath(fittype_dict,"Fit",diag,comp,"param"),"r")
    partial_res = Vector{Dict}()
    full_dict = Dict{String, Any}()
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        push!(partial_res,ALPHAdobs_read_next(fb, size=sz, keys=ks))
    end
    BDIO_close!(fb)
    for dict in partial_res
        merge!(full_dict, dict)
    end
    res = Dict{String, Any}()
    for diag in DIAG
        res[diag] = Dict{String, Dict}()
        for impr_set in ["1","2"]
            res[diag][impr_set] = Dict{String, Vector{Vector{uwreal}}}()
            for key in mykeys
                res[diag][impr_set][key] = []
                for i in collect(1:modelinfo["length"])
                    push!(res[diag][impr_set][key], full_dict["diag$(diag)_$(key)_set$(impr_set):[$i]"])
                end
            end
        end
    end
elseif extract_data == "xydata"
    fb = BDIO_open(joinpath(fittype_dict,"Fit",diag,comp,"xydata"),"r")
    partial_res = Vector{Dict}()
    full_dict = Dict{String, Any}()
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        push!(partial_res,ALPHAdobs_read_next(fb, size=sz, keys=ks))
    end
    BDIO_close!(fb)
    for dict in partial_res
        merge!(full_dict, dict)
    end
    res = Dict{String, Any}()
    res["xdata"] = full_dict["xdata"]
    for (diagIndex,diag) in enumerate(["a","b"])
        res[diag] = Dict{String, Dict}()
        for impr_set in ["1","2"]
            res[diag][impr_set] = Dict{String, Vector{uwreal}}()
            for key in mykeys
                res[diag][impr_set][key] = full_dict["$(key)_set$(impr_set)"][diagIndex,:]
            end
        end
    end
elseif extract_data == "MA"
    res = load(joinpath(fittype_dict,"MA",diag,comp,"MA_$diag.jld2"), "MA")
elseif extract_data == "info"
    res = load(joinpath(fittype_dict,"Fit",diag,comp,"ModelInfo.jld2"), "info")
end

res

##

myres = [subvec[1] for subvec in res["LO"]["1"]["g33_ll"]]; uwerr.(myres)

myres