# Import packages

using HVPobs
using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

using ProgressBars
using Suppressor

# Include Isovector Model (needed constants already included inside)

include("isovModel.jl")

include("HVPtools/Fit&MA.jl")

include("HVPtools/Reader.jl")

include("HVPtools/Utils.jl")

# BDIO path definition

julia_script_directory = @__DIR__

STD_DERIV = true

if STD_DERIV
    @info("Standard derivative has been chosen !!")
    path_bdio = joinpath(julia_script_directory, "..", "ObsBDIOstd")
else
    path_bdio = joinpath(julia_script_directory, "..", "ObsBDIO")
end

# Blind analysis (Simon K.) safe ensambles: 
# SU(3) symm. H101, B450, N202, N300
# H102, N101, C101, S400, N203, N200, D200, N302

ensList = ["A654","B450","C101","C102","D200","D201","D451","D452","E250","H101","H102","H200","J303","J304","J306","J307","J500","J501","N101","N200","N202","N203","N300","N302","N451","S400","F300"] # ,"E300"
ensInfo = EnsInfo.(ensList)

# We do not have charm or disconnected data for some of the ensambles
ensNOcharm = ["J501","N451","D150","D451","J304","C102","D251","D201","J306","J307","F300","H200"]
ensNOdisc  = ["D251","J306","J307","F300"]

DictComptoKey = Dict{String,Vector{String}}(
    "33"      => ["g33_ll","g33_lc"],
    "88"      => ["g88_ll","g88_lc"],

    # only interested in the lc (local-conserved) discr for the cc conn
    # "cc conn" => ["gcc_ll_conn","gcc_lc_conn"],
    # "cc conn" => ["gcc_lc_conn"],
    "cc conn" => ["gcc_lc_conn_beta"],

    # only interested in the cc (conserved-conserved) discr for the cc disc  & c8 disc
    "cc disc" => ["gcc_cc_disc"],
    "c8 disc" => ["gc8_cc_disc"],

    "3333"    => ["g3333_ll","g3333_lc"],
    "8888"    => ["g8888_ll","g8888_lc"],
    "CCCC"    => ["gCCCC_ll","gCCCC_lc"],
    "3388"    => ["g3388_ll","g3388_lc"],
    "33CC"    => ["g33CC_ll","g33CC_lc"],
    "88CC"    => ["g88CC_ll","g88CC_lc"]
)

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> Fits [LO, NLOa&b, NLOa, NLOb, NLOc] <==========================##

wind = "ID"  # NW  SD  ID  LD  ILD
diag = "LO"  # LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+)

base_model = "phi4"  # phi4  simple
perm_model = "All-(a4,a2loga)"  # All  All-(a4)  All-(a4,a2loga)  All-(a4,a2phi4)  ||  Only(phi2inv,logphi2)  Only(a3cutoff,phi2inv,logphi2)

IMPR_SET = ["1old","2"] # ["1","2"] ["1old","2"]

ensExcl = Dict(
    true => Dict(
        # "33"      => ["B450"],
        "33"      => ["D201"],

        # "88"      => ["B450","S400","N302"],
        "88"      => ["S400","N302","D201","J304","C102"],

        # "cc conn" => ["B450","S400","N302"]
        "cc conn" => ["B450"]
    ),
    false => Dict(
        "33"      => ["B450","N451","A654"], # J303
        "88"      => ["B450","S400","N302","J303"],
        "cc conn" => ["B450","S400","N302"]
    )
) # "C102","D450","D451","D251","D201","J304"

INFO = true

@info(" Reading and preparing data [diag. $diag; wind. $wind]")

t0 = Dict(); HVP = Dict(); FVC = Dict(); amu = Dict()
for ens in ensInfo
    println("- Reading data ensemble: $(ens.id)")

    println("   - Reading t0...")

    t0[ens.id] = BDIOread_t0(path_bdio, ens)

    println("   - Reading HVP...    [applying systematics]")

    HVP[ens.id] = Dict()
    for impr_set in IMPR_SET
        hvp, info = BDIOread_HVPens(path_bdio,diag,wind,ens,impr_set,info=true)
    
        HVP[ens.id][impr_set] = apply_syst_HVP(hvp,info["HVPsyst"],diag,wind,ens.id)
    end

    println("   - Reading FVC...    [applying systematics]")

    fvc = BDIOread_FVCens(path_bdio,diag,wind,ens)
    
    FVC[ens.id] = apply_syst_FVC(fvc,diag,wind,ens.id,IMPR_SET=IMPR_SET)

    println("   - aµ = HVP + FVC")

    amu[ens.id] = HVP_VolCorrect(HVP[ens.id],FVC[ens.id],diag,IMPR_SET=IMPR_SET)
end

@info(" ⟹ Data ready to be fitted\n\n")

## Data ready to fit

# comp = [33  88  cc conn  cc disc  c8 disc  ||  3333  8888  CCCC  3388  33CC  88CC]
if diag in ["LO", "NLOa", "NLOb", "NLOa&b", "NLOa&b(+)"]
    if wind in ["NW","SD"]
        COMP = ["33", "88", "cc conn", "cc disc", "c8 disc"]
    elseif wind in ["ID","LD","ILD"]
        COMP = ["33", "88", "cc conn"]
    end
elseif diag == "NLOc"
    COMP = ["3333", "8888", "CCCC", "3388", "33CC", "88CC"]
else
    error("Diag. $diag or wind. $wind not recognised. Please choose between the following options: \n   - diag.: LO  NLOa  NLOb  NLOc  NLOa&b  NLOa&b(+) \n   - wind: NW  SD  ID  LD  ILD")
end

@info(" STARTING FITS [diag. $diag; wind. $wind]")

for comp in [COMP[2]]
    mykeys = DictComptoKey[comp]
    fitIMPR_SET = comp in ["cc disc","c8 disc"] ? [IMPR_SET[1]] : IMPR_SET

    @info(" Fitting for component $comp")

    xdata = []
    ydata = Dict(); [ydata[impr_set] = Dict() for impr_set in IMPR_SET]  # initialize dict

    ensListFit = Vector{String}()

    println("- Creating 'x' & 'y' data points...")
    i = 0
    for ens in ensInfo
        if ((comp ∈ ["cc disc","c8 disc"]) && (ens.kappa_l == ens.kappa_s)) || ((comp == "cc conn") && (ens.id ∈ ensNOcharm)) || (comp == "88" && ens.id ∈ ensNOdisc) || (ens.id ∈ ensExcl[STD_DERIV][comp])
            println("   - [$(ens.id): Either excluded or 0 for contribution $comp]")
        else
            i += 1

            push!(xdata, [1 / (8*t0[ens.id]), 8*t0[ens.id]*m_ens[ens.id]["m_pi"]^2, 8*t0[ens.id]*(m_ens[ens.id]["m_K"]^2+0.5*m_ens[ens.id]["m_pi"]^2)])

            for impr_set in fitIMPR_SET
                for key in mykeys
                    i == 1 ? ydata[impr_set][key] = Vector{uwreal}() : nothing
                    push!(ydata[impr_set][key], amu[ens.id][impr_set][key])
                end
            end
            push!(ensListFit,ens.id)
        end
    end

    mdof = comp in ["cc disc","c8 disc"] ? 3 : 4    # impose a minimum to the dof for all fits 
    nens = length(ensListFit)                       # number of data ensembles

    f_tot_isov, n_par_tot_isov, label_tot_isov = call_models(base_model,perm_model,mdof,nens)

    println("- Fitting points")

    xdata = hcat(xdata...)
    xdata = [xdata[1,:] xdata[2,:] xdata[3,:]]

    fit = Dict{String, Dict}(); [fit[impr_set] = Dict{String, Vector{FitRes}}() for impr_set in IMPR_SET]  # initialize dict
    par = Dict{String, Dict}(); [par[impr_set] = Dict{String, Vector{Vector{uwreal}}}() for impr_set in IMPR_SET]  # initialize dict

    for impr_set in fitIMPR_SET
        println("   - Starting set "*impr_set)
        for key in mykeys
            println("      - Fitting for comp. $key ...")
            fit[impr_set][key] = Vector{FitRes}(); par[impr_set][key] = Vector{uwreal}()
            myFor = INFO ? collect(1:length(f_tot_isov)) : ProgressBar(collect(1:length(f_tot_isov)))
            for i in myFor
                myfit = fit_routine(f_tot_isov[i], value.(xdata), ydata[impr_set][key], n_par_tot_isov[i], pval=true, lineprint=INFO, ensList=ensListFit)
                push!(fit[impr_set][key], myfit)
                push!(par[impr_set][key], myfit.param)
            end
        end
    end

    println("- Printing BDIO & JDL2...")

    pFitMA = joinpath(path_bdio,"Fit&MA")

    pWind = joinpath(pFitMA,wind)
    !ispath(pWind) ? mkdir(pWind) : nothing

    pBase = joinpath(pWind,"base[$base_model]")
    !ispath(pBase) ? mkdir(pBase) : nothing

    pType = joinpath(pBase,perm_model)
    !ispath(pType) ? mkdir(pType) : nothing

    pFit = joinpath(pType,"Fit")
    !ispath(pFit) ? mkdir(pFit) : nothing

    pdiag = joinpath(pFit,diag)
    !ispath(pdiag) ? mkdir(pdiag) : nothing

    pcomp = joinpath(pdiag,comp)
    !ispath(pcomp) ? mkdir(pcomp) : nothing

    pFitRes = joinpath(pcomp,"FitRes.jld2")
    save(pFitRes,"FitRes",fit)

    io = IOBuffer()
    write(io, "parameters")

    fb = ALPHAdobs_create(joinpath(pcomp,"param"), io) 

    for i in collect(1:length(par[IMPR_SET[1]][mykeys[1]]))
        for impr_set in fitIMPR_SET
            parDict = Dict{String,Array{uwreal}}()
            for key in mykeys
                parDict["$(diag)_$(key)_set$(impr_set):[$i]"] = par[impr_set][key][i]
            end
            extra = Dict{String, Any}("impr_set" => impr_set, "diag" => diag, "wind" => wind)
            ALPHAdobs_write(fb, parDict, extra=extra)
        end
    end
    ALPHAdobs_close(fb)

    io = IOBuffer()
    write(io, "xydata")

    fb = ALPHAdobs_create(joinpath(pcomp,"xydata"), io)

    xDict = Dict{String,Array{uwreal}}("xdata" => xdata)
    ALPHAdobs_write(fb, xDict)
    for impr_set in fitIMPR_SET
        yDict = Dict{String,Array{uwreal}}()
        for key in mykeys
            yDict["$(key)_set$(impr_set)"] = ydata[impr_set][key]
        end
        extra = Dict{String, Any}("impr_set" => impr_set, "diag" => diag, "wind" => wind)
        ALPHAdobs_write(fb, yDict, extra=extra)
    end
    ALPHAdobs_close(fb)

    println("- Printing Model information...")

    infoDict = Dict{String,Any}(
        "length" => length(f_tot_isov),
        "nens" => nens,
        "ensList" => ensListFit,
        "n_par_tot_isov" => n_par_tot_isov,
        "label_tot_isov" => label_tot_isov,
    )

    pinfo = joinpath(pcomp,"ModelInfo.jld2")
    save(pinfo,"info",infoDict)

    println("")
end

## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> READING TEST <==========================##

base_model = "phi4"      #  phi4  simple  4-fit
perm_model = "All-(a4,a2loga)"   #  All  All-(a4)  All-(a4,a2loga)  All-(a4,a2phi4)  ||  Only(phi2inv,logphi2)  Only(a3cutoff,phi2inv,logphi2)
# mass_cut       = "<360"        #  all  <360  <300

impr_set = "1"
comp = "33"
diag = "NLOa&b"   # LO  NLOa  NLOb  NLOc  NLOa&b
extract_data = "param"   # FitRes  param  xydata  info

IMPR_SET = ["1old","2"]

pFit = joinpath(path_bdio,"DA",wind,"base[$base_model]",perm_model,"Fit",diag,comp)
if extract_data == "FitRes"
    res = load(joinpath(pFit,"FitRes.jld2"), "FitRes")
elseif extract_data == "param"
    modelinfo = load(joinpath(pFit,"ModelInfo.jld2"), "info")

    fb = BDIO_open(joinpath(pFit,"param"),"r")
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
        for impr_set in IMPR_SET
            res[diag][impr_set] = Dict{String, Vector{Vector{uwreal}}}()
            for key in mykeys
                res[diag][impr_set][key] = []
                for i in collect(1:modelinfo["length"])
                    push!(res[diag][impr_set][key], full_dict["$(diag)_$(key)_set$(impr_set):[$i]"])
                end
            end
        end
    end
elseif extract_data == "xydata"
    fb = BDIO_open(joinpath(pFit,"xydata"),"r")
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
        for impr_set in IMPR_SET
            res[diag][impr_set] = Dict{String, Vector{uwreal}}()
            for key in mykeys
                res[diag][impr_set][key] = full_dict["$(key)_set$(impr_set)"][diagIndex,:]
            end
        end
    end
elseif extract_data == "info"
    res = load(joinpath(pFit,"ModelInfo.jld2"), "info")
end

res

