using ALPHAio, ADerrors
import ADerrors: err

using BDIO
using JLD2

include("Writer.jl")


# Function to extract Mathematica lists from a "contents" string

function TXTread_b33pert(path::String, diag::String; comp::String="g33")
    start_index = findfirst(x -> contains(x, "# $(list_name)"), split(contents, '\n'))
    if start_index === nothing
        error("List of coefficients $(an) has not been found")
    end

    start_index += 1
    end_index = findfirst(x -> isempty(x) || x[1] == '#', split(contents, '\n')[start_index:end])
    end_index = end_index === nothing ? length(contents) : start_index + end_index - 2

    raw_list = string.(split(contents, '\n')[start_index:end_index])
    
    return Baseparse.(string.(raw_list))
end

# General BDIO read structures

function BDIOread_simple(path::String)::Dict
    fb = BDIO_open(path,"r")
    mydict = Dict()
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        mydict = ALPHAdobs_read_next(fb, size=sz, keys=ks)
    end
    return mydict
end

function BDIOread_general(path::String;merge::Bool=true)::Union{Dict,Vector{Dict}}
    fb = BDIO_open(path,"r")
    mydict = Dict()
    mydictvec = Vector{Dict}()
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        sz = tuple(d["size"]...)
        ks = collect(d["keys"])
        push!(mydictvec,ALPHAdobs_read_next(fb, size=sz, keys=ks))
    end
    BDIO_close!(fb)
    if merge
        if length(mydictvec)>1
            for dict in mydictvec
                merge!(mydict, dict)
            end
        else
            mydict = mydictvec[1]
        end
    else
        mydict = mydictvec
    end
    return mydict
end

function BDIOread_scalar(path::String)::uwreal
    fb = BDIO_open(path,"r")
    mydict = Dict()
    while ALPHAdobs_next_p(fb)
        mydict = ALPHAdobs_read_next(fb)
    end
    return mydict
end

function BDIOread_dim0(path::String)::Dict
    fb = BDIO_open(path,"r")
    mydict = Dict()
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        ks = collect(d["keys"])
        mydict =  ALPHAdobs_read_next(fb, keys=ks)
    end
    return mydict
end

# BDIO read for computation phase

function BDIOread_t0(pBDIO::String,ensid::String)
    p = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_t0")
    t0_ = BDIOread_simple(p)
    return t0_["t0"][1]
end
BDIOread_t0(pBDIO::String,ens::EnsInfo) = BDIOread_t0(pBDIO,ens.id)

function BDIOread_TMR(pBDIO::String,ensid::String,diag::String;SU3::Bool=false)
    b = SU3 ? "SU3" : ""
    if diag == "LO"
        p = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_TMR$(b)_LO")
        TMRDict = BDIOread_simple(p)
        TMR = TMRDict["TMR"]
    elseif diag[1:3] == "NLO"
        p = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_TMR$(b)_NLO")
        TMRDict = BDIOread_general(p)
        if diag == "NLOa&b"
            TMR = TMRDict["TMRa"] .+ TMRDict["TMRb"]
        else
            TMR = TMRDict["TMR$(diag[end])"]
        end
    elseif diag == "all"
        pLO  = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_TMR$(b)_LO")
        pNLO = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_TMR$(b)_NLO")
        TMRLO  = BDIOread_simple(pLO)
        TMRNLO = BDIOread_general(pNLO)
        TMR = merge(TMRLO,TMRNLO)
    end
    return TMR
end
BDIOread_TMR(pBDIO::String,ens::EnsInfo,diag::String;SU3::Bool=false) = BDIOread_TMR(pBDIO,ens.id,diag,SU3=SU3)

function BDIOread_corr(pBDIO::String,ensid::String,impr_set::String;STD::Bool=false)
    DERstr = STD ? "_std" : "" 
    p = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_corr_set$(impr_set)$(DERstr)")
    corr = BDIOread_simple(p)
    return corr
end
BDIOread_corr(pBDIO::String,ens::EnsInfo,impr_set::String;STD::Bool=false)= BDIOread_corr(pBDIO,ens.id,impr_set,STD=STD)

function BDIOread_FVCcorr(pBDIO::String,ensid::String)
    p = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_FVC")
    FVC = BDIOread_simple(p)
    return FVC
end
BDIOread_FVCcorr(pBDIO::String,ens::EnsInfo) = BDIOread_FVCcorr(pBDIO,ens.id)

function BDIOread_KappaC_tar(pBDIO::String,ensid::String)
    p = joinpath(pBDIO,"kappaC_tar",ensid,"$(ensid)_kappaC")
    KappaC_tar =  BDIOread_scalar(p)
    return KappaC_tar
end
BDIOread_KappaC_tar(pBDIO::String,ens::EnsInfo) = BDIOread_KappaC_tar(pBDIO,ens.id)

# BDIO read for the Fit phase

function BDIOread_HVPens(pBDIO::String,diag::String,wind::String,ensid::String,impr_set::String;info::Bool=false,STD::Bool=false)
    DERstr = STD ? "_std" : "" 
    p = joinpath(pBDIO,"HVP&FVC",wind,ensid,"$(ensid)_HVP$(diag)_set$(impr_set)$(DERstr)")
    if wind == "SDsub"
        HVPens = BDIOread_general(p)
        for key in collect(keys(HVPens))[collect([length(HVPens[k]) == 1 for k in keys(HVPens)])]
            HVPens[key] = HVPens[key][1]
        end
    else
        HVPens = BDIOread_dim0(p)
    end
    if !info
        return HVPens
    else
        info = load(joinpath(pBDIO,"HVP&FVC",wind,ensid,"$(ensid)_HVP$(diag)_info_set$(impr_set)$(DERstr).jld2"), "HVPinfo")
        return HVPens, info
    end
end
BDIOread_HVPens(pBDIO::String,diag::String,wind::String,ens::EnsInfo,impr_set::String;info::Bool=false,STD::Bool=false) = BDIOread_HVPens(pBDIO,diag,wind,ens.id,impr_set,info=info,STD=STD)

function BDIOread_FVCens(pBDIO::String,diag::String,wind::String,ensid::String;IMPR_SET::Vector{String}=["1","2"],STD::Bool=false)
    if diag != "NLOc"
        p = joinpath(pBDIO,"HVP&FVC",wind,ensid,"$(ensid)_FVC$(diag)")
        if wind == "SDsub"
            FVCens  = BDIOread_general(p, merge=true)
            FVCens["FVCPi"] = FVCens["FVCPi"][1]; FVCens["FVCK"] = FVCens["FVCK"][1]; FVCens["FVC∆ls_amu"] = FVCens["FVC∆ls_amu"][1]
        else
            FVCens = BDIOread_simple(p)
        end
    elseif diag == "NLOc"
        DERstr = STD ? "_std" : "" 
        for impr_set in IMPR_SET
            p = joinpath(pBDIO,"HVP&FVC",wind,ensid,"$(ensid)_FVC$(diag)_set$(impr_set)$(DERstr)")
            FVCens[impr_set] = BDIOread_simple(p)
        end
    end
    return FVCens
end
BDIOread_FVCens(pBDIO::String,diag::String,wind::String,ens::EnsInfo;IMPR_SET::Vector{String}=["1","2"],STD::Bool=false) = BDIOread_FVCens(pBDIO,diag,wind,ens.id,IMPR_SET=IMPR_SET,STD=STD)


function JDL2read_3l0(pBDIO::String,diag::String)
    p = joinpath(pBDIO,"HVP&FVC","SDsub","HVP$(diag)_3l0.jld2")
    HVP3l0 = load(p, "HVP3l0")
    return Float64.(HVP3l0)
end


function BDIOread_XYdata(pBDIO::String,diag::String,wind::String,comp::String,model_str::String,FitCUT::String,impr_set::String;SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Q::Float64=5.0)
    tl_str = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SIMstr = SimpleBase ? "SIMPLE" : ""
    MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
    DERstr = StdDer ? "_std" : "" 
    p = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","xydata$(SUBQstr)$(IMPRstr)$(DERstr)")

    mykeys = DictComptoKey[comp]

    data = BDIOread_general(p)

    Xdata = data["xdata"]
    Ydata = Dict()
    for key in mykeys
        Ydata[key] = data[key]
    end
    return Xdata, Ydata
end
BDIOread_XYdata(pBDIO::String,diag::String,wind::String,comp::String,model_var_list::Vector{Function},FitCUT::String,impr_set::String;SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Q::Float64=5.0) = BDIOread_XYdata(pBDIO,diag,wind,comp,func_str(model_var_list),FitCUT,impr_set,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Q=Q)

function JDL2read_FitRes(pBDIO::String,diag::String,wind::String,comp::String,model_str::String,FitCUT::String,impr_set::String;SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Q::Float64=5.0)
    tl_str = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SIMstr = SimpleBase ? "SIMPLE" : ""
    MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
    DERstr = StdDer ? "_std" : "" 
    p = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","FitRes$(SUBQstr)$(IMPRstr)$(DERstr).jld2")
    FitRes = load(p, "FitRes")
    return FitRes
end
JDL2read_FitRes(pBDIO::String,diag::String,wind::String,comp::String,model_var_list::Vector{Function},FitCUT::String,impr_set::String;SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Q::Float64=5.0) = JDL2read_FitRes(pBDIO,diag,wind,comp,func_str(model_var_list),FitCUT,impr_set,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Q=Q)

function JDL2read_ModelInfo(pBDIO::String,diag::String,wind::String,comp::String,model_str::String,FitCUT::String,impr_set::String;SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Q::Float64=5.0)
    tl_str = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SIMstr = SimpleBase ? "SIMPLE" : ""
    MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
    DERstr = StdDer ? "_std" : "" 
    p = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","ModelInfo$(SUBQstr)$(IMPRstr)$(DERstr).jld2")
    modelinfo = load(p, "info")
    return modelinfo
end
JDL2read_ModelInfo(pBDIO::String,diag::String,wind::String,comp::String,model_var_list::Vector{Function},FitCUT::String,impr_set::String;SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Q::Float64=5.0) = JDL2read_ModelInfo(pBDIO,diag,wind,comp,func_str(model_var_list),FitCUT,impr_set,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Q=Q)

function BDIOread_res(pBDIO::String,diag::String,wind::String,comp::String,model_str::String,FitCUT::String,impr_set::String;SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Q::Float64=5.0,param::Bool=false)
    tl_str = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SIMstr = SimpleBase ? "SIMPLE" : ""
    MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
    DERstr = StdDer ? "_std" : "" 
    p = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","param$(SUBQstr)$(IMPRstr)$(DERstr)")

    mykeys = DictComptoKey[comp]
    
    modelinfo = JDL2read_ModelInfo(pBDIO,diag,wind,comp,model_str,FitCUT,impr_set;SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Q=Q)

    paramDict = BDIOread_general(p)

    res = Dict{String,Array{uwreal}}(); par = Dict{String,Array{Array{uwreal}}}()
    for key in mykeys
        res[key] = []; par[key] = []
        for i in collect(1:modelinfo["length"])
            push!(res[key], paramDict["$(key):[$i]"][1])
            push!(par[key], paramDict["$(key):[$i]"])
        end
    end
    if !isnothing(MultFunc)
        x_ph = [0 phi2_ph phi4_ph]
        for key in mykeys
            res[key] = res[key] .* MultFunc(x_ph)[1]
        end
    end
    param ? (return res, par) : (return res)
end
BDIOread_res(pBDIO::String,diag::String,wind::String,comp::String,model_var_list::Vector{Function},FitCUT::String,impr_set::String;SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Q::Float64=5.0,param::Bool=false) = BDIOread_res(pBDIO,diag,wind,comp,func_str(model_var_list),FitCUT,impr_set,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Q=Q,param=param)

function BDIOread_MA(pBDIO::String,diag::String,wind::String,comp::String,model_str::String,FitCUT::String,impr_set::String;SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Q::Float64=5.0)
    tl_str = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SIMstr = SimpleBase ? "SIMPLE" : ""
    MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
    DERstr = StdDer ? "_std" : "" 
    pbdio = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","MA","MAres$(SUBQstr)$(IMPRstr)$(DERstr)")
    pJDL2 = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","MA","MAinfo$(SUBQstr)$(IMPRstr)$(DERstr).jld2")
    
    res_tot, res = BDIOread_general(pbdio,merge=false)
    MA = Dict("res" => Dict{String,uwreal}() , "res_tot" => Dict{String,Array{uwreal}}() )
    MA["res_tot"] = res_tot
    [MA["res"][key] = res[key][1] for key in keys(res)]
    MAinfo = load(pJDL2,"MAinfo")
    return MA, MAinfo
end
BDIOread_MA(pBDIO::String,diag::String,wind::String,comp::String,model_var_list::Vector{Function},FitCUT::String,impr_set::String;SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Q::Float64=5.0) = BDIOread_MA(pBDIO,diag,wind,comp,func_str(model_var_list),FitCUT,impr_set,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Q=Q)

function BDIOread_MAtot(pBDIO::String,diag::String,wind::String,comp::String;read::String="total",impr_set::Union{Nothing,String}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Q::Float64=5.0)
    tl_str = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    DERstr = StdDer ? "_std" : "" 
    if read == "total"
        pbdio = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"MAtot","MARES$(SUBQstr)$(DERstr)")
        pJDL2 = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"MAtot","MAINFO$(SUBQstr)$(DERstr).jld2")
    elseif read == "impr"
        isnothing(impr_set) ? error("impr_set must be given if 'read=impr'") : nothing
        IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
        pbdio = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"MAtot","MAres$(SUBQstr)$(IMPRstr)$(DERstr)")
        pJDL2 = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"MAtot","MAinfo$(SUBQstr)$(IMPRstr)$(DERstr).jld2")
    else
        error("Read flat '$read' not recognized; please choose between 'total' or 'impr_set'")
    end

    res  = BDIOread_scalar(pbdio)
    info = load(pJDL2,"MAinfo")

    return res, info
end

function BDIOread_mDs(pBDIO::String)
    p = joinpath(pBDIO,"kappaC_tar","mDs_prime")
    mDs = Dict()
    fb = BDIO_open(p,"r")
    i = 0
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        i=i+1
        if i==1
            mDs["mDs ph"] = ALPHAdobs_read_next(fb)
        else
            ks = collect(d["keys"])
            mDs["amDs SU3"] = ALPHAdobs_read_next(fb, keys=ks)
        end
    end
    return mDs
end

