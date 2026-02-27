
# t0, TMR  and corr reader

function BDIOread_t0(pBDIO::String,ensid::String)
    p = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_t0")
    t0_ = BDIOread_simple(p)
    return t0_["t0"][1]
end
BDIOread_t0(pBDIO::String,ens::EnsInfo) = BDIOread_t0(pBDIO,ens.id)

function BDIOread_t0_SU3sym(pBDIO::String,beta::Float64)
    beta_to_ens = Dict(
        3.34 => "A653",
        3.4  => "H101",
        3.46 => "B450",
        3.55 => "N202",
        3.7  => "J307",
        3.85 => "J500",
    )
    t0 = BDIOread_t0(pBDIO,beta_to_ens[beta])
    return t0
end
BDIOread_t0_SU3sym(pBDIO::String,ens::EnsInfo)  = BDIOread_t0_SU3sym(pBDIO,ens.beta)
BDIOread_t0_SU3sym(pBDIO::String,ensid::String) = BDIOread_t0_SU3sym(pBDIO,EnsInfo(ensid))

function BDIOread_TMR(pBDIO::String,ensid::String,diag::String;beta::Bool=false,resc::Bool=false,BLIND::Bool=false)
    SETstr  = resc ? "fPi" : "t0"
    SETstr *= beta ? (resc ? "ph" : "su3") : ""
    BLINstr = BLIND ? "Blind" : ""

    if diag ∉ ["LO","NLO","NLOa","NLOb","NLOa&b","NLOc","1D","NLO_1D","all"]
        error("diag not recognised, please choose between \n - LO, NLO, NLOa, NLOb, NLOa&b, NLOc, 1D, NLO_1D, all")
    end

    if diag == "LO"
        p = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_$(BLINstr)TMR$(SETstr)_LO")
        TMRDict = BDIOread_simple(p)
        TMR = TMRDict["TMR"]
    elseif diag in ["NLOa","NLOb","NLOa&b"]
        p = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_$(BLINstr)TMR$(SETstr)_NLOab")
        TMRDict =  BDIOread_general(p)
        if diag == "NLOa&b"
            TMR = TMRDict["TMRa"] .+ TMRDict["TMRb"]
        else
            TMR = TMRDict["TMR$(diag[end])"]
        end
    elseif diag == "NLOc"
        p = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_$(BLINstr)TMR$(SETstr)_NLOc")
        TMRDict = BDIOread_simple(p)
        TMR = TMRDict["TMRc"]
    elseif diag == "1D"
        pLO  = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_$(BLINstr)TMR$(SETstr)_LO")
        pNLOab = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_$(BLINstr)TMR$(SETstr)_NLOab")
        TMRLO  = BDIOread_simple(pLO)
        TMRNLOab = BDIOread_general(pNLOab)
        TMR = merge(TMRLO,TMRNLOab)
    elseif diag == "NLO_1D"
        pNLOab = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_$(BLINstr)TMR$(SETstr)_NLOab")
        TMRNLOab = BDIOread_general(pNLOab)
        TMR = TMRNLOab
    elseif diag == "all"
        pLO  = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_$(BLINstr)TMR$(SETstr)_LO")
        pNLOab = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_$(BLINstr)TMR$(SETstr)_NLOab")
        pNLOc  = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_$(BLINstr)TMR$(SETstr)_NLOc")
        TMRLO  = BDIOread_simple(pLO)
        TMRNLOab = BDIOread_general(pNLOab)
        TMRNLOc = BDIOread_general(pNLOc)
        TMR = merge(TMRLO,TMRNLOab,TMRNLOc)
    end
    return TMR
end
BDIOread_TMR(pBDIO::String,ens::EnsInfo,diag::String;beta::Bool=false,resc::Bool=false,BLIND::Bool=false) = BDIOread_TMR(pBDIO,ens.id,diag,beta=beta,resc=resc,BLIND=BLIND)

function BDIOread_corr(pBDIO::String,ensid::String,impr_set::String;STD::Bool=false)
    DERstr = STD ? "_std" : "" 
    p = joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_corr_set$(impr_set)$(DERstr)")
    corr = BDIOread_simple(p)
    return corr
end
BDIOread_corr(pBDIO::String,ens::EnsInfo,impr_set::String;STD::Bool=false)= BDIOread_corr(pBDIO,ens.id,impr_set,STD=STD)

function BDIOread_FVCcorr(pBDIO::String,ensid::String;Vref::Bool=false)
    p = Vref ? joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_FVC_Vref") : joinpath(pBDIO,"Corr&Kernel&t0",ensid,"$(ensid)_FVC")
    FVC = BDIOread_simple(p)
    return FVC
end
BDIOread_FVCcorr(pBDIO::String,ens::EnsInfo;Vref::Bool=false) = BDIOread_FVCcorr(pBDIO,ens.id,Vref=Vref)

function TXTread_FVCcorr_GS(pFVC_MLL::String,ensid::String;Vref::Bool=false)
    p = Vref ? joinpath("$pFVC_MLL","$(ensid)_gs_fvc_Vref.txt") : joinpath("$pFVC_MLL","$(ensid)_gs_fvc.txt")
    GS_data_matrix, headers = readdlm(p, '\t', header=true)

    t_gs   = GS_data_matrix[:, 1]
    fvc_gs = [uwreal([GS_data_matrix[:, 2][k],GS_data_matrix[:, 3][k]],ensid*"_gs") for k=1:length(t_gs)]
    return t_gs, fvc_gs
end
TXTread_FVCcorr_GS(pFVC_MLL::String,ens::EnsInfo;Vref::Bool=false) = TXTread_FVCcorr_GS(pFVC_MLL,ens.id,Vref=Vref)

# hvp and fvc (per ens) reader

function BDIOread_HVPens(pBDIO::String,diag::String,wind::String,ensid::String,impr_set::String;info::Bool=false,resc::Bool=false,STD::Bool=false,BLIND::Bool=false)
    DERstr  = STD ? "_std" : "" 
    RESCstr = resc ? "_resc" : ""
    BLINstr = BLIND ? "Blind" : ""
    p = joinpath(pBDIO,"HVP&FVC",wind,ensid,"$(ensid)_$(BLINstr)HVP$(diag)_set$(impr_set)$(RESCstr)$(DERstr)")
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
        info = load(joinpath(pBDIO,"HVP&FVC",wind,ensid,"$(ensid)_$(BLINstr)HVP$(diag)_info_set$(impr_set)$(RESCstr)$(DERstr).jld2"), "HVPinfo")
        return HVPens, info
    end
end
BDIOread_HVPens(pBDIO::String,diag::String,wind::String,ens::EnsInfo,impr_set::String;info::Bool=false,resc::Bool=false,STD::Bool=false,BLIND::Bool=false) = BDIOread_HVPens(pBDIO,diag,wind,ens.id,impr_set,info=info,resc=resc,STD=STD,BLIND=BLIND)

function BDIOread_FVCens(pBDIO::String,diag::String,wind::String,ensid::String;IMPR_SET::Vector{String}=["1","2"],resc::Bool=false,STD::Bool=false,Vref::Bool=false,BLIND::Bool=false)
    RESCstr = resc ? "_resc" : ""
    VREFstr = Vref  ? "_Vref" : ""
    BLINstr = BLIND ? "Blind" : ""
    if diag != "NLOc"
        p = joinpath(pBDIO,"HVP&FVC",wind,ensid,"$(ensid)_$(BLINstr)FVC$(diag)$(RESCstr)$(VREFstr)")
        
        if wind == "SDsub"
            FVCens  = BDIOread_general(p, merge=true)
            FVCens["FVCPi"] = FVCens["FVCPi"][1]; FVCens["FVCK"] = FVCens["FVCK"][1]; FVCens["FVC∆ls_amu"] = FVCens["FVC∆ls_amu"][1]; FVCens["FVC∆ls_amuconn"] = FVCens["FVC∆ls_amuconn"][1]
        else
            FVCens  = BDIOread_dim0(p)
        end
    elseif diag == "NLOc"
        FVCens = Dict()
        DERstr = STD ? "_std" : "" 
        for impr_set in IMPR_SET
            p = joinpath(pBDIO,"HVP&FVC",wind,ensid,"$(ensid)_$(BLINstr)FVC$(diag)_set$(impr_set)$(RESCstr)$(VREFstr)$(DERstr)")
            FVCens[impr_set] = BDIOread_dim0(p)
        end
    end
    return FVCens
end
BDIOread_FVCens(pBDIO::String,diag::String,wind::String,ens::EnsInfo;IMPR_SET::Vector{String}=["1","2"],resc::Bool=false,STD::Bool=false,Vref::Bool=false,BLIND::Bool=false) = BDIOread_FVCens(pBDIO,diag,wind,ens.id,IMPR_SET=IMPR_SET,resc=resc,STD=STD,Vref=Vref,BLIND=BLIND)

# fit reader

function BDIOread_res(pBDIO::String,diag::String,wind::String,comp::String,model_str::String,FitCUT::String,impr_set::String;resc::Bool=false,SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Vref::Bool=false,BLIND::Bool=false,Q::Float64=5.0,param::Bool=false)
    tl_str = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SIMstr = SimpleBase ? "SIMPLE" : ""
    MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
    VREFstr = Vref ? "_Vref" : ""
    RESstr = resc ? "_resc" : ""
    DERstr  = StdDer ? "_std" : "" 
    BLINstr = BLIND ? "Blind_" : ""
    pRes = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","$(BLINstr)Res$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr)")

    res = BDIOread_general(pRes)

    mykeys = DictComptoKey[comp]
    
    # modelinfo = JDL2read_ModelInfo(pBDIO,diag,wind,comp,model_str,FitCUT,impr_set;resc=resc,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Vref=Vref,BLIND=BLIND,Q=Q)

    # res = Dict{String,Array{uwreal}}(); par = Dict{String,Array{Array{uwreal}}}()
    # for key in mykeys
    #     res[key] = []; par[key] = []
    #     for i in collect(1:modelinfo["length"])
    #         push!(res[key], paramDict["$(key):[$i]"][1])
    #         push!(par[key], paramDict["$(key):[$i]"])
    #     end
    # end

    if param
        pPar = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","$(BLINstr)Param$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jdl2")
        parL = load(pPar, "Param")
        par  = Dict()
        for key in keys(parL)
            par[key] = []
            for n=1:length(parL[key])
                push!(par[key], [uwreal([p[1],p[2]],"p[$m] (fit $n, key $key)") for (m,p) in enumerate(parL[key][n])])
            end
        end
    end
    # if !isnothing(MultFunc) # this should not be used since we compute the res in the continuum from the ansatz itself and not the parameters
    #     x_ph = !resc ? [0 phi2_ph phi4_ph] : [0 y_ph z_ph]
    #     for key in mykeys
    #         res[key] = res[key] .* MultFunc(x_ph)[1]
    #     end
    # end
    param ? (return res, par) : (return res)
end
BDIOread_res(pBDIO::String,diag::String,wind::String,comp::String,model_var_list::Vector{Function},FitCUT::String,impr_set::String;resc::Bool=false,SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Vref::Bool=false,BLIND::Bool=false,Q::Float64=5.0,param::Bool=false) = BDIOread_res(pBDIO,diag,wind,comp,func_str(model_var_list),FitCUT,impr_set,resc=resc,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Vref=Vref,BLIND=BLIND,Q=Q,param=param)

function BDIOread_XYdata(pBDIO::String,diag::String,wind::String,comp::String,model_str::String,FitCUT::String,impr_set::String;resc::Bool=false,SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Vref::Bool=false,BLIND::Bool=false,Q::Float64=5.0)
    tl_str  = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SIMstr  = SimpleBase ? "SIMPLE" : ""
    MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
    VREFstr = Vref ? "_Vref" : ""
    RESstr = resc ? "_resc" : ""
    DERstr  = StdDer ? "_std" : ""
    BLINstr = BLIND ? "Blind_" : ""
    p = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","$(BLINstr)XYdata$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr)")

    mykeys = DictComptoKey[comp]

    data = BDIOread_general(p)

    Xdata = data["xdata"]
    Ydata = Dict()
    for key in mykeys
        Ydata[key] = data[key]
    end
    return Xdata, Ydata
end
BDIOread_XYdata(pBDIO::String,diag::String,wind::String,comp::String,model_var_list::Vector{Function},FitCUT::String,impr_set::String;resc::Bool=false,SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Vref::Bool=false,BLIND::Bool=false,Q::Float64=5.0) = BDIOread_XYdata(pBDIO,diag,wind,comp,func_str(model_var_list),FitCUT,impr_set,resc=resc,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Vref=Vref,BLIND=BLIND,Q=Q)

function JDL2read_FitRes(pBDIO::String,diag::String,wind::String,comp::String,model_str::String,FitCUT::String,impr_set::String;resc::Bool=false,SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Vref::Bool=false,BLIND::Bool=false,Q::Float64=5.0)
    tl_str = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SIMstr = SimpleBase ? "SIMPLE" : ""
    MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
    VREFstr = Vref ? "_Vref" : ""
    RESstr = resc ? "_resc" : ""
    DERstr  = StdDer ? "_std" : "" 
    BLINstr = BLIND ? "Blind_" : ""
    p = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","$(BLINstr)FitRes$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jld2")
    FitRes = load(p, "FitRes")
    return FitRes
end
JDL2read_FitRes(pBDIO::String,diag::String,wind::String,comp::String,model_var_list::Vector{Function},FitCUT::String,impr_set::String;resc::Bool=false,SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Vref::Bool=false,BLIND::Bool=false,Q::Float64=5.0) = JDL2read_FitRes(pBDIO,diag,wind,comp,func_str(model_var_list),FitCUT,impr_set,resc=resc,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Vref=Vref,BLIND=BLIND,Q=Q)

function JDL2read_ModelInfo(pBDIO::String,diag::String,wind::String,comp::String,model_str::String,FitCUT::String,impr_set::String;resc::Bool=false,SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Vref::Bool=false,BLIND::Bool=false,Q::Float64=5.0)
    tl_str = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SIMstr = SimpleBase ? "SIMPLE" : ""
    MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
    VREFstr = Vref ? "_Vref" : ""
    RESstr = resc ? "_resc" : ""
    DERstr  = StdDer ? "_std" : "" 
    BLINstr = BLIND ? "Blind_" : ""
    p = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","$(BLINstr)ModelInfo$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jld2")
    modelinfo = load(p, "Info")
    return modelinfo
end
JDL2read_ModelInfo(pBDIO::String,diag::String,wind::String,comp::String,model_var_list::Vector{Function},FitCUT::String,impr_set::String;resc::Bool=false,SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Vref::Bool=false,BLIND::Bool=false,Q::Float64=5.0) = JDL2read_ModelInfo(pBDIO,diag,wind,comp,func_str(model_var_list),FitCUT,impr_set,resc=resc,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Vref=Vref,BLIND=BLIND,Q=Q)

# MA reader

function BDIOread_MA(pBDIO::String,diag::String,wind::String,comp::String,model_str::String,FitCUT::String,impr_set::String;resc::Bool=false,SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Vref::Bool=false,BLIND::Bool=false,Q::Float64=5.0)
    tl_str = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SIMstr = SimpleBase ? "SIMPLE" : ""
    MULTstr = !isnothing(MultFunc) ? "($MultFunc)" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : ""
    VREFstr = Vref ? "_Vref" : ""
    RESstr = resc ? "_resc" : ""
    DERstr = StdDer ? "_std" : "" 
    BLINstr = BLIND ? "Blind" : ""
    pbdio = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","MA","$(BLINstr)MAres$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr)")
    pJDL2 = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"$(SIMstr)$(MULTstr)[$(model_str)]","Fit[$(FitCUT)]","MA","$(BLINstr)MAinfo$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jld2")
    
    res_tot, res = BDIOread_general(pbdio,merge=false)
    MA = Dict("res" => Dict{String,uwreal}() , "res_tot" => Dict{String,Array{uwreal}}() )
    MA["res_tot"] = res_tot
    [MA["res"][key] = res[key][1] for key in keys(res)]
    MAinfo = load(pJDL2,"MAinfo")
    return MA, MAinfo
end
BDIOread_MA(pBDIO::String,diag::String,wind::String,comp::String,model_var_list::Vector{Function},FitCUT::String,impr_set::String;resc::Bool=false,SimpleBase::Bool=false,MultFunc::Union{Nothing,Function}=nothing,StdDer::Bool=false,tlImpr::Bool=false,Vref::Bool=false,BLIND::Bool=false,Q::Float64=5.0) = BDIOread_MA(pBDIO,diag,wind,comp,func_str(model_var_list),FitCUT,impr_set,resc=resc,SimpleBase=SimpleBase,MultFunc=MultFunc,StdDer=StdDer,tlImpr=tlImpr,Vref=Vref,BLIND=BLIND,Q=Q)

function BDIOread_MAtot(pBDIO::String,diag::String,wind::String,comp::String;read::String="total",impr_set::Union{Nothing,String}=nothing,resc::Bool=false,StdDer::Bool=false,tlImpr::Bool=false,Vref::Bool=false,BLIND::Bool=false,Q::Float64=5.0)
    tl_str = (tlImpr && wind in ["SD","SDsub"] && comp in ["g33"]) ? "[tl]" : ""
    SUBQstr = (wind == "SDsub" && comp in ["g33","gCCconn","∆lc_b"]) ? "_Q$Q" : ""
    VREFstr = Vref ? "_Vref" : ""
    RESstr = resc ? "_resc" : ""
    DERstr = StdDer ? "_std" : "" 
    BLINstr = BLIND ? "Blind" : ""
    if read == "total"
        pbdio = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"MAtot","$(BLINstr)MARES$(SUBQstr)$(RESstr)$(VREFstr)$(DERstr)")
        pJDL2 = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"MAtot","$(BLINstr)MAINFO$(SUBQstr)$(RESstr)$(VREFstr)$(DERstr).jld2")

        res  = BDIOread_scalar(pbdio)
        info = load(pJDL2,"MAinfo")

        return res, info
    elseif read == "impr"
        isnothing(impr_set) ? error("impr_set must be given if 'read=impr'") : nothing
        IMPRstr = comp ∉ ["gCCdisc","gC8disc"] ? "_set$impr_set" : "_set"
        pbdio = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"MAtot","$(BLINstr)MAres$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr)")
        pJDL2 = joinpath(pBDIO,"Fit&MA",diag,wind,comp*tl_str,"MAtot","$(BLINstr)MAinfo$(SUBQstr)$(IMPRstr)$(RESstr)$(VREFstr)$(DERstr).jld2")

        res_tot, res = BDIOread_general(pbdio,merge=false)
        info = load(pJDL2,"MAinfo")

        MA = Dict("res" => Dict{String,uwreal}() , "res_tot" => Dict{String,Array{uwreal}}() )
        MA["res_tot"] = res_tot
        [MA["res"][key] = res[key][1] for key in keys(res)]
        return MA, info
    else
        error("Read flat '$read' not recognized; please choose between 'total' or 'impr_set'")
    end
end

# mass and charm kappa reader

function BDIOread_mPP(pBDIO::String,ensid::String)
    p = joinpath(pBDIO,"mass&dec",ensid,"$(ensid)_mPP")
    dict = BDIOread_dim0(p)
end
BDIOread_mPP(pBDIO::String,ens::EnsInfo) = BDIOread_mPP(pBDIO,ens.id)

function BDIOread_fPS(pBDIO::String,ensid::String)
    p = joinpath(pBDIO,"mass&dec",ensid,"$(ensid)_fPS")
    dict = BDIOread_dim0(p)
end
BDIOread_fPS(pBDIO::String,ens::EnsInfo) = BDIOread_fPS(pBDIO,ens.id)

function BDIOread_mDs(pBDIO::String)
    p  = joinpath(pBDIO,"mass&dec","mDs_prime")
    i = 0
    mDs_ph = 0.0
    mDs_SU3 = 0.0
    fb = BDIO_open(p,"r")
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        i=i+1
        if i==1
            mDs_ph = ALPHAdobs_read_next(fb)
        elseif i==2
            ks = collect(d["keys"])
            mDs_SU3 = ALPHAdobs_read_next(fb, keys=ks)
        end
    end
    BDIO_close!(fb)
    return mDs_ph, mDs_SU3
end

function BDIOread_mDs_kappaC(pBDIO::String,ensid::String)
    mDs_ph, mDs_SU3 = BDIOread_mDs(pBDIO)

    p_Ds_dict = joinpath(pBDIO,"mass&dec",ensid,"$(ensid)_mDsKappa")
    Ds_dict = Dict()

    fb = BDIO_open(p_Ds_dict,"r")
    i = 0
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        i=i+1
        if i==1
            sz = tuple(d["size"]...)
            Ds_dict["mDs"] = ALPHAdobs_read_next(fb, size=sz)
        else
            Ds_dict["kappaC"] = ALPHAdobs_read_next(fb)
        end
    end
    BDIO_close!(fb)
    return mDs_ph, mDs_SU3, Ds_dict
end
BDIOread_mDs_kappaC(pBDIO::String,ens::EnsInfo) = BDIOread_mDs_kappaC(pBDIO,ens.id)

# pQCD reader

function TXTread_bQ(pTXT::String,diag::String)
    p = joinpath(pTXT,"$(diag)_b33.txt")
    f_mat = readdlm(p, '\t', '\n', skipstart=3)
    QL  = f_mat[:,1]
    val = f_mat[:,2]
    unc = f_mat[:,3]
    return [uwreal([val[i],unc[i]],"b(Q=$(QL[i]))") for i=1:length(QL)]
end

# ChiPT FVC reader

function JDL2read_FVC_ChPT(pFVCcont::String,diag::String,wind::String;Q::Float64=5.0)
    p = joinpath(pFVCcont,"$(diag)_ChPT.jld2")
    fvc = load(p, "ChPT_FVC")
    if wind != "all"
        if wind != "SDsub"
            return fvc[wind]
        else
            return fvc[wind][Q.==Qlist][1]
        end
    else
        return fvc
    end
end