
function krnl_dα(t::Union{Float64, uwreal}, Q2::Union{Float64,uwreal})

    K = t^2 - 4/Q2 * sin((sqrt(Q2)*t)/2)^2
    return K
end

function krnl_dα_sub(t::Union{Float64, uwreal}, Q2::Union{Float64,uwreal}, Q2m::Union{Float64,uwreal})

    K = t^2 - 4/Q2 * sin((sqrt(Q2)*t)/2)^2 - 4/3 * (Q2/Q2m^2) * sin((sqrt(Q2m)*t)/2)^4
    return K
end

function krnl_dα_qhalf(t::Union{Float64, uwreal}, Q2::Union{Float64,uwreal})

    K = 16. / Q2 * sin((sqrt(Q2)*t)/4)^4
    return K
end

function krnl_dα_qhalf_sub(t::Union{Float64, uwreal}, Q2::Union{Float64,uwreal}, Q2m::Union{Float64,uwreal})
    K = 16. / Q2 * sin((sqrt(Q2)*t)/4)^4 - (Q2/Q2m^2) * sin((sqrt(Q2m)*t)/2)^4
    return K
end

function tmr_integrand(cc::Vector{uwreal}, q2::Union{Float64, uwreal}, krnl::Function; pl::Bool=true, t0ens::Union{uwreal, Nothing}=nothing, wind::Union{Nothing, Window}=nothing, data::Bool=false)
    
    if isnothing(t0ens) && !isnothing(wind)
        error("You assigned a window without passign t0ens. \n Conversion to t[fm] failed.")
    end
    T = length(cc)
    x0 = Float64.(collect(0:T-1)) 
    k = krnl.(x0, q2) #.* value.(sqrt.(t0ens) ./ t0sqrt_ph ).^(-2)
    
    integrand = cc .* k  #abs.(cc .* k)
    x0 = x0 .* value(t0sqrt_ph / sqrt(t0ens))
    if !isnothing(wind)
        integrand = integrand .* wind(x0)
    end


    Thalf = try
        Int64(T/2)
   catch
        Int64((T+1)/2)
   end
    if pl        
        uwerr.(integrand)
        if !isnothing(t0ens)
            # x0 = x0 .* value(t0sqrt_ph / sqrt(t0ens))
            xlabel(L"$t \ [\mathrm{fm}]$")
        else
            xlabel(L"$t/a$")
        end
        errorbar(x0[1:Thalf], value.(integrand)[1:Thalf], err.(integrand)[1:Thalf], fmt="s", mfc="none", color="forestgreen", capsize=2)
        ylabel(L"$G(t) \cdot K(t, Q^2)$")
        xlim(-0.05, 4)
        display(gcf())
        close("all")
    end
    !data ? (return sum(integrand[1:Thalf])) : (return sum(integrand[1:Thalf]), integrand[1:Thalf])
end
tmr_integrand(cc::Corr, q2::Union{Float64, uwreal}, krnl::Function; pl::Bool=true, t0ens::Union{uwreal, Nothing}=nothing, wind::Union{Nothing, Window}=nothing, data::Bool=false) = tmr_integrand(cc.obs, q2, krnl; pl=pl, t0ens=t0ens, wind=wind, data=data)

function tmr_integrand(cc::Vector{uwreal}, q2::Union{Float64, uwreal}, q2m::Union{Float64, uwreal}, krnl::Function; pl::Bool=true, t0ens::Union{uwreal, Nothing}=nothing, wind::Union{Nothing, Window}=nothing, data::Bool=false)
    if isnothing(t0ens) && !isnothing(wind)
        error("You assigned a window without passign t0ens. \n Conversion to t[fm] failed.")
    end
    T = length(cc)
    x0 = Float64.(collect(0:T-1)) 
    k = krnl.(x0, q2, q2m) #.* value.(sqrt.(t0ens) ./ t0sqrt_ph ).^(-2)
    
    integrand = cc .* k # abs.(cc .* k)
    x0 = x0 .* value(t0sqrt_ph / sqrt(t0ens))
    if !isnothing(wind)
        integrand = integrand .* wind(x0)
    end


    Thalf = try
        Int64(T/2)
   catch
        Int64((T+1)/2)
   end
    if pl        
        uwerr.(integrand)
        if !isnothing(t0ens)
            # x0 = x0 .* value(t0sqrt_ph / sqrt(t0ens))
            xlabel(L"$t \ [\mathrm{fm}]$")
        else
            xlabel(L"$t/a$")
        end
        errorbar(x0[1:Thalf], value.(integrand)[1:Thalf], err.(integrand)[1:Thalf], fmt="s", mfc="none", color="forestgreen", capsize=2)
        ylabel(L"$G(t) \cdot K(t, Q^2)$")
        xlim(-0.05, 4)
        display(gcf())
        close("all")
    end
    !data ? (return sum(integrand[1:Thalf])) : (return sum(integrand[1:Thalf]), integrand[1:Thalf])
end
tmr_integrand(cc::Corr, q2::Union{Float64, uwreal}, q2m::Union{Float64, uwreal}, krnl; pl::Bool=true, t0ens::Union{uwreal, Nothing}=nothing, wind::Union{Nothing, Window}=nothing, data::Bool=false) = tmr_integrand(cc.obs, q2, q2m, krnl; pl=pl, t0ens=t0ens, wind=wind, data=data)


function tmr_integrand_3l(cc::Vector{Tu}, q2::Union{Tu, Tf}, q2m::Union{Tu, Tf}, krnl::Function; pl::Bool=true, t0ens::Union{Tu, Nothing}=nothing, wind::Bool=false, d=0.4, delta=0.15 ) where {Tu, Tf}
    
    T = length(cc)

    x0 = Float64.(collect(0:T-1)) 
    k = krnl.(x0, q2, q2m) #.* value.(sqrt.(t0ens) ./ t0sqrt_ph ).^(-2)
    x0 = Float64.(collect(0:T-1)) .* value(t0sqrt_ph / sqrt(t0ens))
    
    integrand = cc .* k
    
    Thalf = try
        Int64(T/2)
   catch
        Int64((T+1)/2)
   end

    if pl        
        uwerr.(integrand)

        cumtmr = cumsum(integrand)
        plot(x0, value.(cumtmr))
        xlabel(L"$t \ [\mathrm{fm}]$")
        ylabel(L"$\sum G(t) \cdot K(t, Q^2)$")
        if wind
            axvline(d, ls="dashed", color="black", lw=0.8)
        end
        display(gcf())
        close()

        errorbar(x0[1:Thalf], abs.(value.(integrand))[1:Thalf], err.(integrand)[1:Thalf], fmt="s", mfc="none", color="forestgreen", capsize=2)
        xlabel(L"$t/a$")
        ylabel(L"$G(t) \cdot K(t, Q^2)$")
        if wind
            axvline(d, ls="dashed", color="black", lw=0.8)
        end
        display(gcf())
        close("all")
    end

    if wind
        w = window_krnl.(x0, d=d, delta=delta)
        integrand .*=  w
    end
    return   integrand[1:Thalf] * sign(value(integrand[5]))
end

function window_krnl(x0; d=0.4, delta=0.15)
    w = 1 - 0.5 * (1 + tanh((x0-d)/delta))
    return w
end

function treelevel_continuum_correlator(t)
    return 1. / (2* pi^2 * t^3)
end

function get_meff_BMA(corr::Vector{Corr}, ens::EnsInfo; pl::Bool=true, 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, path_plt::Union{String,Nothing}=nothing, ll::LaTeXString=L"$m_{\pi}$")
   
    M = Vector{uwreal}(undef, length(corr))
    if !isnothing(path_plt) 
        path_plt = joinpath(path_plt, ens.id)  
    end

    
 
    const_fit = 0
    if ens.id in ["A653", "A654", "D150", "B450", "N451", "N452", "D450", "D451", "D452", "D251", "E250"]
        taux = collect(10:40) 
        Thalf = Int(length(corr[1].obs)/2)+1
        T = length(corr[1].obs)
        tmin = Thalf .- taux
        tmax = Thalf .+ taux
        tmax = tmax[1:1:end]
        tmin = tmin[1:1:end] 
        @. const_fit_pbc(x,p) = p[2] * (exp(-p[1]*(x)) + exp(-p[1]*(T-x)) )# + p[4] * (exp(-p[3]*(x)) + exp(-p[3]*(T-x)) )
        param=2
        pbc = true
        const_fit = const_fit_pbc
    else
        @. const_fit_obc(x, p) = p[1] + 0 * x
        param = 1

        # tmin and tmax tuned on H101 and then scaled 
        tminaux = collect(50:70) #* round(Int64, sqrt(value(t0(ens.beta) / t0(3.4))))
        T = length(corr[1].obs)
        tmaxaux = collect(50:70) #* round(Int64, sqrt(value(t0(ens.beta) / t0(3.4))))
        tmax = T .- collect(reverse(tmaxaux))
        tmax = tmax[1:2:end]
        tmin = tminaux[1:2:end] 
        pbc = false
        const_fit = const_fit_obc
    end

    println("tmax: ", tmax)
    println("tmin: ", tmin)

    for k in eachindex(corr)
        t = string(L"$\kappa_1 = $" ,L" $\kappa_2 = $" )   
        if ens.id in ["A653", "A654","D150", "B450", "N451", "N452", "D450", "D451", "D452", "D251", "E250"]
            # data = -1. .* corr[k].obs
            data = -1. .* corr[k].obs
        else
            _, data = meff(corr[k], [tmin[end-4],tmax[1]], pl=false, wpm=wpm, data=true)
        end

        name_bma = ll*"_BMA.pdf"
        isnothing(path_plt) ? tt = nothing : tt = joinpath(path_plt, name_bma)

        M[k], syst, all_res, FitDict = bayesian_av(const_fit, data, tmin, tmax, param, pl=pl, wpm=wpm, data=true, path_plt=tt, plt_title=t , label=ll, pbc=pbc)
        println(FitDict["AIC"])
        if pl
            fig = figure(figsize=(10,7.5))
            modstot = FitDict["mods"]
            AIC  = FitDict["AIC"]
            Wtot = FitDict["W"]
            IDX  = partialsortperm(AIC, 1:2)
            mods = modstot[IDX]
            main_res = all_res[IDX]
            W = Wtot[IDX]
            isnothing(wpm) ? uwerr.(data) : [uwerr(data[i],wpm) for i in eachindex(data)]
            vv = value.(data)
            ee = err.(data)
            xx = collect(0:length(vv)-1) 
            errorbar(xx, vv, ee, fmt="d", mfc="none", color="black", capsize=2)
            CC = ["forestgreen", "royalblue", "tomato"]
            for (ii,trange) in enumerate(mods)
                isnothing(wpm) ? uwerr(main_res[ii]) : uwerr(main_res[ii], wpm)
                vaux = value(main_res[ii])
                eaux = err(main_res[ii])
                regex = r"([0-9]+),([0-9]+)"
                aux = match(regex, trange)
                fill_between(parse(Float64,aux[1]):parse(Float64,aux[2]), vaux-eaux, vaux+eaux, color=CC[ii], alpha=0.5, label=round(W[ii], digits=3))
            end
            
            axhline(y=value(main_res[1])-14*err(main_res[1]), xmin=(tmin[1])/length(xx), xmax=(tmin[end])/length(xx), ls="dashed", color="green")
            axhline(y=value(main_res[1])-14.3*err(main_res[1]), xmin=(tmax[1])/length(xx), xmax=(tmax[end])/length(xx), ls="dashed", color="red")
            text(tmin[1],value(main_res[1])-13*err(main_res[1]), L"$\{t_{\mathrm{min}}\}$" )
            text(tmax[1],value(main_res[1])-13*err(main_res[1]), L"$\{t_{\mathrm{max}}\}$" )

            legend(title="Weights")
            xlabel(L"$t/a$")
            ylabel(ll)
            xlim(xx[1], xx[end])
            ylim(value(main_res[1]) -20*err(main_res[1]), value(main_res[1]) +20*err(main_res[1]))
            twiny()
            xlim(xx[1], xx[end])
            fm_time = round.((xx)[1:10:end] .* value(a(ens.beta)), digits=2)
            xticks(xx[1:10:end], fm_time)
            # xticks()
            xlabel(L"$t-t_0\ [\mathrm{fm}]$", fontsize=14)
            display(fig)
            if !isnothing(path_plt)
                savefig(joinpath(path_plt, ll*"_plateau.pdf"))
            end
            close("all")
        end
    end
    return M
end


function get_dec_const_BMA(corr_pp::Corr, corr_a0p::Corr, mps::uwreal, ens::EnsInfo; pl::Bool=true, 
    wpm::Union{Dict{Int64,Vector{Float64}},Dict{String,Vector{Float64}}, Nothing}=nothing, path_plt::Union{String,Nothing}=nothing, ll::LaTeXString=L"$m_{\pi}$")
   
    if !isnothing(path_plt) 
        path_plt = joinpath(path_plt, ens.id)  
    end

    
 
    if ens.id in ["A653", "A654","B450", "E250"]
        taux = collect(15:40) 
        Thalf = Int(length(corr_pp.obs)/2)+1
        T = length(corr_pp.obs)
        tmin = Thalf .- taux
        tmax = Thalf .+ taux
        tmax = tmax[1:2:end]
        tmin = tmin[1:2:end] 
        @. fit_corr(x,p) = p[1] * (exp(-p[2]*(x)) + exp(-p[2]*(T-x)) )
        @. fit_corr_a0p(x,p) = p[1] * (exp(-p[2]*(x)) - exp(-p[2]*(T-x)) )

        param=2
        pbc = true
    else
        #@. const_fit(x, p) = p[1] + 0 * x
        param = 1 

        # tmin and tmax tuned on H101 and then scaled 
        tminaux = collect(28:55) #* round(Int64, sqrt(value(t0(ens.beta) / t0(3.4))))
        T = length(corr_pp.obs)
        tmaxaux = collect(30:60) #* round(Int64, sqrt(value(t0(ens.beta) / t0(3.4))))
        tmax = T .- collect(reverse(tmaxaux))
        tmax = tmax[1:2:end]
        tmin = tminaux[1:2:end] 
        pbc = false
    end

    println("tmax: ", tmax)
    println("tmin: ", tmin)

    t = string(L"$\kappa_1 = $" ,L" $\kappa_2 = $" )   
    name_bma = ll*"_BMA.pdf"
    isnothing(path_plt) ? tt = nothing : tt = joinpath(path_plt, name_bma)
    if ens.id in [ "A653", "A654","B450", "E250"]
        data_pp = -1. .* corr_pp.obs
        data_a0p =  corr_a0p.obs
        idx = findall(x->x==0.0, value.(data_a0p))
        data_a0p[idx[2]] = uwreal([0.0, 1e-10], "test")
        idx = findall(x->x==0.0, value.(data_a0p))
        uwerr.(data_a0p)
        errorbar(collect(1:T), value.(data_a0p), err.(data_a0p), fmt="s")
        display(gcf())
        close("all")
        cp, syst, all_res, FitDict = bayesian_av(fit_corr, data_pp, tmin, tmax, param, pl=pl, wpm=wpm, data=true, path_plt=tt, plt_title=t , label=ll, pbc=pbc)
        caa, syst, all_res, FitDict = bayesian_av(fit_corr_a0p, data_a0p, tmin, tmax, param, pl=pl, wpm=wpm, data=true, path_plt=tt, plt_title=t , label=ll, pbc=pbc)

        f =  sqrt(2) / sqrt(mps) * caa / sqrt(cp)
    else
        _, data = meff(corr[k], [tmin[end-4],tmax[1]], pl=false, wpm=wpm, data=true)
    end


    if pl
        data = data_pp
        fig = figure(figsize=(10,7.5))
        modstot = FitDict["mods"]
        AIC  = FitDict["AIC"]
        Wtot = FitDict["W"]
        IDX  = partialsortperm(AIC, 1:3)
        mods = modstot[IDX]
        main_res = all_res[IDX]
        W = Wtot[IDX]
        isnothing(wpm) ? uwerr.(data) : [uwerr(data[i],wpm) for i in eachindex(data)]
        vv = value.(data)
        ee = err.(data)
        xx = collect(0:length(vv)-1) 
        errorbar(xx, vv, ee, fmt="d", mfc="none", color="black", capsize=2)
        CC = ["forestgreen", "royalblue", "tomato"]
        for (ii,trange) in enumerate(mods)
            isnothing(wpm) ? uwerr(main_res[ii]) : uwerr(main_res[ii], wpm)
            vaux = value(main_res[ii])
            eaux = err(main_res[ii])
            regex = r"([0-9]+),([0-9]+)"
            aux = match(regex, trange)
            fill_between(parse(Float64,aux[1]):parse(Float64,aux[2]), vaux-eaux, vaux+eaux, color=CC[ii], alpha=0.5, label=round(W[ii], digits=3))
        end
        
        axhline(y=value(main_res[1])-14*err(main_res[1]), xmin=(tmin[1])/length(xx), xmax=(tmin[end])/length(xx), ls="dashed", color="green")
        axhline(y=value(main_res[1])-14.3*err(main_res[1]), xmin=(tmax[1])/length(xx), xmax=(tmax[end])/length(xx), ls="dashed", color="red")
        text(tmin[1],value(main_res[1])-13*err(main_res[1]), L"$\{t_{\mathrm{min}}\}$" )
        text(tmax[1],value(main_res[1])-13*err(main_res[1]), L"$\{t_{\mathrm{max}}\}$" )
        legend(title="Weights")
        xlabel(L"$t/a$")
        ylabel(ll)
        xlim(xx[1], xx[end])
        ylim(value(main_res[1]) -20*err(main_res[1]), value(main_res[1]) +20*err(main_res[1]))
        twiny()
        xlim(xx[1], xx[end])
        fm_time = round.((xx)[1:12:end] .* value(a(ens.beta)), digits=2)
        xticks(xx[1:12:end], fm_time)
        # xticks()
        xlabel(L"$t-t_0\ [\mathrm{fm}]$", fontsize=14)
        display(fig)
        if !isnothing(path_plt)
            savefig(joinpath(path_plt, ll*"_plateau.pdf"))
        end
        close("all")
    end
    return f
end


function get_TIC(fitcat::FitCat)
    chi2tot = getfield.(fitcat.fit, :chi2)
    chi2exptot = getfield.(fitcat.fit, :chi2exp)
    ic = chi2tot .- 2 .* chi2exptot
    return ic
end

function get_bma_weight(ic::Vector{Float64})
    w = exp.(-0.5 .* ic)
    w = w ./ sum(w)
    return w
end

function get_w_from_fitcat(catvec::Vector{FitCat}; norm::Bool=false)
    ic = get_TIC.(catvec)
    ic = vcat(ic...)
    w = get_bma_weight(ic)
    if norm
        w_min = minimum(w)
        w_max = maximum(w)
        w = 1.0 .* (w .- w_min) ./ (w_max - w_min) .+ 0.0 # normalise between 0.1 and 1
        replace!(x -> x>1 ? 1 : x, w)
    end
    return w
end

function model_average(results::Vector{uwreal}, ww)
    ww = ww ./ sum(ww)
    fin_res = sum(ww .* results)  
    aux1 = sum( ww .* results.^2) 
    aux2 = sum( ww .* results)^2 
    syst = sqrt(abs(value(aux1 - aux2)))
    #syst = sum( )
    return fin_res, syst
    
end

function set_fluc_to_zero!(a::uwreal, id_str::String)
    ws = ADerrors.wsg
    id_int = ws.str2id[id_str]
    idx = ADerrors.find_mcid(a, id_int)
    if isnothing(idx)
        error("No error available... maybe run uwerr")
    else
        nd = ws.fluc[ws.map_ids[a.ids[idx]]].nd
        for j in 1:length(a.prop)
            if (a.prop[j] && ((ws.map_nob[j] == a.ids[idx])))
                ws.fluc[j].delta[:] .= 0.0
            end
        end
        return nothing
    end
    
end

function corrBounding(corr::Vector{uwreal}, mUp::uwreal, mLow::uwreal, tcut::Int64, ens::EnsInfo)
    T = HVPobs.Data.get_T(ens.id)

    upperBound = corr[:]
    lowerBound = corr[:]

    # @. G_PBC(x0, p) = exp(-p*(x0-tcut)) + exp(-p*(T-x0+tcut))
    # @. G_OBC(x0, p) = exp(-p*(x0-tcut))
    @. G_PBC(x0, p) = exp(-p*(x0)) + exp(-p*(T-x0))
    @. G_OBC(x0, p) = exp(-p*(x0))


    if ens.bc == "obc"
        upperBound[tcut:end] = corr[tcut] ./ G_OBC(tcut, mUp) .* G_OBC(collect(tcut:T), mUp)
        lowerBound[tcut:end] = corr[tcut] ./ G_OBC(tcut, mLow) .* G_OBC(collect(tcut:T), mLow)
    elseif ens.bc == "pbc"
        # println("G_PBC(tcut, mUp): ", G_PBC(tcut, mUp))
        upperBound[tcut:end] = corr[tcut] ./G_PBC(tcut, mUp) .* G_PBC(collect(tcut:T), mUp)
        lowerBound[tcut:end] = corr[tcut] ./G_PBC(tcut, mLow) .* G_PBC(collect(tcut:T), mLow)
    else
        error("Boundary conditions for ensemble $(ens.id) are not recognised. Supported values are \"obc\" or \"pbc\" ")
    end
    return lowerBound, upperBound
end
function corrBounding(corr::Vector{uwreal}, mUp33::uwreal, mUp88, mLow::uwreal, tcut::Int64, ens::EnsInfo)
    T = HVPobs.Data.get_T(ens.id)

    upperBound = corr[:]
    lowerBound = corr[:]

    @. G_PBC(x0, p) = exp(-p*(x0-tcut)) + exp(-p*(T-x0+tcut))
    @. G_OBC(x0, p) = exp(-p*(x0-tcut))

    if ens.bc == "obc"
        upperBound[tcut:end] = corr[tcut] .* (G_OBC(collect(tcut:T), mUp33) - G_OBC(collect(tcut:T), mUp88))
        lowerBound[tcut:end] = corr[tcut] .* G_OBC(collect(tcut:T), mLow)
    elseif ens.bc == "pbc"
        upperBound[tcut:end] = corr[tcut] .* (G_PBC(collect(tcut:T), mUp33) - G_PBC(collect(tcut:T), mUp88))
        lowerBound[tcut:end] = corr[tcut] .* G_PBC(collect(tcut:T), mLow)
    else
        error("Boundary conditions for ensemble $(ens.id) are not recognised. Supported values are \"obc\" or \"pbc\" ")
    end
    return lowerBound, upperBound
end

function boundingMethod(corr::Vector{uwreal}, ens::EnsInfo, krnl::Vector{Union{Tu, Tf}}, sector::String; pl::Bool=true, wind::Union{Window,Nothing}=nothing, path_pl::Union{Nothing,String}=nothing, qval::String="", data::Bool=false ) where {Tu, Tf}

    T = HVPobs.Data.get_T(ens.id)
    thalf = Int64(T/2)
    tcut_arr = collect(8: thalf - ceil(Int64, 0.2 / HVPobs.Data.get_a(ens.beta) ))
    
    mpi = m_ens[ens.id]["m_pi"] 
    mrho = m_ens[ens.id]["m_rho"]
    E2pi = 2*sqrt(mpi^2 + (2*pi/ens.L)^2)
    E3pi = 2*sqrt(mpi^2 + (2*pi/ens.L)^2) + sqrt(mpi^2 + 2(2*pi/ens.L)^2)

    meffdata = meff(corr)

    store_tmrUp = []
    store_tmrLow = []
    if sector == "33"
        mUp = mrho < E2pi ? mrho : E2pi
        ll = mrho < E2pi ? L"$M_\rho$" : L"$E_{\pi\pi}$"
        yll = L"$\mathit{\bar\Pi}^{3,3}(-Q^2)$"
    elseif sector == "88"
        mUp = mrho
        ll =  L"$M_\omega$"
        yll = L"$\mathit{\bar\Pi}^{8,8}(-Q^2)$"
    elseif sector == "3388"
        mUp33 = mrho < E2pi ? mrho : E2pi 
        mUp88 = mrho
        println("mUp33: ", value(mUp33), " mUp88: ", value(mUp88) )
        ll =  L"$M_\rho$"
        yll = L"$\Delta_{ls}(-Q^2)$"
    elseif sector == "08"
        mUp = mrho 
        ll =  L"$M_\rho$"
        yll = L"$\mathit{\bar\Pi}^{0,8}(-Q^2) $"
    else
        error("Sector $(sector) not recognised. Supported values are \"33\", \"88\", \"08\" ")
    end

    for tcut in tcut_arr
        if  HVPobs.Data.get_a(ens.beta)*tcut < 1.8
            mLow =  meffdata[tcut]
        else
            tcutaux = ceil(Int64, 1.8 / HVPobs.Data.get_a(ens.beta) )
            mLow = meffdata[tcutaux]
        end
        #println("M_eff low bound at t=$(tcut): ", mLow)

        if sector == "3388"
            lowerBound, upperBound = corrBounding(corr, mUp33, mUp88, mLow, tcut, ens)
        else
            lowerBound, upperBound = corrBounding(corr, mUp, mLow, tcut, ens)
        end
    
        tmrUp_tmp  = upperBound .* krnl 
        tmrLow_tmp = lowerBound .* krnl

        if !isnothing(wind)
            t_fm = Float64.(collect(0:T-1)) .* HVPobs.Data.get_a(ens.beta)
            tmrUp_tmp .*= wind(t_fm)
            tmrLow_tmp .*= wind(t_fm)
        end
        push!(store_tmrUp, sum(tmrUp_tmp[1:thalf]))
        push!(store_tmrLow, sum(tmrLow_tmp[1:thalf]))
    end
    uwerr.(store_tmrLow)
    uwerr.(store_tmrUp)

    tmin = 0#tcut_arr[1]
    for k in eachindex(store_tmrLow)
        difff = store_tmrUp[k] - store_tmrLow[k]
        # println(abs(value(difff)))
        # println(0.5*max(err(store_tmrUp[k]), err(store_tmrLow[k])))
        if abs(value(difff)) < 0.5 * max(err(store_tmrUp[k]), err(store_tmrLow[k]))
            tmin += k
            break
        end
    end
    res = 0.0
    tmax = 0 
    try
        tmax = tmin + ceil(Int64, 0.8 / HVPobs.Data.get_a(ens.beta)) # 0.8fm is the length of the averaged interval
        res = plat_av(0.5*(store_tmrLow[tmin:tmax] .+ store_tmrUp[tmin:tmax])) 
    catch
        tmax = tmin + ceil(Int64, 0.4 / HVPobs.Data.get_a(ens.beta)) # 0.4fm is the length of the averaged interval
        res = plat_av(0.5*(store_tmrLow[tmin:tmax] .+ store_tmrUp[tmin:tmax])) 
    end

    if pl 
        tmin += tcut_arr[1]
        tmax += tcut_arr[1]
        uwerr(res)

        fill_between(tcut_arr, value(res)-err(res), value(res)+err(res), color="lightskyblue", alpha=0.8)

        errorbar(tcut_arr, value.(store_tmrLow), err.(store_tmrLow), fmt="d", capsize=2, ms=4, color="navy", label=L"$M_{\mathrm{eff}}$")
        errorbar(tcut_arr, value.(store_tmrUp), err.(store_tmrUp), fmt="d", capsize=2, ms=4, color="darkred", label=ll )
 
        axvline(tmin, ls=":", color="gray")
        axvline(tmax, ls=":", color="gray")
        ylim(value(res) - 20* err(res), value(res) + 20* err(res))
        xlim(tcut_arr[1], tcut_arr[end])
        ylabel(yll)
        xlabel(L"$t_{\mathrm{cut}}/a$")
        legend()
        tight_layout()
        display(gcf())
        if !isnothing(path_pl)
            tt = joinpath(path_pl, "boundMethod_$(ens.id)_$(sector)_Q$(qval).pdf")
            savefig(tt)
        end
        close("all")
    end
    if !data
        return res
    else
        return res, store_tmrLow, store_tmrUp, tcut_arr
    end
end

function reconstr_corr(
    ens::EnsInfo,
    E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};
    nmax::Union{Nothing,Int64}=nothing,
    impr_set::String="1",
    IMPR::Bool=true,
    std::Bool=true,
    RENORM::Bool=true,
    FreeEN::Bool=false,
    total::Bool=false,
)
    length(Z) != length(Zimpr) ? error("Length of Z and Zimpr are not the same") : nothing
    # nmax chooses the tower of states to be added
    # if nmax is not defined, then the max number of available states are used
    # FreeEN ensures that the last energy level is not used when saturating the corr; usefull for impr BM
    if isnothing(nmax)
        FreeEN ? nmax = minimum([length(E)-1,length(Z)]) : nmax = minimum([length(E),length(Z)])
    else
        if FreeEN
            nmax > length(E)-1 ? error("Input E data not large enough to tower $nmax (-1) states") : nothing
        else
            nmax > length(E) ? error("Input E data not large enough to tower $nmax states") : nothing
        end
        nmax > length(Z) ? error("Input Z data not large enough to tower $nmax states") : nothing
    end

    # t = collect(1:Int64(HVPobs.Data.get_T(ens.id))/2+1)
    t = collect(1:Int64(HVPobs.Data.get_T(ens.id)))
    corrVec_PiPi = [corr_n(1,t.-1,E,Z.^2 ,ens.L)]
    corrVec_Impr = [corr_n(1,t.-1,E,Zimpr.*Z,ens.L)]
    for n=2:nmax
        push!(corrVec_PiPi,corrVec_PiPi[n-1] .+ corr_n(n,t.-1,E,Z.^2 ,ens.L))
        push!(corrVec_Impr,corrVec_Impr[n-1] .+ corr_n(n,t.-1,E,Zimpr.*Z,ens.L))
    end

    # corrVec = deepcopy(corrVec_PiPi)
    if IMPR
        if impr_set == "1"
            cv_l = cv_loc(ens.beta)
        elseif impr_set == "1old"
            cv_l = cv_loc_old(ens.beta)
        elseif impr_set == "2"
            cv_l = cv_loc_set2(ens.beta)
        else
            error("Impr Set $impr_set not recoginsed, please choose between '1', '1old' and '2'.")
        end
        # [1:Int64(HVPobs.Data.get_T(ens.id)/2+1)]
        # [improve_corr_vkvk!(corrVec[n], corrVec_JPi[n], 2*cv_l, std=std, treelevel=true) for n in collect(1:nmax)]
        #corrVec = [improve_corr_vkvk!(corrVec_PiPi[n], corrVec_Impr[n], 2*cv_l, std=std) for n in 1:nmax]
        corrVec = [corrVec_PiPi[n] .+ (2*cv_l).*corrVec_Impr[n] for n=1:nmax]
    else
        corrVec = corrVec_PiPi
    end

    if RENORM
        Z3 = get_Z3(ens, impr_set=impr_set)
        [renormalize!(corr, Z3^2) for corr in corrVec]
        if !total
            [renormalize!(corr, Z3^2) for corr in corrVec_PiPi]
            [renormalize!(corr, Z3^2) for corr in corrVec_Impr]
        end
    end

    !total ? (return corrVec) : (return corrVec, corrVec_PiPi, corrVec_Impr)
end

reconstr_corr(ensid::String,E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};nmax::Union{Nothing,Int64}=nothing,impr_set::String="1",IMPR::Bool=true,RENORM::Bool=true,FreeEN::Bool=false,total::Bool=false,) = reconstr_corr(EnsInfo(ensid),E::Vector{uwreal},Z::Vector{uwreal},Zimpr::Vector{uwreal};nmax=nmax,impr_set=impr_set,IMPR=IMPR,RENORM=RENORM,FreeEN=FreeEN,total=total,)

function corr_n(n::Int64,t::Union{Vector{Int64},Vector{Float64}},E::Vector{uwreal},Z::Vector{uwreal},L::Int64)
    return (1/2) * Z[n] / (L^3) .* exp.(-E[n] .* t) # the 1/2 factor converts light-light corr to the iso-vector corr
end

function findfirst_uninterrupted(bools)
    i = 0
    stop = false
    while !stop
        arg = findfirst(bools[i+1:end])
        i = !isnothing(arg) ? (i+arg) : nothing
        if isnothing(i) || isnothing(findfirst(.!bools[i+1:end]))  # stops if there are unterrupted trues or if it reached the end
            stop = true
        end
    end
    return i
end

function get_spectr_data(path::String,ens::EnsInfo)
    ph5 = joinpath(path,"$(ens.id)-pSq0-T1up.h5")
    file = h5open(ph5, "r")

    # file["summary"]

    # Read the energies
    EnKeys = HDF5.keys(file["pSq0-T1up/spectrum"])
    E = uwreal[]
    for n in EnKeys
        En = read(file["pSq0-T1up/spectrum/$n"])
        push!(E,uwreal(jackknife_err(En),n))
    end

    # Read the overlaps
    ZKeys = HDF5.keys(file["pSq0-T1up/overlaps"])
    Z = uwreal[]
    Z_impr = uwreal[]
    for p in ZKeys
        Zp = read(file["pSq0-T1up/overlaps/$p"])
        Zp_impr = read(file["pSq0-T1up/overlaps_imp/$p"])
        push!(Z,uwreal(jackknife_err(Zp),p))
        push!(Z_impr,uwreal(jackknife_err(Zp_impr),p))
    end

    return E, Z, Z_impr
end

function jackknife_resampling(data)
    N = length(data)
    means = zeros(N)
    
    for i in 1:N
        sample = deleteat!(copy(data), i)  # Leave-one-out sample
        means[i] = mean(sample)
    end
    
    mean_jack = mean(means)
    std_jack = sqrt((N - 1) * mean((means .- mean_jack) .^ 2))
    
    return [mean_jack, std_jack]
end

function bootstrap_resampling(data, B=1000)
    N = length(data)
    means = zeros(B)
    
    for i in 1:B
        resample = data[rand(1:N, N)]  # Resample with replacement
        means[i] = mean(resample)
    end
    
    mean_boot = mean(means)
    std_boot = std(means)  # Standard deviation of bootstrap samples
    
    return [mean_boot, std_boot]
end

function jackknife_err(jakknife_sample,∆dof=0)
    N    = length(jakknife_sample)
    Mean = mean(jakknife_sample)
    Err  = sqrt(N*(N-1)/(N-∆dof)) * sqrt(mean((jakknife_sample.-Mean).^2))
    return [Mean, Err]
end