wpmm = Dict{String, Vector{Float64}}()
# wpmm["H101"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["H102r002"] = [5.0, -2.0, -1.0, -1.0]
# wpmm["H400"]     = [5.0, -1.5, -1.0, -1.0]
# wpmm["N202"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["N200"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["N203"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["N300"]     = [5.0, -1.5, -1.0, -1.0]
# wpmm["J303"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["J304"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["F300"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["J306"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["J307"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["J500"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["E300"]     = [5.0, -2.0, -1.0, -1.0]
wpmm["A654"]     = [5.0, -2.0, -1.0, -1.0]
# wpmm["E300"]     = [-1.0, 2.0, -1.0, -1.0]


function meff_MA(corr::Corr; 
    pl2state0::Vector{Float64}=[0.5,0.6], plconst0::Vector{Float64}=[0.7,1.0], plf::Float64=1.0, plstep::Int64=1, 
    mdof::Int64=4, 
    state_fit::Bool=true, 
    AIC::Bool=true, 
    returnfitMA::Bool=false, 
    plot::Bool=false, 
    pval::Bool=false,
    fitinfo::Bool=false,
    bc::String="obc"
    )
    
    obs = corr.obs
    T = HVPobs.Data.get_T(corr.id)

    if bc == "obc"
        m_obs = meff(obs[1:Int64(T/2+1)]) # first and last data points are lost (-2) and derivative is taken (-1); length(meff_) = length(obs) - 3
        # m_obs = meff(obs) # first and last data points are lost (-2) and derivative is taken (-1); length(meff_) = length(obs) - 3
        len = length(m_obs)

        @. const_model(x,p) = p[1] + 0*x

        plconst0_vec = collect(floor(Int64,plconst0[1]*len):plstep:ceil(Int64,min(plf*len-mdof,plconst0[2]*len)))
        pl_f = ceil(Int64,plf*len)

        if isempty(plconst0_vec)
            error("Not enough data to fit with d.o.f. ≥ $mdof. Decrease $mdof or increase the plateau search range. ")
        end

        wpm = corr.id in keys(wpmm) ? wpmm : nothing

        fitconst_vec = Vector{FitRes}()
        fitstate_vec = Vector{FitRes}()
        for p0 in plconst0_vec
            # plateau = [p0,pl_f]
            m_data = m_obs[p0:pl_f] 
            fit = fit_routine(const_model,collect(p0:pl_f).+1.5, m_data, 1, pval=pval, wpm=wpm, info=fitinfo, lineprint=fitinfo)
            push!(fitconst_vec,fit)
        end

        if state_fit
            @. state_model(x,p) = p[1] + p[2] * exp(- p[3] * x)
            # @. state_model(x,p) = (p[1] + (p[1] + p[3]) * p[2] * exp(- p[3] * x)) / (1 + p[2] * exp(- p[3] * x))

            plstate0_vec = collect(max(floor(Int64,pl2state0[1]*len),1):plstep:ceil(Int64,pl2state0[end]*len))

            p0erratics = []
            for p0 in plstate0_vec
                try
                    m_data = m_obs[p0:pl_f]
                    fit = fit_routine(state_model,collect(p0:pl_f).+1.5, m_data, 3, pval=pval, wpm=wpm, info=fitinfo, lineprint=fitinfo)
                    push!(fitstate_vec,fit)
                catch
                    push!(p0erratics,p0)
                    # filter!(x -> x != p0, plstate0_vec)
                end
            end

            plstate0_vec = filter(x -> !(x in p0erratics), plstate0_vec)

            p0_vec  = [plstate0_vec,plconst0_vec]
            fit_vec = [fitstate_vec,fitconst_vec]
        else
            p0_vec  = plconst0_vec
            fit_vec = fitconst_vec
        end
    elseif bc == "pbc"
        len = length(obs)/2+1
        obs = -obs[1:Int(len)]

        if state_fit
            error("Fit to ground+exited stated still not implemented for pbc")
        end

        @. cosh_model(x,p) = p[2] * (exp(-p[1]*x) + exp(-p[1]*(T-x))) +  p[4] * (exp(-p[3]*x) + exp(-p[3]*(T-x)))

        plconst0_vec = collect(floor(Int64,plconst0[1]*len):plstep:ceil(Int64,min(plf*len-mdof,plconst0[2]*len)))
        pl_f = ceil(Int64,plf*len)

        if isempty(plconst0_vec)
            error("Not enough data to fit with d.o.f. ≥ $mdof. Decrease $mdof or increase the plateau search range. ")
        end

        wpm = corr.id in keys(wpmm) ? wpmm : nothing

        fitconst_vec = Vector{FitRes}()
        fitstate_vec = Vector{FitRes}()
        for p0 in plconst0_vec
            # plateau = [p0,pl_f]
            data = obs[p0:pl_f] 
            fit = fit_routine(cosh_model,collect(p0:pl_f).-1, data, 4, pval=pval, wpm=wpm, info=fitinfo, lineprint=fitinfo)
            push!(fitconst_vec,fit)
        end

        p0_vec  = plconst0_vec
        fit_vec = fitconst_vec
    else
        error("bc = $bc not recognised")
    end

    w = get_w_from_fitres(vcat(fit_vec...), AIC=AIC)
    res_vec = [par[1] for par in getfield.(vcat(fit_vec...),:param)]
    meff_res, meff_sys = model_average(res_vec, w)
    meff_res = meff_res[1]

    res = meff_res + uwreal([0.0,meff_sys],"meff MA syst"); uwerr(res)

    if plot
        fig = figure(figsize=(16,12))
        # subplots_adjust(hspace=0.1)
        gs = fig.add_gridspec(4, 1, height_ratios=[4, 1, 1, 1])  # Adjust the height_ratios as needed

        if bc == "obc"
            ax1 = fig.add_subplot(gs[1, 1])
            if state_fit
                x0 = collect(max(floor(Int64,len*(pl2state0[1]))+1-2,1):len+1) .+ 0.5
                m_vec = m_obs[max(floor(Int64,len*(pl2state0[1]))-2,1):end]; uwerr.(m_vec)
            else
                x0 = collect(max(floor(Int64,len*(plconst0[1]))+1-2,1):len+1) .+ 0.5
                m_vec = m_obs[max(floor(Int64,len*(plconst0[1]))-2,1):end]; uwerr.(m_vec)
            end

            title("$(corr.id)")
            errorbar(x0, value.(m_vec), err.(m_vec), fmt="o", capsize=2, color="black")
            fill_between(x0, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")
            if state_fit
                maxw_arg = argmax(w[1:length(fit_vec[1])])
                par = getfield.(fit_vec[1],:param)[maxw_arg]
                x_ = collect(max(floor(Int64,len*(pl2state0[1]))+1-2,1):0.1:len+1) .+ 0.5
                y2st_fit = state_model(x_,par); uwerr.(y2st_fit)
                fill_between(x_, value.(y2st_fit)+err.(y2st_fit), value.(y2st_fit)-err.(y2st_fit), alpha=0.3, color="orange")

                axvline(x=plstate0_vec[1]+1+0.1, color="red", linestyle="--")
                axvline(x=plstate0_vec[end]+2-0.1, color="red", linestyle="--")
            end
            axvline(x=plconst0_vec[1]+1+0.1, color="cyan", linestyle="--")
            axvline(x=plconst0_vec[end]+2-0.1, color="cyan", linestyle="--")
            if plf != 1.0
                axvline(x=ceil(Int64,plf*len)+2, color="gray", linestyle=":")
            end
            axis("tight")
            ylabel(L"$m_{\rm{eff}}$")
            ylim(res.mean-3*res.err,maximum(value.(m_vec) .+ 1.5 .* err.(m_vec)))
            setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels
        elseif bc == "pbc"
            ax1 = fig.add_subplot(gs[1, 1])

            x0 = collect(max(floor(Int64,len*(plconst0[1])),1)-2:len) .- 1
            vec = obs[Int(x0[1])+1:end]; uwerr.(vec)


            title("$(corr.id)")
            errorbar(x0, value.(vec), err.(vec), fmt="o", capsize=2, color="black")
            maxw_arg = argmax(w[1:length(fit_vec)])
            par = getfield.(fit_vec,:param)[maxw_arg]

            x_ = collect(max(floor(Int64,len*(plconst0[1])),1)-0.5:0.1:len).-1
            ycosh_fit = cosh_model(x_,par); uwerr.(ycosh_fit)
            fill_between(x_, value.(ycosh_fit)+err.(ycosh_fit), value.(ycosh_fit)-err.(ycosh_fit), alpha=0.3, color="orange")

            axvline(x=plconst0_vec[1]-1.5, color="red", linestyle="--")
            axvline(x=plconst0_vec[end]-0.5, color="red", linestyle="--")

            axis("tight")
            ylabel(L"$C(t)$")

            plconst0_vec = Float64.(plconst0_vec) .- (1.5 + 1.0)
        end

        ax2 = fig.add_subplot(gs[2, 1])
        uwerr.(res_vec)
        if state_fit
            errorbar(plstate0_vec .+ 1.5, value.(res_vec[1:length(plstate0_vec)]), err.(res_vec[1:length(plstate0_vec)]), fmt="d", mfc="none", capsize=2, color="orange")
            errorbar(plconst0_vec .+ 1.5, value.(res_vec[end-length(plconst0_vec)+1:end]), err.(res_vec[end-length(plconst0_vec)+1:end]), fmt="d", mfc="none", capsize=2, color="blue")
        else
            errorbar(plconst0_vec .+ 1.5, value.(res_vec), err.(res_vec), fmt="d", mfc="none", color="blue")
        end
        fill_between(x0, value.(res)+err.(res), value.(res)-err.(res), alpha=0.4, color="green")
        # ylabel(L"m_{D_s}")
        ylabel(L"m_0")
        setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels

        ax3 = fig.add_subplot(gs[3, 1])
        fill_between(x0, maximum(w)/2, maximum(w)/2, alpha=0.0, color="white")
        if state_fit
            PyPlot.plot(plstate0_vec .+ 1.5, w[1:length(plstate0_vec)], linestyle="none", marker="o", mfc="none", color="orange")
            PyPlot.plot(plconst0_vec .+ 1.5, w[end-length(plconst0_vec)+1:end], linestyle="none", marker="o", mfc="none", color="blue")
        else
            PyPlot.plot(plconst0_vec .+ 1.5, w, linestyle="none", marker="o", mfc="none", color="blue")
        end
        ylabel(latexstring("\\rm{w}\\ [\\rm{TIC}]"))

        if pval
            setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels
    
            ax4 = fig.add_subplot(gs[4, 1])
            if state_fit
                pval_vec = [getfield.(fit_vec_,:pval) for fit_vec_ in fit_vec]
                fill_between(x0, maximum(vcat(pval_vec...))/2, maximum(vcat(pval_vec...))/2, alpha=0.0, color="white")
                PyPlot.plot(plstate0_vec .+ 1.5, pval_vec[1], linestyle="none", marker="o", mfc="none", color="orange")
                PyPlot.plot(plconst0_vec .+ 1.5, pval_vec[2], linestyle="none", marker="o", mfc="none", color="blue")
            else
                pval_vec = getfield.(fit_vec,:pval)
                fill_between(x0, maximum(vcat(pval_vec...))/2, maximum(vcat(pval_vec...))/2, alpha=0.0, color="white")
                PyPlot.plot(plconst0_vec .+ 1.5, pval_vec, linestyle="none", marker="o", mfc="none", color="blue")
            end
            ylabel(L"$\rm{p-values}$")
        end
        xlabel(L"$t/a$")


        tight_layout()
        display(fig)
        close("all")

    end

    if returnfitMA
        return [meff_res, meff_sys], [p0_vec, fit_vec, w]
    else
        return [meff_res, meff_sys]
    end
end