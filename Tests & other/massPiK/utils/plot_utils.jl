function plot_tmr(tmr::Vector{uwreal}, ens::EnsInfo, q2::Float64; p::Union{String,Nothing}=nothing, info::String="_")
    tvals = collect(0:Int64(length(tmr))-1)
    yy = tmr ./ sqrt(8*t0(ens.beta))
    uwerr.(yy)
    errorbar(tvals, abs.(value.(yy)), err.(yy), fmt="s", mfc="none", color="royalblue", label=L"$Q^2= $"*"$(q2) "*L"$\mathrm{GeV}^2$" )
    xlabel(L"$t/a$")
    ylabel(L"$G(t) \cdot K(t, Q^2)$")
    legend()
    tight_layout()
    display(gcf())
    if !isnothing(p)
        tt = string("tmr_$(ens.id)_Q", Int64(q2), "_", info, ".pdf")
        savefig(joinpath(p, ens.id, tt))
    end
    close("all")

end

function plot_tmr_multiq(tmr::Vector{Vector{uwreal}}, ens::EnsInfo, q2::Vector{Float64}; p::Union{String,Nothing}=nothing, info::String="_")
    tvals = collect(0:Int64(length(tmr[1]))-1)
    cc=["forestgreen", "royalblue", "red"]
    for k in eachindex(tmr)
        yy = tmr[k] ./ sqrt.(8 .*t0.(ens.beta)); uwerr.(yy)
        errorbar(tvals, abs.(value.(yy)), err.(yy), fmt="o", color=cc[k], label=L"$Q^2=\ $"*"$(q2[k]) "*L"$\mathrm{GeV}^2$" )
        fill_between(tvals, abs.(value.(yy))- err.(yy), abs.(value.(yy))+ err.(yy), color=cc[k], alpha=0.2 )
    end
    xlabel(L"$t/a$")
    ylabel(L"$G(t) \cdot K(t, Q^2) / \sqrt{8t_0}$")
    legend()
    twiny()
    xx_fm =  tvals #.* value.(t0sqrt_ph ./ sqrt.(t0.(ens.beta)))
    fm_time = round.((xx_fm[1:6:end] ) .* value.(t0sqrt_ph ./ sqrt.(t0.(ens.beta))) , digits=1)
    xticks(xx_fm[1:6:end], fm_time)
    xlabel(L"$t\ [\mathrm{fm}]$", size=14)
    tight_layout()
    display(gcf())
    if !isnothing(p)
        tt = string("tmr$(info)_$(ens.id).pdf")
        savefig(joinpath(p, ens.id, tt))
    end
    close("all")
end


@doc raw"""
    plot_cl_all_set(fc_ll_s1, fc_ll_s2, fc_lc_s1, fc_lc_s2; ylab::LaTeXString=L"$\bar{\Pi}^{88,\mathrm{sub}}(-Q^2)$", nmom::Int64=3, path_plot::Union{String,Nothing}=nothing)

This function takes as input four Vector{Vector{FitCat}} (for each discretization and set of improvement coefficients) and plots
the continuum limit extrapolations for all the sets simultaneously. The opacity of each line is associated to the 
weights given by the model average. One plot for each value of the momenta is generated.
"""
function plot_cl_all_set(fc_ll_s1, fc_ll_s2, fc_lc_s1, fc_lc_s2; ylab::LaTeXString=L"$\bar{\Pi}^{88,\mathrm{sub}}(-Q^2)$", nmom::Int64=3, path_plot::Union{String,Nothing}=nothing, f_tot_isov=f_tot_isov)
    
    for q in [6,8,12]#1:nmom
        fig = figure(figsize=(10,7.))
        for k_cat in eachindex(fc_lc_s1[q])

            #ww_ll_s1 = get_w_from_fitcat(fc_ll_s1[q], norm=true)
            #ww_ll_s2 = get_w_from_fitcat(fc_ll_s2[q], norm=true)
            ww_lc_s1 = get_w_from_fitcat(fc_lc_s1[q], norm=true)
            ww_lc_s2 = get_w_from_fitcat(fc_lc_s2[q], norm=true)

            # xarr = [Float64.(range(1e-8, 0.05, length=100)) fill(phi2_ph, 100) fill(phi4_ph, 100)]
            xarr = [Float64.(range(1e-8, 0.008, length=100)) fill(phi2_ph, 100) fill(phi4_ph, 100)]
            count=0
            model_tot = vcat([f_tot_isov for _ in 1:length(fc_lc_s1[q])]...)
            for (k_mod, model) in enumerate(model_tot)
                k_mod_fit = mod(k_mod,length(f_tot_isov))
                k_mod_fit == 0.0 ? k_mod_fit=length(f_tot_isov) : k_mod_fit
                # if q in [1, 2, 3, 4, 5, 6 ] &&  "a3" ∉ label_tot_dltiso[k_modaux]
                    # println("model $(label_tot_dltiso[k_modaux]) skipped")
                    # count +=1
                    # k_mod = count
                # elseif q in [1, 2, 3, 4, 5, 6 ] &&  "a3" ∈ label_tot_dltiso[k_modaux]
                    # continue
                # else
                    # k_mod = k_modaux
                # end
            
                # fit_param_ll_s1 = fc_ll_s1[q][k_cat].fit[k_mod_fit].param
                # fit_param_ll_s2 = fc_ll_s2[q][k_cat].fit[k_mod_fit].param
                fit_param_lc_s1 = fc_lc_s1[q][k_cat].fit[k_mod_fit].param
                fit_param_lc_s2 = fc_lc_s2[q][k_cat].fit[k_mod_fit].param
    
                # yy_ll_s1  = model(xarr, fit_param_ll_s1)
                # yy_ll_s2  = model(xarr, fit_param_ll_s2)
                yy_lc_s1  = model(xarr, fit_param_lc_s1)
                yy_lc_s2  = model(xarr, fit_param_lc_s2)
                
                # plot(xarr[:,1], value.(yy_ll_s1),alpha=ww_ll_s1[k_mod], color="#009E73")
                # plot(xarr[:,1], value.(yy_ll_s2),alpha=ww_ll_s2[k_mod], color="#0072B2")
                plot(xarr[:,1], value.(yy_lc_s1),alpha=ww_lc_s1[k_mod], color="#D55E00")
                plot(xarr[:,1], value.(yy_lc_s2),alpha=ww_lc_s2[k_mod], color="#E69F00")
    
            end
        end
        axvline.(unique(value.(fc_lc_s1[q][1].xdata[:,1])), linewidth=1.5, alpha=0.5, color="gray", ls=":")
        xlim(0.0,0.008)
        L2D = PyPlot.matplotlib.lines.Line2D
        custom_lines = [#L2D([0], [0], color="#009E73", lw=2),
                #L2D([0], [0], color="#0072B2", lw=2),
                L2D([0], [0], color="#D55E00", lw=2),
                L2D([0], [0], color="#E69F00", lw=2)]
        # legend( custom_lines,  [L"$\mathrm{LL,\ set\ 1}$", L"$\mathrm{LL,\ set\ 2}$", L"$\mathrm{LC,\ set\ 1}$", L"$\mathrm{LC,\ set\ 2}$"], ncol=2)
        legend( custom_lines,  [ L"$\mathrm{LC,\ set\ 1}$", L"$\mathrm{LC,\ set\ 2}$"], ncol=2)

        # xlabel(L"$a^2/(8t_0)$")
        xlabel(L"$a^2 \ \mathrm{[fm]}$")
        ylabel(ylab)
        tight_layout()
        display(fig)
        if !isnothing(path_plot)
            fname = string("cont_lim_tot_pi88_q", q, ".pdf")
            savefig(joinpath(path_plot, "continuum-lim", fname))
        end
        close("all")
    end
end

function plot_chiral_best_fit(fc::Vector{Vector{FitCat}}; nmom::Int64=3, nfit::Int64=0, ylab::LaTeXString=L"$\bar{\Pi}^{88,\mathrm{sub}}(-Q^2)$", tt::Union{Nothing,Vector{String}}=nothing, path_plot::Union{String,Nothing}=nothing, f_tot_isov=f_tot_isov)

    for q in [6,8,12]#1:nmom
        fig = figure(figsize=(10,7.))
        # println("using the first three momenta")
        println("\n- Momentum: $(q)")
        fccat = vcat(fc[q]...)
        w_tot  = get_w_from_fitcat(fccat)
        wmax, wmax_idx = findmax(w_tot)

        if nfit == 0 
            model_idx = mod(wmax_idx, length(f_tot_isov))
            if model_idx == 0 
                model_idx =length(f_tot_isov)
            end
            cat_idx = Int((wmax_idx - model_idx ) / length(f_tot_isov))+1
        else
            model_idx = nfit
            # cat_idx = 1
            cat_idx = Int((wmax_idx - model_idx ) / length(f_tot_isov))+1
        end
        println("   - highest weigth: ", wmax)
        println("   - idx max weigth: ", wmax_idx)
        println("   - idx best model: ", model_idx)
        println("   - idx best cat  : ", cat_idx)
        println("   - best χ2/χ2exp : ", fccat[cat_idx].fit[model_idx].chi2 / fccat[cat_idx].fit[model_idx].chi2exp )

        xdata = fccat[cat_idx].xdata; uwerr.(xdata)
        ydata = fccat[cat_idx].ydata; uwerr.(xdata)
        fit_par = fccat[cat_idx].fit[model_idx].param
        model = f_tot_isov[model_idx]

        xdata_proj = [xdata[:,1] xdata[:,2] fill(phi4_ph, length(xdata[:,1]))]
        ydata_proj = ydata - model(xdata, fit_par) + model(xdata_proj, fit_par) ; uwerr.(ydata_proj)

        if maximum(value.(xdata[:,1])) >  0.006#0.04
            betatot =  [3.4, 3.46, 3.55, 3.7, 3.85] 
            # color = ["tomato", "forestgreen", "violet", "orange", "navy"]
            color = ["#E69F00", "#56B4E9",  "#009E73", "#CC79A7", "#D55E00"] 
            fmttot = ["d", "s", "^", "h", "8"]
        else
            betatot =  [3.46, 3.55, 3.7, 3.85]
            # color = ["forestgreen", "violet", "orange", "navy"]
            color = [ "#56B4E9",  "#009E73", "#CC79A7", "#D55E00"] 
            fmttot = ["s", "^", "h", "8"]
        end
        # plot data points
        # for (k, b) in enumerate(betatot)
        for (k,b) in enumerate(sort(unique(value.(xdata[:,1])), rev=true))    
            # if b == 3.85
                # println("beta 3.85 removed!")
                # continue
            # end
            # n_ = findall(x->x.beta == b, ensinfo)
            n_ = findall(x->x == b, value.(xdata[:,1]))
            a2_aux  = mean(value.(xdata[:,1][n_]))
            errorbar(value.(xdata[n_,2]), value.(ydata_proj[n_]), xerr=err.(xdata[n_,2]), yerr=err.(ydata_proj[n_]), fmt=fmttot[k], label=string(L"$\beta = \ $", betatot[k]), color=color[k], capsize=2, ms=10, mfc="none" )
            # dashed lines
            xxx_tot = [fill(a2_aux, 100) Float64.(range(maximum(value.(xdata[:,2][n_])), 0.8, length=100)) fill(phi4_ph, 100)]
            xxx_chir = [fill(a2_aux, 100) Float64.(range(0.001, maximum(value.(xdata[:,2][n_])), length=100)) fill(phi4_ph, 100)]
            yy_chir = model(xxx_chir, fit_par)
            yy_tot = model(xxx_tot, fit_par)
            replace!(value.(yy_chir), Inf=>0.0)
            replace!(value.(yy_tot), Inf=>0.0)
            uwerr.(yy_chir)
            uwerr.(yy_tot)
            fill_between(xxx_chir[:,2], value.(yy_chir)-err.(yy_chir), value.(yy_chir)+err.(yy_chir), color=color[k], alpha=0.5)
            plot(xxx_tot[:,2], value.(yy_tot), ls="--", color=color[k], lw=0.9)
        end

        # cont lim band
        xarr = [fill(1e-8, 100) Float64.(range(0.01, 0.8, length=100)) fill(value(phi4_ph),100)]
        yarr = model(xarr, fit_par); uwerr.(yarr)
        fill_between(xarr[:,2], value.(yarr) .- err.(yarr), value.(yarr) .+ err.(yarr), alpha=0.2, color="gray")
        # phys res
        ph_res = model([0.0 value(phi2_ph) value(phi4_ph)], fit_par)[1]; uwerr(ph_res)
        # ph_res = fit_par[1]; uwerr(ph_res) #model([0.0 value(phi2_ph) value(phi4_ph)], fit_par)[1]; uwerr(ph_res)
        println(ph_res)
        errorbar(value(phi2_ph), value(ph_res), err(ph_res), fmt="o", capsize=2, color="red", ms=10, mfc="none")
        axvline(value(phi2_ph), ls="dashed", color="black", lw=0.2, alpha=0.7) 
        
        xlim(0.04, 0.8)
        #ylim(0.0097, 0.0105)
        legend(ncol=2)#, loc="upper right")
        xlabel(L"$\phi_2$")
        ylabel(ylab)
        if !isnothing(tt)
            title(join(tt, " "))
        end
        tight_layout()
        display(fig)
        if !isnothing(path_plot)
            fname = string("chiral_lim_", join(tt,"_"), "q_$(q).pdf") 
            savefig(joinpath(path_plot, "chiral-lim", fname ))
        end
        close("all")
    end
end

function plot_cl_best_fit(fc::Vector{Vector{FitCat}}; nmom=3, ylab::LaTeXString=L"$\bar{\Pi}^{88,\mathrm{sub}}(-Q^2)$", tt::Union{Nothing,Vector{String}}=nothing, path_plot::Union{String,Nothing}=nothing, f_tot_isov=f_tot_isov )
    
    for q in [6,8,12]#1:nmom
        fig = figure(figsize=(10,7.))
        println("\n- Momentum: $(q)")
        fccat = vcat(fc[q]...)
        w_tot  = get_w_from_fitcat(fccat)
        wmax, wmax_idx = findmax(w_tot)

        model_idx = mod(wmax_idx, length(f_tot_isov))
        if model_idx == 0 
            model_idx =length(f_tot_isov)
        end
        cat_idx = Int((wmax_idx - model_idx ) / length(f_tot_isov))+1

        println("   - highest weigth: ", wmax)
        println("   - idx max weigth: ", wmax_idx)
        println("   - idx best model: ", model_idx)
        println("   - idx best cat  : ", cat_idx)
        println("   - best χ2/χ2exp : ", fccat[cat_idx].fit[model_idx].chi2 / fccat[cat_idx].fit[model_idx].chi2exp )

        xdata = fccat[cat_idx].xdata; uwerr.(xdata)
        ydata = fccat[cat_idx].ydata; 
        fit_par = fccat[cat_idx].fit[model_idx].param
        model = f_tot_isov[model_idx]

        xdata_proj = [xdata[:,1] fill(phi2_ph, length(xdata[:,1])) fill(phi4_ph, length(xdata[:,1]))]
        ydata_proj = model(xdata_proj, fit_par); uwerr.(ydata_proj)

        errorbar(value.(xdata[:,1]), xerr=err.(xdata[:,1]), value.(ydata_proj), yerr=err.(ydata_proj), fmt="p", capsize=2, color="royalblue")

        # phys res
        ph_res = model([0.0 phi2_ph phi4_ph], fit_par)[1]; uwerr(ph_res)
        println(ph_res)
        errorbar(0.0, value(ph_res), err(ph_res), fmt="o", capsize=2, color="red", ms=10, mfc="none")

        # cont lim band
        # xarr = [Float64.(range(0.00, 0.047, length=100)) fill(phi2_ph, 100)  fill(phi4_ph,100)]
        xarr = [Float64.(range(0.00, 0.008, length=100)) fill(phi2_ph, 100)  fill(phi4_ph,100)]
        yarr = model(xarr, fit_par); uwerr.(yarr)
        fill_between(xarr[:,1], value.(yarr) .- err.(yarr), value.(yarr) .+ err.(yarr), alpha=0.2, color="royalblue")

        # xlabel(L"$a^2/(8t_0)$")
        xlabel(L"$a^2 \ [\mathrm{fm}^2]$")
        ylabel(ylab)
        if !isnothing(tt)
            title(join(tt, " "))
        end
        tight_layout()
        display(fig)
        if !isnothing(path_plot)
            fname = string("cont_lim_", join(tt,"_"), "q_$(q).pdf") 
            savefig(joinpath(path_plot, "continuum-lim", fname ))
        end
        close("all")
    end
end

function plot_mAve_summary(fc::Vector{Vector{FitCat}}; nmom=3, ylab::Union{Nothing, LaTeXString}=nothing, xlab::Union{Nothing, Vector{Vector{String}}}=nothing, charge_factor::Float64=1., models::Union{Vector{Function}, Nothing}=nothing, path_plot::Union{String,Nothing}=nothing, tt::Union{Nothing,Vector{String}}=nothing)
    for q in [6,8,12]#1:nmom
        println("\n- Momentum: $(q)")
        fccat = vcat(fc[q]...)
        w_tot  = get_w_from_fitcat(fccat)
        wmax, wmax_idx = findmax(w_tot)
        if isnothing(models)
            all_res = vcat([getindex.(getfield.(getfield(fccat[k], :fit), :param), 1) for k in eachindex(fccat)]...)
        else
            all_res = Vector{uwreal}()
            for (k, c) in enumerate(fccat)
                for (j,mod) in enumerate(models)
                    push!(all_res, mod([0.0 phi2_ph phi4_ph], c.fit[j].param)[1])
                end
            end
        end
        chi2exp = vcat([getfield.(getfield(fccat[k], :fit), :chi2exp) for k in eachindex(fccat)]...)
        chi2    = vcat([getfield.(getfield(fccat[k], :fit), :chi2) for k in eachindex(fccat)]...)
        pval    = vcat([getfield.(getfield(fccat[k], :fit), :pval) for k in eachindex(fccat)]...)
        chi2_over_chi2exp = chi2 ./ chi2exp
        
        final_res, syst = model_average(all_res, w_tot) .* charge_factor; uwerr(final_res)
        final_err = sqrt(err(final_res)^2 + syst^2)
        println("model ave: ", final_res )

        # sort weights and discard 5% tail
        idxW = sortperm(w_tot, rev=true)
        cumulative_w = cumsum(w_tot[idxW])
        idxcumw = findfirst(x -> x > 0.95, cumulative_w)
        idxW = sort(idxW[1:idxcumw])
        
        all_res = all_res[idxW] .* charge_factor; uwerr.(all_res)
        w_tot   = w_tot[idxW]
        chi2_over_chi2exp = chi2_over_chi2exp[idxW]
        pval = pval[idxW]


        fig = figure(figsize=(10,7.5))
        subplots_adjust(hspace=0.1) 

        subplot(411)  
        if !isnothing(tt)
            title(join(tt, " "))
        end
        ax1 = gca()
        setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels
        setp(ax1.get_xticklines(),visible=false) # Disable x tick lines
        x = collect(1:length(all_res))
        v = value.(all_res)
        e = err.(all_res)

        axhline(value(final_res), ls="-.", linewidth=1, color="red")
        axhline(value(final_res)-final_err, ls=":", linewidth=1, color="red")
        axhline(value(final_res)+final_err, ls=":", linewidth=1, color="red")

        errorbar(x, v, e, fmt="d", color="orange", capsize=2)
        if !isnothing(ylab)
            ylabel(ylab)
        end

        subplot(412)
        ax2 = gca()
        setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels
        bar(x, w_tot, alpha=0.4, color="royalblue", edgecolor="blue", linewidth=1.5)
        ylim(0,1)
        ylabel(L"$W$")

        subplot(413)
        ax3 = gca()
        setp(ax3.get_xticklabels(),visible=false) # Disable x tick labels
        pval = fill(1., length(w_tot))
        bar(x, pval, alpha=0.4, color="forestgreen", edgecolor="darkgreen", linewidth=1.5)
        ylim(0,1)
        ylabel(L"$p-values$")

        subplot(414)
        bar(x, chi2_over_chi2exp, alpha=0.4, color="tomato", edgecolor="darkred", linewidth=1.5)
        ylabel(L"$\chi^2 / \chi^2_{\mathrm{exp}}$")
        xticks(x, join.(xlab[idxW], " "), rotation=75, size=10)
        
        tight_layout()
        display(fig)
        if !isnothing(path_plot)
            fname = string("summaryAv_", join(tt,"_"), "q_$(q).pdf") 
            savefig(joinpath(path_plot, "modelAve", fname))
        end

        close("all")

    end
end