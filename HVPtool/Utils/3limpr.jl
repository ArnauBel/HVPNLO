
using QuadGK

function treelevel_continuum_correlator(t)
    return 1. / (2* pi^2 * t^3)
end


function compute_HVPtl0(diag::String,wind::String,Qlist::Union{Vector{Float64},Vector{Int64}},path_coef::String)
    exp_diag = diag == "LO" ? 2 : 3
    Tildef = Dict("LO" => Tildef2, "NLOa" => Tildef4a, "NLOb" => Tildef4b, "NLOa&b" => (x,path) -> Tildef4a(x,path) + Tildef4b(x,path))

    f = nothing  # initialize f

    if wind == "SDsub"
        tl_cont = Float64[]
        for Q in Qlist
            f = x0 -> treelevel_continuum_correlator(x0) .* (Window("SD")(x0) * hbarc^2 * Tildef[diag](massmu/hbarc * x0,path_coef) - (Window("SD")(0) * (16/(Q/hbarc)^2)^2 * π^2 * (massmu/hbarc)^2 * C4[diag](massmu/hbarc * x0) * sin((Q/hbarc/4) * x0)^4))
            res, _ = quadgk(f, 0, 5,rtol=1e-3)
            amu = (alpha/pi)^exp_diag * res * 1e10 / 2
            push!(tl_cont, amu)
        end
    elseif wind == "SD"
        f = x0 -> treelevel_continuum_correlator(x0) .* (Window("SD")(x0) * hbarc^2 * Tildef[diag](massmu/hbarc * x0,path_coef))
        res, _ = quadgk(f, 0, 5,rtol=1e-3)
        tl_cont = (alpha/pi)^exp_diag * res * 1e10 / 2
    else
        error("3l improvement cannot be applied to wind $wind")
    end
    return tl_cont
end
