module HVPtool

include("uwConst.jl")

export artificial_err
export sqrtt0_ph_Regensburg, sqrtt0_ph_Madrid, sqrtt0_ph_Bruno, sqrtt0_ph_CLS
export fPi_ph_PDGFLAG
export sqrtt0_ph, fPi_ph, fK_ph, MD_ph, mpi_ph, mK_ph
export phi2_ph, phi4_ph, y_ph, z_ph
export b_values, t0sym, fPiph, asym, a2_rescaling
export Zvc_l, E0_ens, meson_ens # m_ens

include("Const/Const.jl")

export hbarc, masse, massmu, alpha, GammaEuler, Qlist, C4, kcd_in, p0_smallpbc_dict, wpmm, NLOab_zero

include("Corr/Corr.jl")

export corr33,  corr88_conn, corr08_conn, corrC_conn

include("TMR/TMR.jl")

using .TMR
export f2, Tildef2_num, Tildef2, Tildef2Inner
export f4a, Tildef4a_num, Tildef4a, Tildef4aInner
export f4b, Tildef4b_num, Tildef4b, Tildef4bInner
export Tildef4c, Tildef4cInner

include("isovModel.jl")

export a2, a2loga, a3, a4
export a2phi2, a2phi4, a3phi2, phi2, phi2sqr, phi2inv, phi2log, logphi2, phi4, phi4sqr, phi4inv, phi4log, logphi4
export a2y, a2z, a3y, y, ysqr, yinv, ylog, logy, z, zsqr, zinv, zlog, logz
export deltaphi
export call_models

include("ModAver/ModAver.jl")

using. ModAver
export FitCat
export get_w_from_fitres, get_w_from_fitcat
export model_average, renorm_w!

include("WriteRead/WriteRead.jl")

using .WriteRead
export paste_str, func_str, create_path
export BDIOread_t0, BDIOread_TMR, BDIOread_corr, BDIOread_FVCcorr
export BDIOread_HVPens, BDIOread_FVCens
export BDIOread_res, BDIOread_XYdata, JDL2read_FitRes, JDL2read_ModelInfo
export BDIOread_MA, BDIOread_MAtot
export BDIOread_mPP, BDIOread_fPS, BDIOread_mDs_kappaC, BDIOread_mDs
export TXTread_FVCcorr_GS, JDL2read_FVC_ChPT, TXTread_bQ

include("Utils/Utils.jl")

using .Utils
export ensCheck, Window2D
export meff_MA
export treelevel_continuum_correlator, compute_HVPtl0
export apply_syst_HVP, apply_syst_HVP!, apply_syst_FVC, apply_syst_FVC!, HVP_VolCorrect, HVP_VolCorrect!, HVP_3limpr!
export Eeff, corr_bound, bounding_method
export findfirst_uninterrupted, get_spectr_data, corr_n, reconstr_corr
export jackknife_resampling, bootstrap_resampling, jackknife_err
export set_fluc_to_zero!, add_t0_err!, get_t0err

include("UwOper/UwOper.jl")

end