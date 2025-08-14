## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> 2206.06582v2: Intermediate window paper (Tables VIII and IX) <==========================##


using ALPHAio, ADerrors
import ADerrors: err

using BDIO

using Printf

julia_script_directory = @__DIR__

# ensList = ["H101", "B450", "N202", "N300", "H102", "N101", "C101", "N203", "N200", "D200"]
ensList = ["H101","B450","N202","N300","H102","N101","C101","S400","N203","N200","D200","N302","E250","J303","E300","J500","A654","N451","D452","J501"]
ensInfo = EnsInfo.(ensList)

# Iso-vector contributions (Table VIII)
IsoVecRes = Dict{String, Dict{String, Dict{String, uwreal}}}()

IsoVecRes["H101"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([172.10, 0.39], "result"),
        "lc" => uwreal([173.35, 0.39], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([150.16, 0.39], "result"),
        "lc" => uwreal([155.36, 0.39], "result")
    )
)

IsoVecRes["H102"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([178.54, 0.52], "result"),
        "lc" => uwreal([179.75, 0.52], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([157.27, 0.53], "result"),
        "lc" => uwreal([162.26, 0.53], "result")
    )
)

IsoVecRes["N101"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([186.31, 0.43], "result"),
        "lc" => uwreal([187.56, 0.42], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([165.61, 0.44], "result"),
        "lc" => uwreal([170.48, 0.43], "result")
    )
)

IsoVecRes["C101"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([192.19, 0.41], "result"),
        "lc" => uwreal([193.40, 0.41], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([172.25, 0.43], "result"),
        "lc" => uwreal([176.94, 0.42], "result")
    )
)

IsoVecRes["B450"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([168.12, 0.38], "result"),
        "lc" => uwreal([168.82, 0.38], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([152.53, 0.38], "result"),
        "lc" => uwreal([155.68, 0.38], "result")
    )
)

IsoVecRes["N202"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([168.14, 0.68], "result"),
        "lc" => uwreal([168.45, 0.69], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([158.36, 0.67], "result"),
        "lc" => uwreal([159.92, 0.68], "result")
    )
)

IsoVecRes["N203"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([173.75, 0.43], "result"),
        "lc" => uwreal([174.11, 0.43], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([164.22, 0.43], "result"),
        "lc" => uwreal([165.77, 0.43], "result")
    )
)

IsoVecRes["N200"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([180.17, 0.43], "result"),
        "lc" => uwreal([180.43, 0.42], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([171.02, 0.44], "result"),
        "lc" => uwreal([172.41, 0.43], "result")
    )
)

IsoVecRes["D200"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([188.37, 0.38], "result"),
        "lc" => uwreal([188.69, 0.37], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([179.52, 0.39], "result"),
        "lc" => uwreal([180.91, 0.38], "result")
    )
)

IsoVecRes["N300"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([160.99, 0.59], "result"),
        "lc" => uwreal([161.08, 0.59], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([156.34, 0.59], "result"),
        "lc" => uwreal([156.89, 0.59], "result")
    )
)

IsoVecRes["A653"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([173.94, 0.36], "result"),
        "lc" => uwreal([176.25, 0.37], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([142.15, 0.35], "result"),
        "lc" => uwreal([151.27, 0.37], "result")
    )
)

IsoVecRes["H105"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([184.82, 0.50], "result"),
        "lc" => uwreal([186.01, 0.49], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([164.28, 0.53], "result"),
        "lc" => uwreal([169.09, 0.51], "result")
    )
)

IsoVecRes["H200"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([165.17, 0.91], "result"),
        "lc" => uwreal([165.44, 0.91], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([155.70, 0.89], "result"),
        "lc" => uwreal([157.21, 0.89], "result")
    )
)

IsoVecRes["D450"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([189.36, 0.26], "result"),
        "lc" => uwreal([190.03, 0.27], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([174.95, 0.26], "result"),
        "lc" => uwreal([177.79, 0.26], "result")
    )
)

IsoVecRes["D452"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([194.96, 0.33], "result"),
        "lc" => uwreal([195.61, 0.33], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([181.21, 0.34], "result"),
        "lc" => uwreal([183.97, 0.34], "result")
    )
)

IsoVecRes["E250"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([194.75, 0.26], "result"),
        "lc" => uwreal([194.96, 0.26], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([186.36, 0.27], "result"),
        "lc" => uwreal([187.61, 0.26], "result")
    )
)

IsoVecRes["J303"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([179.51, 0.54], "result"),
        "lc" => uwreal([179.57, 0.55], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([175.24, 0.55], "result"),
        "lc" => uwreal([175.67, 0.55], "result")
    )
)

IsoVecRes["E300"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([188.05, 0.49], "result"),
        "lc" => uwreal([188.13, 0.49], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([183.96, 0.49], "result"),
        "lc" => uwreal([184.38, 0.50], "result")
    )
)

IsoVecRes["J500"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([162.00, 0.72], "result"),
        "lc" => uwreal([162.04, 0.72], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([159.69, 0.72], "result"),
        "lc" => uwreal([159.97, 0.72], "result")
    )
)

IsoVecRes["J501"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([170.16, 0.98], "result"),
        "lc" => uwreal([170.15, 0.98], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([167.92, 0.98], "result"),
        "lc" => uwreal([168.13, 0.98], "result")
    )
)

IsoVecRes["N451"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([183.40, 0.28], "result"),
        "lc" => uwreal([184.05, 0.28], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([168.49, 0.27], "result"),
        "lc" => uwreal([171.40, 0.27], "result")
    )
)
    


# Iso-scalar contributions (Table IX)
IsoScaRes = Dict{String, Dict{String, Dict{String, uwreal}}}()

IsoScaRes["H101"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([57.36, 0.13], "result"),
        "lc" => uwreal([57.78, 0.13], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([50.05, 0.13], "result"),
        "lc" => uwreal([51.78, 0.13], "result")
    )
)

IsoScaRes["H102"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([55.30, 0.16], "result"),
        "lc" => uwreal([55.71, 0.16], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([47.94, 0.16], "result"),
        "lc" => uwreal([49.70, 0.16], "result")
    )
)

IsoScaRes["N101"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([53.55, 0.11], "result"),
        "lc" => uwreal([53.97, 0.11], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([46.18, 0.11], "result"),
        "lc" => uwreal([47.99, 0.11], "result")
    )
)

IsoScaRes["C101"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([52.67, 0.11], "result"),
        "lc" => uwreal([53.08, 0.11], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([45.39, 0.11], "result"),
        "lc" => uwreal([47.18, 0.11], "result")
    )
)

IsoScaRes["B450"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([56.04, 0.13], "result"),
        "lc" => uwreal([56.27, 0.13], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([50.84, 0.13], "result"),
        "lc" => uwreal([51.89, 0.13], "result")
    )
)

IsoScaRes["N202"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([56.05, 0.23], "result"),
        "lc" => uwreal([56.15, 0.23], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([52.79, 0.22], "result"),
        "lc" => uwreal([53.31, 0.23], "result")
    )
)

IsoScaRes["N203"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([53.41, 0.13], "result"),
        "lc" => uwreal([53.50, 0.13], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([50.12, 0.13], "result"),
        "lc" => uwreal([50.65, 0.13], "result")
    )
)

IsoScaRes["N200"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([51.61, 0.11], "result"),
        "lc" => uwreal([51.70, 0.10], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([48.35, 0.11], "result"),
        "lc" => uwreal([48.88, 0.10], "result")
    )
)

IsoScaRes["D200"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([50.36, 0.10], "result"),
        "lc" => uwreal([50.46, 0.10], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([47.11, 0.10], "result"),
        "lc" => uwreal([47.67, 0.09], "result")
    )
)

IsoScaRes["N300"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([53.66, 0.20], "result"),
        "lc" => uwreal([53.69, 0.20], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([52.11, 0.20], "result"),
        "lc" => uwreal([52.30, 0.20], "result")
    )
)

IsoScaRes["A653"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([57.98, 0.12], "result"),
        "lc" => uwreal([58.75, 0.12], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([47.38, 0.12], "result"),
        "lc" => uwreal([50.42, 0.12], "result")
    )
)

IsoScaRes["H105"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([53.16, 0.16], "result"),
        "lc" => uwreal([53.57, 0.15], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([45.83, 0.15], "result"),
        "lc" => uwreal([47.61, 0.15], "result")
    )
)

IsoScaRes["D450"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([51.47, 0.06], "result"),
        "lc" => uwreal([51.70, 0.07], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([46.23, 0.06], "result"),
        "lc" => uwreal([47.33, 0.06], "result")
    )
)

IsoScaRes["D452"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([50.90, 0.10], "result"),
        "lc" => uwreal([51.12, 0.10], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([45.74, 0.10], "result"),
        "lc" => uwreal([46.84, 0.10], "result")
    )
)

IsoScaRes["H200"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([55.06, 0.30], "result"),
        "lc" => uwreal([55.15, 0.30], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([51.90, 0.30], "result"),
        "lc" => uwreal([52.40, 0.30], "result")
    )
)

IsoScaRes["E250"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([49.65, 0.09], "result"),
        "lc" => uwreal([49.76, 0.09], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([46.45, 0.09], "result"),
        "lc" => uwreal([47.01, 0.09], "result")
    )
)

IsoScaRes["J303"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([49.80, 0.12], "result"),
        "lc" => uwreal([49.82, 0.12], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([48.25, 0.12], "result"),
        "lc" => uwreal([48.43, 0.12], "result")
    )
)

IsoScaRes["E300"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([48.77, 0.08], "result"),
        "lc" => uwreal([48.80, 0.08], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([47.24, 0.08], "result"),
        "lc" => uwreal([47.44, 0.08], "result")
    )
)

IsoScaRes["J500"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([54.00, 0.24], "result"),
        "lc" => uwreal([54.01, 0.24], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([53.23, 0.24], "result"),
        "lc" => uwreal([53.32, 0.24], "result")
    )
)

IsoScaRes["N451"] = Dict(
    "1old" => Dict(
        "ll" => uwreal([52.80, 0.06], "result"),
        "lc" => uwreal([53.01, 0.06], "result")
    ),
    "2" => Dict(
        "ll" => uwreal([47.51, 0.06], "result"),
        "lc" => uwreal([48.60, 0.06], "result")
    )
)


# Finite-Volume corrections (Table VI)
FVCRes = Dict(
    "A653" => Dict(
        "HP(t<t*)" => uwreal([0.80, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([1.19, 0.04], "HP(t>t*)")
    ),

    "H101" => Dict(
        "HP(t<t*)" => uwreal([0.73, 0.02], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.13, 0.00], "HP(t>t*)")
    ),
    
    "H102" => Dict(
        "HP(t<t*)" => uwreal([0.62, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.57, 0.02], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.19, 0.00], "Kaon loop")
    ),
    
    "H105" => Dict(
        "HP(t<t*)" => uwreal([0.54, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([2.10, 0.06], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.14, 0.00], "Kaon loop")
    ),
    
    "N101" => Dict(
        "HP(t<t*)" => uwreal([0.29, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.00, 0.00], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.01, 0.00], "Kaon loop")
    ),
    
    "C101" => Dict(
        "HP(t<t*)" => uwreal([0.73, 0.02], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.03, 0.00], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.01, 0.00], "Kaon loop")
    ),
    
    "B450" => Dict(
        "HP(t<t*)" => uwreal([0.63, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([1.21, 0.03], "HP(t>t*)")
    ),
    
    "S400" => Dict(
        "HP(t<t*)" => uwreal([0.50, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([1.87, 0.04], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.34, 0.00], "Kaon loop")
    ),
    
    "N451" => Dict(
        "HP(t<t*)" => uwreal([0.54, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.01, 0.00], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.02, 0.00], "Kaon loop")
    ),
    
    "D450" => Dict(
        "HP(t<t*)" => uwreal([0.32, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.00, 0.00], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.00, 0.00], "Kaon loop")
    ),
    
    "D452" => Dict(
        "HP(t<t*)" => uwreal([0.88, 0.02], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.07, 0.00], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.00, 0.00], "Kaon loop")
    ),
    
    "H200" => Dict(
        "HP(t<t*)" => uwreal([0.45, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([4.14, 0.09], "HP(t>t*)")
    ),
    
    "N202" => Dict(
        "HP(t<t*)" => uwreal([0.44, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.01, 0.00], "HP(t>t*)")
    ),
    
    "N203" => Dict(
        "HP(t<t*)" => uwreal([0.57, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.11, 0.00], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.09, 0.00], "Kaon loop")
    ),
    
    "N200" => Dict(
        "HP(t<t*)" => uwreal([0.66, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.81, 0.02], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.07, 0.00], "Kaon loop")
    ),
    
    "D200" => Dict(
        "HP(t<t*)" => uwreal([0.96, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.12, 0.00], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.01, 0.00], "Kaon loop")
    ),
    
    "E250" => Dict(
        "HP(t<t*)" => uwreal([0.53, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.00, 0.00], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.00, 0.00], "Kaon loop")
    ),
    
    "N300" => Dict(
        "HP(t<t*)" => uwreal([0.63, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([1.37, 0.03], "HP(t>t*)")
    ),
    
    "N302" => Dict(
        "HP(t<t*)" => uwreal([0.45, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([2.29, 0.05], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.33, 0.00], "Kaon loop")
    ),
    
    "J303" => Dict(
        "HP(t<t*)" => uwreal([0.81, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.93, 0.02], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.05, 0.00], "Kaon loop")
    ),
    
    "E300" => Dict(
        "HP(t<t*)" => uwreal([0.76, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.02, 0.00], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.00, 0.00], "Kaon loop")
    ),
    
    "J500" => Dict(
        "HP(t<t*)" => uwreal([0.74, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([0.92, 0.03], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.00, 0.00], "Kaon loop")
    ),
    
    "J501" => Dict(
        "HP(t<t*)" => uwreal([0.43, 0.01], "HP(t<t*)"),
        "HP(t>t*)" => uwreal([2.02, 0.05], "HP(t>t*)"),
        "Kaon loop" => uwreal([0.29, 0.00], "Kaon loop")
    )
)


##------ 

diag = "LO"
comp = "88"
wind = "ID"

STD_DERIV = true
FVCbool   = false
t0SHIFT   = true

IMPR_SET = ["1old","2"]

if comp == "33"
    paperRes = IsoVecRes
    mult = 1
elseif comp == "88"
    paperRes = IsoScaRes
    mult = 1/3
end

ensBool = [ensid in keys(paperRes) for ensid in getfield.(ensInfo,:id)]

if STD_DERIV
    path_bdio = joinpath(julia_script_directory, "..", "..", "ObsBDIOstd")
else
    path_bdio = joinpath(julia_script_directory, "..", "..", "ObsBDIO")
end

Res = Dict()
for ens in ensInfo[ensBool]
    println("- Reading data ensemble: $(ens.id)")

    println("   - Reading t0...")

    t0 = BDIOread_t0(path_bdio, ens)

    println("   - Reading HVP...    [applying systematics]")

    HVP = Dict()
    for impr_set in IMPR_SET
        hvp, info = BDIOread_HVPens(path_bdio,diag,wind,ens,impr_set,info=true)
    
        HVP[impr_set] = apply_syst_HVP(hvp,info["HVPsyst"],diag,wind,ens.id)
    end

    if FVCbool
        println("   - Reading FVC...    [applying systematics]")

        fvc = BDIOread_FVCens(path_bdio,diag,wind,ens)
        
        FVC = apply_syst_FVC(fvc,diag,wind,ens.id,IMPR_SET=IMPR_SET)

        println("   - aµ = HVP + FVC")

        Res_ = HVP_VolCorrect(HVP,FVC,diag,IMPR_SET=IMPR_SET)
    else
        Res_ = HVP
    end
    
    if t0SHIFT
        println("   - Performing t0 shift...")

        for impr_set in IMPR_SET
            for key in ["g$(comp)_ll","g$(comp)_lc"]
                uwerr(Res_[impr_set][key]); der = mchist(Res_[impr_set][key], "sqrtt0 [fm]")[1] / artificial_err
                Res_[impr_set][key] = Res_[impr_set][key] + value(sqrtt0_ph_CLS - sqrtt0_ph_Regensburg) * der
            end
        end
    end

    Res[ens.id] = Res_
end

println("\n\n\n")

# Initialize the LaTeX table string
latex_table = """
\\begin{table}[h]
\\centering
\\begin{tabular}{|c|c|c|c|c|}
\\hline
ens & Set 1 ll & Set 1 lc & Set 2 ll & Set 2 lc \\\\
\\hline
\\hline

"""

# Iterate over each ens.id
for ens in ensInfo[ensBool]
    row = ens.id * " & "
    for impr_set in IMPR_SET
        for discr in ["ll","lc"]
            res_value = mult*Res[ens.id][impr_set]["g$(comp)_$(discr)"]; uwerr(res_value)
            if FVCbool
                pap_value = paperRes[ens.id][impr_set][discr]; uwerr(pap_value)
            else
                if comp == "33"
                    if ens.kappa_l == ens.kappa_s
                        pap_value = paperRes[ens.id][impr_set][discr] - (FVCRes[ens.id]["HP(t<t*)"]+FVCRes[ens.id]["HP(t>t*)"]); uwerr(pap_value)
                    else
                        pap_value = paperRes[ens.id][impr_set][discr] - (FVCRes[ens.id]["HP(t<t*)"]+FVCRes[ens.id]["HP(t>t*)"]+FVCRes[ens.id]["Kaon loop"]); uwerr(pap_value)
                    end
                elseif comp == "88"
                    if ens.kappa_l == ens.kappa_s
                        pap_value = paperRes[ens.id][impr_set][discr] - (1/3) * (FVCRes[ens.id]["HP(t<t*)"]+FVCRes[ens.id]["HP(t>t*)"]); uwerr(pap_value)
                    else
                        pap_value = paperRes[ens.id][impr_set][discr] - 3/2 * (FVCRes[ens.id]["Kaon loop"]); uwerr(pap_value) # 3/2*
                    end
                end
            end
            difference = res_value - pap_value; uwerr(difference)
            sigma = abs(value(difference))/err(difference)
            cell = "\\begin{tabular}{@{}c@{}} $(round(value(pap_value),digits=2)) \$\\pm\$ $(round(err(pap_value),digits=2)) \\\\ $(round(value(res_value),digits=2)) \$\\pm\$ $(round(err(res_value),digits=2)) \\\\ $(round(sigma,digits=2)) \\end{tabular}"
            row *= cell * " & "
        end
    end
    # Remove the last " & " and add a line break
    row = row[1:end-3] * " \\\\\n \\hline \n"
    latex_table *= row
end

# Close the LaTeX table
if comp == "33"
    latex_table *= """
    \\hline
    \\end{tabular}
    \\caption{Comparison iso-vector}
    \\end{table}
    """
elseif comp == "88"
    latex_table *= """
    \\hline
    \\end{tabular}
    \\caption{Comparison iso-scalar}
    \\end{table}
    """
end

# Print or save the LaTeX table
println(latex_table)



## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##
## <<------------------------------------------------------------------------------------------------------------------------->> ##

##==========================> 1904.03120v1: 2019 (g-2)µ paper (Tables VIII and IX) <==========================##


using ALPHAio, ADerrors
import ADerrors: err

using BDIO


julia_script_directory = @__DIR__
path_bdio = joinpath(julia_script_directory, "..", "..", "ObsBDIO")

ensList = ["H101", "H102", "N101", "C101", "B450", "S400", "N202", "N203", "N200", "D200", "E250", "N300", "N302", "J303"]
ensInfo = EnsInfo.(ensList)

kcd_in = Dict(
    "A653" => Dict("kappaC" => 0.119743, "kappaC_err" => 0.000017, "kappaC_sim" => 0.119743, "kappaC_sim_plus" => 0.119793),
    "A654" => Dict("kappaC" => 0.120079, "kappaC_err" => 0.000025, "kappaC_sim" => 0.120177, "kappaC_sim_plus" => 0.120227),
    "H101" => Dict("kappaC" => 0.122897, "kappaC_err" => 0.000018, "kappaC_sim" => 0.122908, "kappaC_sim_plus" => 0.122938),
    "H102" => Dict("kappaC" => 0.123041, "kappaC_err" => 0.000026, "kappaC_sim" => 0.123050, "kappaC_sim_plus" => 0.123080),
    "U101" => Dict("kappaC" => 0.123244, "kappaC_err" => 0.000019, "kappaC_sim" => 0.123251, "kappaC_sim_plus" => 0.123281),
    "H105" => Dict("kappaC" => 0.123244, "kappaC_err" => 0.000019, "kappaC_sim" => 0.123251, "kappaC_sim_plus" => 0.123281),
    "N101" => Dict("kappaC" => 0.123244, "kappaC_err" => 0.000019, "kappaC_sim" => 0.123251, "kappaC_sim_plus" => 0.123281),
    "C101" => Dict("kappaC" => 0.123362, "kappaC_err" => 0.000015, "kappaC_sim" => 0.123367, "kappaC_sim_plus" => 0.123397),
    "B450" => Dict("kappaC" => 0.125093, "kappaC_err" => 0.000017, "kappaC_sim" => 0.125089, "kappaC_sim_plus" => 0.125129),
    "S400" => Dict("kappaC" => 0.125252, "kappaC_err" => 0.000020, "kappaC_sim" => 0.125267, "kappaC_sim_plus" => 0.125317),
    "N401" => Dict("kappaC" => 0.125439, "kappaC_err" => 0.000015, "kappaC_sim" => 0.125447, "kappaC_sim_plus" => 0.125477),
    "N451" => Dict("kappaC" => 0.125439, "kappaC_err" => 0.000015, "kappaC_sim" => 0.125447, "kappaC_sim_plus" => 0.125477),
    "D450" => Dict("kappaC" => 0.125585, "kappaC_err" => 0.000007, "kappaC_sim" => 0.125585, "kappaC_sim_plus" => 0.125635),
    "D452" => Dict("kappaC" => 0.125645, "kappaC_err" => 0.000005, "kappaC_sim" => 0.125640, "kappaC_sim_plus" => 0.125690),
    "H200" => Dict("kappaC" => 0.127579, "kappaC_err" => 0.000016, "kappaC_sim" => 0.127626, "kappaC_sim_plus" => 0.127666),
    "N202" => Dict("kappaC" => 0.127579, "kappaC_err" => 0.000016, "kappaC_sim" => 0.127626, "kappaC_sim_plus" => 0.127666),
    "N203" => Dict("kappaC" => 0.127714, "kappaC_err" => 0.000011, "kappaC_sim" => 0.127713, "kappaC_sim_plus" => 0.127733),
    "N200" => Dict("kappaC" => 0.127858, "kappaC_err" => 0.000007, "kappaC_sim" => 0.127859, "kappaC_sim_plus" => 0.127879),
    "D200" => Dict("kappaC" => 0.127986, "kappaC_err" => 0.000006, "kappaC_sim" => 0.127986, "kappaC_sim_plus" => 0.127956),
    "E250" => Dict("kappaC" => 0.128052, "kappaC_err" => 0.000005, "kappaC_sim" => 0.128054, "kappaC_sim_plus" => 0.128064),
    "N300" => Dict("kappaC" => 0.130099, "kappaC_err" => 0.000018, "kappaC_sim" => 0.130099, "kappaC_sim_plus" => 0.130149),
    "N302" => Dict("kappaC" => 0.130247, "kappaC_err" => 0.000009, "kappaC_sim" => 0.130243, "kappaC_sim_plus" => 0.130263),
    "J303" => Dict("kappaC" => 0.130362, "kappaC_err" => 0.000009, "kappaC_sim" => 0.130362, "kappaC_sim_plus" => 0.130382),
    "E300" => Dict("kappaC" => 0.130432, "kappaC_err" => 0.000010, "kappaC_sim" => 0.130421, "kappaC_sim_plus" => 0.130400),
    "J500" => Dict("kappaC" => 0.131663, "kappaC_err" => 0.000016, "kappaC_sim" => 0.131644, "kappaC_sim_plus" => 0.131600)
)

data_dict = Dict(
    "H101" => Dict(
        "interpolation" => 0.8615,
        "simulated" => [
            uwreal([0.8513, 0.0006], "sim H101"),
            uwreal([0.8614, 0.0006], "sim H101"),
            uwreal([0.8714, 0.0006], "sim H101"),
            uwreal([0.8813, 0.0006], "sim H101")
        ]
    ),
    "H102" => Dict(
        "interpolation" => 0.8615,
        "simulated" => [
            uwreal([0.8528, 0.0009], "sim H102"),
            uwreal([0.8629, 0.0009], "sim H102"),
            uwreal([0.8729, 0.0009], "sim H102"),
            uwreal([0.8828, 0.0009], "sim H102")
        ]
    ),
    "N101" => Dict(
        "interpolation" => 0.8615,
        "simulated" => [
            uwreal([0.8495, 0.0006], "sim N101"),
            uwreal([0.8563, 0.0006], "sim N101"),
            uwreal([0.8630, 0.0006], "sim N101"),
            uwreal([0.8697, 0.0006], "sim N101")
        ]
    ),
    "C101" => Dict(
        "interpolation" => 0.8615,
        "simulated" => [
            uwreal([0.8467, 0.0004], "sim C101"),
            uwreal([0.8534, 0.0004], "sim C101"),
            uwreal([0.8602, 0.0004], "sim C101"),
            uwreal([0.8669, 0.0004], "sim C101")
        ]
    ),
    "B450" => Dict(
        "interpolation" => 0.7615,
        "simulated" => [
            uwreal([0.7543, 0.0008], "sim B450"),
            uwreal([0.7614, 0.0008], "sim B450"),
            uwreal([0.7683, 0.0008], "sim B450"),
            uwreal([0.7752, 0.0008], "sim B450")
        ]
    ),
    "S400" => Dict(
        "interpolation" => 0.7615,
        "simulated" => [
            uwreal([0.7457, 0.0007], "sim S400"),
            uwreal([0.7528, 0.0007], "sim S400"),
            uwreal([0.7599, 0.0007], "sim S400"),
            uwreal([0.7669, 0.0007], "sim S400")
        ]
    ),
    "N202" => Dict(
        "interpolation" => 0.6410,
        "simulated" => [
            uwreal([0.6347, 0.0006], "sim N202"),
            uwreal([0.6384, 0.0006], "sim N202"),
            uwreal([0.6421, 0.0006], "sim N202"),
            uwreal([0.6458, 0.0006], "sim N202")
        ]
    ),
    "N203" => Dict(
        "interpolation" => 0.6410,
        "simulated" => [
            uwreal([0.6341, 0.0004], "sim N203"),
            uwreal([0.6416, 0.0004], "sim N203"),
            uwreal([0.6490, 0.0004], "sim N203"),
            uwreal([0.6563, 0.0004], "sim N203")
        ]
    ),
    "N200" => Dict(
        "interpolation" => 0.6410,
        "simulated" => [
            uwreal([0.6320, 0.0003], "sim N200"),
            uwreal([0.6395, 0.0003], "sim N200"),
            uwreal([0.6469, 0.0003], "sim N200"),
            uwreal([0.6543, 0.0003], "sim N200")
        ]
    ),
    "D200" => Dict(
        "interpolation" => 0.6410,
        "simulated" => [
            uwreal([0.6330, 0.0002], "sim D200"),
            uwreal([0.6367, 0.0002], "sim D200"),
            uwreal([0.6405, 0.0002], "sim D200"),
            uwreal([0.6442, 0.0002], "sim D200")
        ]
    ),
    "E250" => Dict(
        "interpolation" => 0.6410,
        "simulated" => [
            uwreal([0.6317, 0.0002], "sim E250"),
            uwreal([0.6392, 0.0002], "sim E250"),
            uwreal([0.6467, 0.0002], "sim E250"),
            uwreal([0.6541, 0.0002], "sim E250")
        ]
    ),
    "N300" => Dict(
        "interpolation" => 0.4969,
        "simulated" => [
            uwreal([0.5126, 0.0007], "sim N300"),
            uwreal([0.5008, 0.0007], "sim N300"),
            uwreal([0.4888, 0.0007], "sim N300"),
            uwreal([0.4766, 0.0007], "sim N300")
        ]
    ),
    "N302" => Dict(
        "interpolation" => 0.4969,
        "simulated" => [
            uwreal([0.5066, 0.0003], "sim N302"),
            uwreal([0.4947, 0.0003], "sim N302"),
            uwreal([0.4827, 0.0003], "sim N302"),
            uwreal([0.4704, 0.0003], "sim N302")
        ]
    ),
    "J303" => Dict(
        "interpolation" => 0.4969,
        "simulated" => [
            uwreal([0.5033, 0.0004], "sim J303"),
            uwreal([0.4954, 0.0004], "sim J303"),
            uwreal([0.4873, 0.0004], "sim J303"),
            uwreal([0.4792, 0.0003], "sim J303")
        ]
    )
)


fb = BDIO_open(joinpath(path_bdio,"kappaC_tar","mDs_prime"),"r")
amDs_SU3 = Dict{String,uwreal}()
i = 0
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    nobs = d["nobs"]
    dims = d["dimensions"]
    i=i+1
    if i==1
        mDs_ph = ALPHAdobs_read_next(fb)
    else
        ks = collect(d["keys"])
        amDs_SU3 = ALPHAdobs_read_next(fb, keys=ks)
    end
end


KappaC_tar = Dict{String,uwreal}(); myKappaC_tar = Dict{String,uwreal}()
mDs_sim = Dict{String,Vector{uwreal}}(); mymDs_sim = Dict{String,Vector{uwreal}}()
for ens in ensInfo
    fb = BDIO_open(joinpath(path_bdio,"kappaC_tar",ens.id,"$(ens.id)_kappaC"),"r")
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        nobs = d["nobs"]
        dims = d["dimensions"]
        myKappaC_tar[ens.id] = ALPHAdobs_read_next(fb)
    end
    uwerr(myKappaC_tar[ens.id])

    fb = BDIO_open(joinpath(path_bdio,"kappaC_tar",ens.id,"$(ens.id)_mDs"),"r")
    while ALPHAdobs_next_p(fb)
        d = ALPHAdobs_read_parameters(fb)
        nobs = d["nobs"]
        dims = d["dimensions"]
        sz = tuple(d["size"]...)
        mymDs_sim[ens.id] = ALPHAdobs_read_next(fb, size=sz)
    end
    uwerr.(mymDs_sim[ens.id])

    KappaC_tar[ens.id] = uwreal([kcd_in[ens.id]["kappaC"],kcd_in[ens.id]["kappaC_err"]],"kappaC_int"); uwerr(KappaC_tar[ens.id])
    mDs_sim[ens.id] = data_dict[ens.id]["simulated"]; uwerr.(mDs_sim[ens.id])
end


##


# Initialize the LaTeX table string
latex_table = """
\\begin{table}[h]
\\centering
\\begin{tabular}{|c|c|c c c c|c|}
\\hline
ens & a mDs (int) & a mDs (sim) &  &  &  & kappaC (sim) \\\\
\\hline
\\hline

"""

# Iterate over each ens.id
for ens in ensInfo
    row = ens.id * " & "

    cell = "\\begin{tabular}{@{}c@{}} $(data_dict[ens.id]["interpolation"]) \\\\ $(round(value(amDs_SU3["$(ens.beta)"]),digits=4)) \\\\ ... \\end{tabular}"
    row *= cell * " & "

    difference = mDs_sim[ens.id] .- mymDs_sim[ens.id]; uwerr.(difference)
    sigma = abs.(value.(difference))./err.(difference)
    # cell = "\\begin{tabular}{@{}c@{}} $(round(value(mDs_sim[ens.id][1]),digits=4)) +/- $(round(err(mDs_sim[ens.id][1]),digits=4)) & $(round(value(mDs_sim[ens.id][2]),digits=4)) +/- $(round(err(mDs_sim[ens.id][2]),digits=4)) & $(round(value(mDs_sim[ens.id][3]),digits=4)) +/- $(round(err(mDs_sim[ens.id][3]),digits=4)) & $(round(value(mDs_sim[ens.id][4]),digits=4)) +/- $(round(err(mDs_sim[ens.id][4]),digits=4)) \\\\ $(round(value(mymDs_sim[ens.id][1]),digits=4)) +/- $(round(err(mymDs_sim[ens.id][1]),digits=4)) & $(round(value(mymDs_sim[ens.id][2]),digits=4)) +/- $(round(err(mymDs_sim[ens.id][2]),digits=4)) & $(round(value(mymDs_sim[ens.id][3]),digits=4)) +/- $(round(err(mymDs_sim[ens.id][3]),digits=4)) & $(round(value(mymDs_sim[ens.id][4]),digits=4)) +/- $(round(err(mymDs_sim[ens.id][4]),digits=4)) \\\\ $(round(sigma[1],digits=2)) & $(round(sigma[2],digits=2)) & $(round(sigma[3],digits=2)) & $(round(sigma[4],digits=2)) \\end{tabular}"
    
    for i=collect(1:4)
        cell = "\\begin{tabular}{@{}c@{}} $(round(value(mDs_sim[ens.id][i]),digits=4)) +/- $(round(err(mDs_sim[ens.id][i]),digits=4)) \\\\ $(round(value(mymDs_sim[ens.id][i]),digits=4)) +/- $(round(err(mymDs_sim[ens.id][i]),digits=4)) \\\\ $(round(sigma[i],digits=2)) \\end{tabular}"
        row *= cell * " & "
    end

    difference = KappaC_tar[ens.id] - myKappaC_tar[ens.id]; uwerr(difference)
    sigma = abs(value(difference))/err(difference)
    cell = "\\begin{tabular}{@{}c@{}} $(round(value(KappaC_tar[ens.id]),digits=6)) +/- $(round(err(KappaC_tar[ens.id]),digits=6)) \\\\ $(round(value(myKappaC_tar[ens.id]),digits=6)) +/- $(round(err(myKappaC_tar[ens.id]),digits=6)) \\\\ $(round(sigma,digits=2)) \\end{tabular}"
    row *= cell

    row = row * " \\\\\n \\hline \n"
    latex_table *= row
end

# Close the LaTeX table
latex_table *= """
\\end{tabular}
\\caption{Comparison with 2019: results strange-heavy data for kappaC}
\\end{table}
"""


# Print or save the LaTeX table
println(latex_table)



## Data from window paper


amu = Dict()
for ensid in keys(IsoVecRes)
    amu[ensid] = Dict()
    for impr_set = ["1old","2"]
        amu[ensid][impr_set] = Dict()
        for discr in ["ll","lc"]
            amu[ensid][impr_set]["g33_$discr"] = IsoVecRes[ensid][impr_set][discr]
            if ensid != "J501"
                amu[ensid][impr_set]["g88_$discr"] = IsoScaRes[ensid][impr_set][discr]
            end
        end
    end
end

ensList = []
for ensid in keys(IsoVecRes)
    ensid != "H200" ? push!(ensList,ensid) : nothing
end
ensInfo = EnsInfo.(ensList)

