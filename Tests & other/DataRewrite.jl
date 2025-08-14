
using DelimitedFiles

julia_script_directory = @__DIR__

##

ens = "N452" # N452  H200
# gamma = "V1V1"

GAMMA = ["V1V1","V2V2","V3V3","V1T10","V2T20","V3T30"]

for comp in ["light","strange"]
    println("Starting $comp")
    for gamma in GAMMA
        println("- gamma: $gamma")
        pathr = joinpath(julia_script_directory,"..","..","LMEData","HVP_data","$ens","$(comp)_","raw_data","data_HVP_$gamma")
        pathw = joinpath(julia_script_directory,"..","..","LMEData","HVP_data","$ens","$(comp)","raw_data","data_HVP_$gamma")

        f = readdlm(pathr, ' ', '\n', skipstart=3)

        open(pathw, "w") do file
            write(file, "nb replicas  : 1\nnb confs replica 0 : 1000\n\n")
            for i=collect(1:1000)
                write(file, "# rep 0 conf $(f[(i-1)*129+1,5])\n\n")
                for j=collect(0:127)
                    write(file, "$(f[(i-1)*129+1+j+1,1])\t$(f[(i-1)*129+1+j+1,2])\t$(f[(i-1)*129+1+j+1,3])\n")
                end
                write(file, "\n")
            end
        end
    end
end

##

gamma = "V1V1"
pathr = joinpath(julia_script_directory,"..","..","LMEData","HVP_data","N452","light_lowstat","raw_data","data_HVP_$gamma")

f = readdlm(pathr, '\t', '\n', skipstart=3)

parse.(Float64, split(f[2,2]))
f[1,:]

##

GAMMA = ["V1V1","V2V2","V3V3","V1V1c","V2V2c","V3V3c","V1T10","V2T20","V3T30","V1cT10","V2cT20","V3cT30"]

for comp in ["light","strange"]
    println("Starting $comp")
    for gamma in GAMMA
        println("- gamma: $gamma")
        pathr = joinpath(julia_script_directory,"..","..","LMEData","HVP_data","$ens","$(comp)_lowstat","raw_data","data_HVP_$gamma")
        pathw = joinpath(julia_script_directory,"..","..","LMEData","HVP_data","$ens","$(comp)","raw_data","data_HVP_$gamma")

        f = readdlm(pathr, '\t', '\n', skipstart=3)

        open(pathw, "w") do file
            write(file, "nb replicas  : 1\nnb confs replica 0 : 1000\n\n")
            for i=collect(1:1000)
                write(file, "$(f[(i-1)*129+1,1])\n\n")
                for j=collect(0:127)
                    f_ = parse.(Float64, split(f[(i-1)*129+1+j+1,2]))
                    write(file, "$(f[(i-1)*129+1+j+1,1])\t$(f_[1])\t$(f_[2])\n")
                end
                write(file, "\n")
            end
        end
    end
end
