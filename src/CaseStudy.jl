"""
LiionBatteryModels core, allows access to the following submodules:
- SenthilModel
- AshwiniModel
"""
module LiionBatteryModels
    
    include("senthil/senthil.jl")
    include("ashwini/ashwini.jl")

    export SenthilModel, AshwiniModel
end

export LiionBatteryModels
