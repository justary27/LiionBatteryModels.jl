include("../constants/constants.jl")

using ..Constants

"""
## Interfacial flux

### Fields
- q₂ᵢₙ: Interfacial flux in negative electrode.
- q₂ᵢₚ: Interfacial flux in positive electrode.
"""
struct InterfacialFlux
    q₂ᵢₙ::Vector{Float64}
    q₂ᵢₚ::Vector{Float64}
end

