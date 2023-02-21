include("inerfacial_flux.jl")

# Refer equations 77 & 78.

"""
## Interfacial electrolyte concentration

### Fields
- c₂ᵢₙ: Interfacial electrolyte concentration in negative electrode.
- c₂ᵢₚ: Interfacial electrolyte concentration in positive electrode.
"""
struct InterfacialConc
    c₂ᵢₙ::Vector{Float64}
    c₂ᵢₚ::Vector{Float64}
end

function InterfacialConc(q₂ᵢₖ::InterfacialFlux)
    c₂ᵢₚ = c₂₀ + αᵢₙ * q₂ᵢₖ.q₂ᵢₙ + αᵢₚ * q₂ᵢₖ.q₂ᵢₚ
    c₂ᵢₙ = c₂ᵢₚ + lₛ * (q₂ᵢₖ.q₂ᵢₙ + q₂ᵢₖ.q₂ᵢₚ) / (2 * D₂ₛ)

    InterfacialConc(c₂ᵢₙ, c₂ᵢₚ)
end

