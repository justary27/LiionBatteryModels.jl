include("solid_potential.jl")

# Cell voltage
# Refer equation 158.

"""Cell voltage across the battery"""
function V(I, c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    ϕ₁(c₁, c₁ᵣ, I, L, c₂ᵢₖ, q₂ᵢₖ) - ϕ₁(c₁, c₁ᵣ, I, 0, c₂ᵢₖ, q₂ᵢₖ)
end
