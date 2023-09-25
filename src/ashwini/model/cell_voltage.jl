include("solid_potential.jl")

"""Cell voltage across the battery"""
function V(I, c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, dₙ₁, dₙ₂, dₚ₁, dₚ₂, params)
    ϕ₁(L, I, c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, dₙ₁, dₙ₂, params) - ϕ₁(0, I, c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, dₚ₁, dₚ₂, params)
end
