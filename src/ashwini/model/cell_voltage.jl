include("solid_potential.jl")

"""Cell voltage across the battery"""
function V(I, c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, paramsₙ, paramsₚ)
    ϕ₁(L, I, c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, paramsₚ) - ϕ₁(0, I, c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, paramsₙ)
end
