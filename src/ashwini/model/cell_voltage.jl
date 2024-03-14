include("solid_potential.jl")

"""Cell voltage across the battery"""
function V(I, c₁::cₛ, c₁ᵣ::cₛᵣ, paramsₙ, paramsₚ, dₙ₁, dₚ₁)
    ϕ₁(L, I, c₁::cₛ, c₁ᵣ::cₛᵣ, paramsₚ, dₙ₁) - ϕ₁(0, I, c₁::cₛ, c₁ᵣ::cₛᵣ, paramsₙ, dₚ₁)
end
