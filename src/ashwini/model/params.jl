include("cell_voltage.jl")

# The model parameter dₖ₁

# TODO: Fix these functions

function dₙ₁(c₁::cₛ, c₁ᵣ::cₛᵣ, I, dₙ₁ₚ, dₙ₂ₚ, params)
    F*lₙ*(1/k₂ₙ - 1/σ₁ₙ - 2*θ*log(c₂ᵢₙ(params)./c₂₀(params))/lₙ)/(2*R*T) + log(jₙ₀(c₁, c₁ᵣ, I, lₙ, params)./jₙ₀(c₁, c₁ᵣ, I, 0, params))
end

function dₚ₁(c₁::cₛ, c₁ᵣ::cₛᵣ, I, dₚ₁ₚ, dₚ₂ₚ, params)
    F*lₚ*[I/(2*k₂ₚ) - I/(2*σ₁ₚ) + 2*θ*log(c₂ᵢₚ(params)./c₂ₗ(params))/lₚ]/(2*R*T) + log(jₚ₀(c₁, c₁ᵣ, I, lₙ + lₛ, params)./jₚ₀(c₁, c₁ᵣ, I, L, params))
end

# The model parameter dₖ₂

function dₙ₂(dₙ₁)
    dₙ₁/(log((exp(dₙ₁) - 1)/dₙ₁))
end

function dₚ₂(dₚ₁)
    dₚ₁/(log((exp(dₚ₁) - 1)/dₚ₁))
end
