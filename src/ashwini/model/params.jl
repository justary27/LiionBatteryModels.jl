include("cell_voltage.jl")

# The model parameter dₖ₁

# TODO: Fix these functions

function dₙ₁(I, cₛₙ₀::Float64, cₛₙₗ::Float64, params)
    F*lₙ*(I/k₂ₙ(params) - I/σ₁ₙ - 2*θ*log(c₂ᵢₙ(params)./c₂ₙ₀(params))/lₙ - (Uₙ(cₛₙₗ) - Uₙ(cₛₙ₀))/lₙ)/(2*R*T) + log(jₙ₀(lₙ, cₛₙₗ, params)./jₙ₀(0, cₛₙ₀ , params))
end

function dₚ₁(I, cₛₚ₀::Float64, cₛₚₗ::Float64, params)
    F*lₚ(I/(2*k₂ₚ(params)) - I/(2*σ₁ₚ) + 2*θ*log(c₂ᵢₚ(params)./c₂ₚₗ(params))/lₚ - (Uₚ(cₛₚₗ) - Uₚ(cₛₚ₀))/lₚ)/(2*R*T) + log(jₚ₀(lₙ + lₛ, cₛₚ₀, params)./jₚ₀(L, cₛₚₗ, params))
end

# The model parameter dₖ₂

function dₙ₂(dₙ₁)
    dₙ₁/(log((exp(dₙ₁) - 1)/dₙ₁))
end

function dₚ₂(dₚ₁)
    dₚ₁/(log((exp(dₚ₁) - 1)/dₚ₁))
end
