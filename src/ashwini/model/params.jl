include("cell_voltage.jl")

# The model parameter dₖ₁

function dₙ₁(I, cₛₙ₀::Float64, cₛₙₗ::Float64, c₂ᵢₙ::Float64, c₂ₙ₀::Float64, c₂ₙₐ::Float64)
    F*lₙ*(I/(2*k₂ₙ(c₂ₙₐ)) - I/(2*σ₁ₙ) - 2*θ*log(c₂ᵢₙ/c₂ₙ₀)/lₙ - (Uₙ(cₛₙₗ) - Uₙ(cₛₙ₀))/lₙ)/(2*R*T) + log(jₙ₀(cₛₙₗ, c₂ᵢₙ)./jₙ₀(cₛₙ₀ , c₂ₙ₀))
end

function dₙ₁(I, cₛₙ₀::Float64, cₛₙₗ::Float64, paramsₙ)
    F*lₙ*(I/(2*k₂ₙ(paramsₙ)) - I/(2*σ₁ₙ) - 2*θ*log(c₂ᵢₙ(paramsₙ)./c₂ₙ₀(paramsₙ))/lₙ - (Uₙ(cₛₙₗ) - Uₙ(cₛₙ₀))/lₙ)/(2*R*T) + log(jₙ₀(lₙ, cₛₙₗ, paramsₙ)./jₙ₀(0, cₛₙ₀ , paramsₙ))
end

function dₚ₁(I, cₛₚ₀::Float64, cₛₚₗ::Float64, c₂ᵢₚ::Float64, c₂ₚₗ::Float64, c₂ₚₐ::Float64)
    F*lₚ*(I/(2*k₂ₚ(c₂ₚₐ)) - I/(2*σ₁ₚ) + 2*θ*log(c₂ᵢₚ/c₂ₚₗ)/lₚ - (Uₚ(cₛₚₗ) - Uₚ(cₛₚ₀))/lₚ)/(2*R*T) + log(jₚ₀(cₛₚ₀, c₂ᵢₚ)/jₚ₀(cₛₚₗ, c₂ₚₗ))
end

function dₚ₁(I, cₛₚ₀::Float64, cₛₚₗ::Float64, paramsₚ)
    F*lₚ*(I/(2*k₂ₚ(paramsₚ)) - I/(2*σ₁ₚ) + 2*θ*log(c₂ᵢₚ(paramsₚ)./c₂ₚₗ(paramsₚ))/lₚ - (Uₚ(cₛₚₗ) - Uₚ(cₛₚ₀))/lₚ)/(2*R*T) + log(jₚ₀(lₙ + lₛ, cₛₚ₀, paramsₚ)./jₚ₀(L, cₛₚₗ, paramsₚ))
end

# The model parameter dₖ₂

function dₙ₂(dₙ₁)
    dₙ₁/(log((exp(dₙ₁) - 1)/dₙ₁))
end

function dₚ₂(dₚ₁)
    dₚ₁/(log((exp(dₚ₁) - 1)/dₚ₁))
end
