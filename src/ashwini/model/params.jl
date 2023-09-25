include("cell_voltage.jl")

# The model parameter dₖ₁

function dₙ₁(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I, dₙ₁ₚ, dₙ₂ₚ, params)
    F*lₙ*(1/k₂ₙ - 1/σ₁ₙ - 2*θ*log(c₂ᵢₙ(params)./c₂₀(params))/lₙ)/(2*R*T) + log(jₙ₀(c₁, c₁ᵣ, I, lₙ, dₙ₁ₚ, dₙ₂ₚ, params)./jₙ₀(c₁, c₁ᵣ, I, 0, dₙ₁ₚ, dₙ₂ₚ, params))
end

function dₚ₁(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I, dₚ₁ₚ, dₚ₂ₚ, params)
    F*lₚ*(1/k₂ₚ - 1/σ₁ₚ + 2*θ*log(c₂ᵢₚ(params)./c₂ₗ(params))/lₚ)/(2*R*T) + log(jₚ₀(c₁, c₁ᵣ, I, lₙ + lₛ, dₚ₁ₚ, dₚ₂ₚ, params)./jₚ₀(c₁, c₁ᵣ, I, L, dₚ₁ₚ, dₚ₂ₚ, params))
end

# The model parameter dₖ₂

function dₙ₂(dₙ₁)
    dₙ₁/(log((exp(dₙ₁) - 1)/dₙ₁))
end

function dₚ₂(dₚ₁)
    dₚ₁/(log((exp(dₚ₁) - 1)/dₚ₁))
end
