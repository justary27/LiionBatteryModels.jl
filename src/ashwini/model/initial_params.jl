include("conductivity.jl")

# Initial value of dₖ₁

function dₙ₁₀(I)
    F*lₙ*I*(1/k₂ₙ(c₂₀) - 1/σ₁ₙ)/(4*R*T)
end

function dₚ₁₀(I)
    F*lₚ*I*(1/k₂ₚ(c₂₀) - 1/σ₁ₚ)/(4*R*T)
end
