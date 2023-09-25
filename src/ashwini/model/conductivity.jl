include("electrolyte_conc.jl")

# Effective electronic conductivity

"""Effective electronic conductivity in negative electode"""
σ₁ₙ = k₁ₙ * (ϵ₁ₙ ^ brugₙ)

"""Effective electronic conductivity in positive electode"""
σ₁ₚ = k₁ₚ * (ϵ₁ₚ ^ brugₚ)

# Ionic conductivity

"""Electrloyte conductivity"""
function k₂(c₂)
    1e-4*c₂*(-10.5 + 0.074*T - 6.69e-5*T^2 + 6.68e-4*c₂ - 1.78e-5*c₂*T + 2.8e-8*c₂*T^2 + 4.94e-7*c₂.^2 - 8.86e-10*T*c₂.^2).^2
end

"""Electrloyte conductivity in negative electrode"""
function k₂ₙ(params)
    k₂(c₂̄ₙ(params)) * (ϵ₂ₙ ^ brugₙ)
end

function k₂ₙ(c₂::Float64)
    k₂(c₂) * (ϵ₂ₙ ^ brugₙ)
end

"""Electrloyte conductivity in separator"""
function k₂ₛ(params)
    k₂(c₂̄ₛ(params)) * (ϵ₂ₛ ^ brugₛ)
end

function k₂ₛ(c₂::Float64)
    k₂(c₂) * (ϵ₂ₛ ^ brugₛ)
end

"""Electrloyte conductivity in positive electrode"""
function k₂ₚ(params)
    k₂(c₂̄ₚ(params)) * (ϵ₂ₚ ^ brugₚ)
end

function k₂ₚ(c₂::Float64)
    k₂(c₂) * (ϵ₂ₚ ^ brugₚ)
end
