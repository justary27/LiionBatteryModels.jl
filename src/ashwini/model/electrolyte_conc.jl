include("../constants/constants.jl")

using QuadGK

"""Electrloyte concentration in negative electrode"""
function c₂ₙ(x, params)
    n₀ = params[1]; n₂ = params[2]; n₃ = params[3]

    n₀ + n₂ * (lₙ^2 - x^2) + n₃ * (lₙ^3 - x^3)
end

"""Electrloyte concentration in separator"""
function c₂ₛ(x, params)
    s₀ = params[1]; s₁ = params[2]; s₂ = params[3]

    s₀ + s₁ * (x - lₙ) + s₂ * (x - lₙ)^2
end

"""Electrloyte concentration in positive electrode"""
function c₂ₚ(x, params)
    p₀ = params[1]; p₂ = params[2]; p₃ = params[3]

    p₀ + p₂ * (lₚ^2 - (L-x)^2) + p₃ * (lₚ^3 - (L-x)^3)
end

# Overall electrolyte concentration @ x

"""Overall electrolyte concentration function"""
function c₂(x, params)
    if x <= lₙ
        c₂ₙ(x, params)
    elseif x > lₙ && x <= lₙ + lₛ
        c₂ₛ(x, params)
    else
        c₂ₚ(x, params)
    end
end

# Averaged electrolyte concentration in different battery domains.

"""Average electrolyte concentration in negative electrode"""
function c₂̄ₙ(params)
    integral, errest = quadgk(x->c₂(x, params), 0, lₙ)
    
    integral/lₙ
end

"""Average electrolyte concentration in separator"""
function c₂̄ₛ(params)
    integral, errest = quadgk(x->c₂(x, params), lₙ, lₛ + lₙ)
    
    integral/lₛ 
end

"""Average electrolyte concentration in positive electrode"""
function c₂̄ₚ(params)
    integral, errest = quadgk(x->c₂(x, params), lₛ + lₙ, L)
    
    integral/lₚ
end

"""Electrolyte concentration @ x=0 in negative electrode"""
function c₂ₙ₀(params)
    c₂(0, params)
end

"""Interfacial electrolyte concentration in negative electrode"""
function c₂ᵢₙ(params)
    c₂(lₙ, params)
end

"""Electrolyte concentration in the middle of battery"""
function c₂ₘ(params)
    c₂(lₙ + lₛ/2, params)
end

"""Interfacial electrolyte concentration in positive electrode"""
function c₂ᵢₚ(params)
    c₂(lₙ + lₛ, params)
end

"""Electrolyte concentration @ x = L in positive electrode."""
function c₂ₗ(params)
    c₂(L, params)
end
