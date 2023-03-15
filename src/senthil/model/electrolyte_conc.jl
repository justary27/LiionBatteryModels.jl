include("interfacial_conc.jl")

using QuadGK

# Electrloyte concentrations @ x = 0, L and (lₙ + lₛ/2)
# Refer equations 60, 67 & 76

function c₂₀ₜ(c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    c₂ᵢₖ.c₂ᵢₙ + q₂ᵢₖ.q₂ᵢₙ*lₙ/(2*D₂ₙ)
end

function c₂ₗₜ(c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    c₂ᵢₖ.c₂ᵢₚ - q₂ᵢₖ.q₂ᵢₚ*lₚ/(2*D₂ₚ)
end

function c₂ₘₜ(c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    c₂ᵢₖ.c₂ᵢₙ - 3*q₂ᵢₖ.q₂ᵢₙ*lₛ/(8*D₂ₛ) - q₂ᵢₖ.q₂ᵢₚ*lₛ/(8*D₂ₛ)
end

# Electrloyte concentrations in different domains of battery
# Refer equations 55, 63 & 70

"""Electrloyte concentration in negative electrode"""
function c₂ₙ(x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    c₂ᵢₖ.c₂ᵢₙ + (lₙ^2 - x^2) * q₂ᵢₖ.q₂ᵢₙ/(2*lₙ*D₂ₙ)
end

"""Electrloyte concentration in separator"""
function c₂ₛ(x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    c₂ᵢₖ.c₂ᵢₙ - ((x-lₙ).*q₂ᵢₖ.q₂ᵢₙ)/D₂ₛ + ((x-lₙ)^2 * (q₂ᵢₖ.q₂ᵢₙ - q₂ᵢₖ.q₂ᵢₚ))/(2*lₛ.*D₂ₛ)
end

"""Electrloyte concentration in positive electrode"""
function c₂ₚ(x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    c₂ᵢₖ.c₂ᵢₚ - (lₚ^2 - (L-x)^2) * q₂ᵢₖ.q₂ᵢₚ/(2*lₚ*D₂ₚ)
end

# Overall electrolyte concentration @ x

"""Overall electrolyte concentration function"""
function c₂(x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    if x <= lₙ
        c₂ₙ(x, c₂ᵢₖ, q₂ᵢₖ)
    elseif x > lₙ && x <= lₙ + lₛ
        c₂ₛ(x, c₂ᵢₖ, q₂ᵢₖ)
    else
        c₂ₚ(x, c₂ᵢₖ, q₂ᵢₖ)
    end
end

# Averaged electrolyte concentration in different battery domains.
# TODO: Correct these functions

"""Average electrolyte concentration in negative electrode"""
function c₂̄ₙ(c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    integral, errest = quadgk(x->c₂(x, c₂ᵢₖ, q₂ᵢₖ), 0, lₙ)
    
    integral/lₙ
end

"""Average electrolyte concentration in separator"""
function c₂̄ₛ(c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    integral, errest = quadgk(x->c₂(x, c₂ᵢₖ, q₂ᵢₖ), lₙ, lₛ + lₙ)
    
    integral/lₛ 
end

"""Average electrolyte concentration in positive electrode"""
function c₂̄ₚ(c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    integral, errest = quadgk(x->c₂(x, c₂ᵢₖ, q₂ᵢₖ), lₛ + lₙ, L)
    
    integral/lₚ
end

