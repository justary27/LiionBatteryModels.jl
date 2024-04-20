include("../constants/constants.jl")

using QuadGK

"""
Electrloyte concentration in negative electrode [68]
- Parameters: 
    - x: position in -ve electrode.
    - n0: params[1]
    - n2: params[2]
    - n3: params[3]
"""
function c₂ₙ(x, params::Vector{Float64})
    n₀ = params[1]; n₂ = params[2]; n₃ = params[3]

    n₀ + n₂ * (lₙ^2 - x^2) + n₃ * (lₙ^3 - x^3)
end

"""
Electrloyte concentration in separator [70]
- Parameters: 
    - x: position in the separator.
    - s0: params[1]
    - s1: params[2]
    - s2: params[3]
"""
function c₂ₛ(x, params::Vector{Float64})
    s₀ = params[1]; s₁ = params[2]; s₂ = params[3]

    s₀ + s₁ * (x - lₙ) + s₂ * (x - lₙ)^2
end

"""
Electrloyte concentration in positive electrode [69]
- Parameters: 
    - x: position in +ve electrode.
    - p0: params[1]
    - p2: params[2]
    - p3: params[3]
"""
function c₂ₚ(x, params::Vector{Float64})
    p₀ = params[1]; p₂ = params[2]; p₃ = params[3]

    p₀ + p₂ * (lₚ^2 - (L-x)^2) + p₃ * (lₚ^3 - (L-x)^3)
end

# Overall electrolyte concentration @ x

"""
Overall electrolyte concentration function
    """
function c₂(x, params::Vector{Float64})
    if x <= lₙ
        c₂ₙ(x, params)
    elseif x > lₙ && x < lₙ + lₛ
        c₂ₛ(x, params)
    else
        c₂ₚ(x, params)
    end
end

# Averaged electrolyte concentration in different battery domains.

"""
Average electrolyte concentration in negative electrode
"""
function c₂̄ₙ(params::Vector{Float64})
    integral, errest = quadgk(x->c₂ₙ(x, params), 0, lₙ)
    
    integral/lₙ
end

"""
Average electrolyte concentration in separator
"""
function c₂̄ₛ(params::Vector{Float64})
    integral, errest = quadgk(x->c₂ₛ(x, params), lₙ, lₛ + lₙ)
    
    integral/lₛ 
end

"""
Average electrolyte concentration in positive electrode
"""
function c₂̄ₚ(params::Vector{Float64})
    integral, errest = quadgk(x->c₂ₚ(x, params), lₛ + lₙ, L)
    
    integral/lₚ
end

"""
Electrolyte concentration @ x=0 in negative electrode
"""
function c₂ₙ₀(params::Vector{Float64})
    c₂(0, params)
end

"""
Interfacial electrolyte concentration in negative electrode
"""
function c₂ᵢₙ(params::Vector{Float64})
    c₂(lₙ, params)
end

"""
Electrolyte concentration in the middle of battery
"""
function c₂ₘ(params::Vector{Float64})
    c₂(lₙ + lₛ/2, params)
end

"""
Interfacial electrolyte concentration in positive electrode
"""
function c₂ᵢₚ(params::Vector{Float64})
    c₂(lₙ + lₛ, params)
end

"""
Electrolyte concentration @ x = L in positive electrode.
"""
function c₂ₚₗ(params::Vector{Float64})
    c₂(L, params)
end
