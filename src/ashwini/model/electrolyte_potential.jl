include("initial_params.jl")

"""
Interfacial potential at negative electrode
"""
function ϕ₂ᵢₙ(I, paramsₙ, paramsₛ)
    2*θ*log(c₂ᵢₙ(paramsₙ)./c₂ₘ(paramsₛ)) + I*lₛ/(2*k₂ₛ(lₙ, paramsₙ))
end

function ϕ₂ᵢₙ(I, c₂ᵢₙ::Float64, c₂ₘ::Float64)
    2*θ*log(c₂ᵢₙ./c₂ₘ) + I*lₛ/(2*k₂ₛ(c₂ᵢₙ))
end

"""
Interfacial potential at positive electrode
"""
function ϕ₂ᵢₚ(I, paramsₚ, paramsₛ)
    2*θ*log(c₂ᵢₚ(paramsₚ)./c₂ₘ(paramsₛ)) - I*lₛ/(2*k₂ₛ(lₙ + lₛ, paramsₚ))
end

function ϕ₂ᵢₚ(I, c₂ᵢₚ::Float64, c₂ₘ::Float64)
    2*θ*log(c₂ᵢₚ./c₂ₘ) - I*lₛ/(2*k₂ₛ(c₂ᵢₚ))
end

# Electrloyte potentials in different domains of battery

# TODO: Fix below 3 fns

"""
Electrolyte potential in negative electrode
"""
function ϕ₂ₙ(x, I, paramsₙ::Vector{Float64}, paramsₛ::Vector{Float64}, args...)
    dₙ₁ = args[1]

    ϕ₂ᵢₙ(I, paramsₙ, paramsₛ) + 2*θ*log(c₂(x, paramsₙ)./c₂ᵢₙ(paramsₙ)) + I*(lₙ-x)/k₂ₙ(x, paramsₙ) - I*exp(dₙ₁)*((lₙ-x) - lₙ*(1-exp(-dₙ₁+dₙ₁*x/lₙ))/dₙ₁)/(k₂ₙ(x, paramsₙ)*(exp(dₙ₁) - 1))
end

function ϕ₂ₙ(x, I::Float64, c₂::Float64, c₂ᵢₙ::Float64, ϕ₂ᵢₙ::Float64, dₙ₁::Float64)
    ϕ₂ᵢₙ + 2*θ*log(c₂/c₂ᵢₙ) + I*(lₙ-x)/k₂ₙ(c₂) - I*exp(dₙ₁)*((lₙ-x) - lₙ*(1-exp(-dₙ₁+dₙ₁*x/lₙ))/dₙ₁)/(k₂ₙ(c₂)*(exp(dₙ₁) - 1))
end

"""
Electrolyte potential in separator
"""
function ϕ₂ₛ(x, I, paramsₛ, args...)
    2*θ*log(c₂(x, paramsₛ)./c₂ₘ(paramsₛ)) - I*(x - (lₙ + lₛ/2))/k₂ₛ(x, paramsₛ)
end

function ϕ₂ₛ(x, I, c₂::Float64, c₂ₘ::Float64)
    2*θ*log(c₂/c₂ₘ) - I*(x - (lₙ + lₛ/2))/k₂ₛ(c₂)
end

"""
Electrolyte potential in positive electrode
"""
function ϕ₂ₚ(x, I, paramsₚ, paramsₛ, args...)
    dₚ₁ = args[1]

    # println(c₂ᵢₚ(params))
    # println(c₂(x, params))

    ϕ₂ᵢₚ(I, paramsₚ, paramsₛ) + 2*θ*log(c₂(x, paramsₚ)./c₂ᵢₚ(paramsₚ)) + I*(lₙ+lₛ-x)/(k₂ₚ(x, paramsₚ)) + I*exp(dₚ₁)*((x-lₙ-lₛ) + lₚ*(-1+exp(-dₚ₁*(x-lₙ-lₛ)/lₚ))/dₚ₁)/((k₂ₚ(x, paramsₚ))*(exp(dₚ₁) - 1))
end

function ϕ₂ₚ(x, I, c₂::Float64, c₂ᵢₚ::Float64, ϕ₂ᵢₚ::Float64, dₚ₁::Float64)
    ϕ₂ᵢₚ + 2*θ*log(c₂/c₂ᵢₚ) + I*(lₙ+lₛ-x)/(k₂ₚ(c₂)) + I*exp(dₚ₁)*((x-lₙ-lₛ) + lₚ*(-1+exp(-dₚ₁*(x-lₙ-lₛ)/lₚ))/dₚ₁)/((k₂ₚ(c₂))*(exp(dₚ₁) - 1))
end

# Overall electrolyte potential @x

"""
Overall electrolyte potential function
"""
function ϕ₂(x, I, paramsₖ, paramsₛ, args...)
    if x <= lₙ
        ϕ₂ₙ(x, I, paramsₖ, paramsₛ, args...)
    elseif x > lₙ && x <= lₙ + lₛ
        ϕ₂ₛ(x, I, paramsₖ, args...)
    else
        ϕ₂ₚ(x, I, paramsₖ, paramsₛ, args...)
    end
end
