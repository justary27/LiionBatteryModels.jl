include("initial_params.jl")

"""
Interfacial potential at negative electrode
"""
function ϕ₂ᵢₙ(I, params)
    2*θ*log(c₂ᵢₙ(params)./c₂ₘ(params)) + I*lₛ/(2*k₂ₛ(params))
end

"""
Interfacial potential at positive electrode
"""
function ϕ₂ᵢₚ(I, params)
    2*θ*log(c₂ᵢₚ(params)./c₂ₘ(params)) - I*lₛ/(2*k₂ₛ(params))
end

# Electrloyte potentials in different domains of battery

# TODO: Fix below 3 fns

"""
Electrolyte potential in negative electrode
"""
function ϕ₂ₙ(x, I, params, args...)
    dₙ₁ = args[1]

    ϕ₂ᵢₙ(I, params) + 2*θ*log(c₂(x, params)./c₂ᵢₙ(params)) + I*(lₙ-x)/k₂ₙ(params) - I*exp(dₙ₁)*[(lₙ-x) - lₙ*(1-exp(-dₙ₁+dₙ₁*x/lₙ))/dₙ₁]/((k₂ₙ(params))*(exp(dₙ₁) - 1))
end

"""
Electrolyte potential in separator
"""
function ϕ₂ₛ(x, I, params, args...)
    2*θ*log(c₂(x, params)./c₂ₘ(params)) - I*(x - (lₙ + lₛ/2))/k₂ₛ(params)
end

"""
Electrolyte potential in positive electrode
"""
function ϕ₂ₚ(x, I, params, args...)
    dₚ₁ = args[1]
    
    ϕ₂ᵢₚ(I, params) + 2*θ*log(c₂(x, params)./c₂ᵢₚ(params)) + I*(lₙ+lₛ-x)/(k₂ₚ(params)) + I*exp(dₚ₁)*[(x-lₙ-lₛ) + lₚ*(-1+exp(-dₚ₁*(x-lₙ-lₛ)/lₚ))/dₚ₁]/((k₂ₚ(params))*(exp(dₚ₁) - 1))
end

# Overall electrolyte potential @x

"""
Overall electrolyte potential function
"""
function ϕ₂(x, I, params, args...)
    if x <= lₙ
        ϕ₂ₙ(x, I, params, args)
    elseif x > lₙ && x <= lₙ + lₛ
        ϕ₂ₛ(x, I, params, args)
    else
        ϕ₂ₚ(x, I, params, args)
    end
end
