include("conductivity.jl")

# Interfacial potentials
# Refer equations 105 & 106

"""Interfacial potential at negative electrode"""
function ϕ₂ᵢₙ(c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux, I)
    2*θ*log.(c₂ᵢₖ.c₂ᵢₙ./c₂ₘₜ(c₂ᵢₖ, q₂ᵢₖ)) .+ I*lₛ./(2*k₂ₛ(c₂ᵢₖ, q₂ᵢₖ))
end

"""Interfacial potential at positive electrode"""
function ϕ₂ᵢₚ(c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux, I)
    2*θ*log.(c₂ᵢₖ.c₂ᵢₚ./c₂ₘₜ(c₂ᵢₖ, q₂ᵢₖ)) .- I*lₛ./(2*k₂ₛ(c₂ᵢₖ, q₂ᵢₖ))
end

# Electrloyte potentials in different domains of battery
# Refer equations 104, 113 & 116

"""Electrloyte potential in negative electrode"""
function ϕ₂ₙ(x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux, I)
    ϕ₂ᵢₙ(c₂ᵢₖ, q₂ᵢₖ, I) .+ 2*θ*log.(c₂ₙ(x, c₂ᵢₖ, q₂ᵢₖ)/c₂ᵢₖ.c₂ᵢₙ) .+ (I*(lₙ-x))./k₂ₙ(c₂ᵢₖ, q₂ᵢₖ) .- (I*(lₙ-x)^2)./(2*k₂ₙ(c₂ᵢₖ, q₂ᵢₖ)*lₙ)
end

"""Electrloyte potential in separator"""
function ϕ₂ₛ(x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux, I)
    2*θ*log.(c₂ₛ(x, c₂ᵢₖ, q₂ᵢₖ)/c₂ₘₜ(c₂ᵢₖ, q₂ᵢₖ)) .- I*(x-(lₙ+lₛ/2))./k₂ₛ(c₂ᵢₖ, q₂ᵢₖ)
end

"""Electrloyte potential in positive electrode"""
function ϕ₂ₚ(x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux, I)
    ϕ₂ᵢₚ(c₂ᵢₖ, q₂ᵢₖ, I) .+ 2*θ*log.(c₂ₚ(x, c₂ᵢₖ, q₂ᵢₖ)/c₂ᵢₖ.c₂ᵢₚ) .- (I*(x-lₙ-lₛ))./k₂ₚ(c₂ᵢₖ, q₂ᵢₖ) .+ (I*(x-lₙ-lₛ)^2)./(2*k₂ₚ(c₂ᵢₖ, q₂ᵢₖ)*lₚ)
end

# Overall electrolyte potential @x

"""Overall electrolyte potential function"""
function ϕ₂(x::Float64,c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux, I)
    if x <= lₙ
        ϕ₂ₙ(x, c₂ᵢₖ, q₂ᵢₖ, I)
    elseif x > lₙ && x <= lₙ + lₛ
        ϕ₂ₛ(x, c₂ᵢₖ, q₂ᵢₖ, I)
    else
        ϕ₂ₚ(x, c₂ᵢₖ, q₂ᵢₖ, I)
    end
end

