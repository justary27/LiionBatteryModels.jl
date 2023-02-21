include("solid_conc.jl")

# Over-potential independent rate pre-factor.
# Refer equation 156.

"""Over-potential in negative electrode"""
function jₙ₀(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I, x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    kₙ*(c₁ₙₘₐₓ - cₛₙ(c₁, c₁ᵣ, I)).^(0.5).*cₛₙ(c₁, c₁ᵣ, I).^(0.5).*c₂(x, c₂ᵢₖ, q₂ᵢₖ).^(0.5)
end

"""Over-potential in positive electrode"""
function jₚ₀(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I, x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    kₚ*(c₁ₙₘₐₓ - cₛₚ(c₁, c₁ᵣ, I)).^(0.5).*cₛₚ(c₁, c₁ᵣ, I).^(0.5).*c₂(x, c₂ᵢₖ, q₂ᵢₖ).^(0.5)
end

# State of Charge

"""State of Charge in negative electrode"""
function SOCₙ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    cₛₙ(c₁, c₁ᵣ, I)/c₁ₙₘₐₓ
end

"""State of Charge in positive electrode"""
function SOCₚ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    cₛₚ(c₁, c₁ᵣ, I)/c₁ₚₘₐₓ
end

function DoD(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    (SOCₚ(c₁, c₁ᵣ, I) - SOCₚₘᵢₙ)/(SOCₚₘₐₓ - SOCₚₘᵢₙ)
end

# Electrode open circuit potential

"""Open circuit potential in negative electrode"""
function Uₙ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    SOCₙ = SOCₙ(c₁, c₁ᵣ, I)
    0.13966+0.68920.*exp(-49.20361.*SOCₙ)+0.41903.*exp(-254.40067.*SOCₙ)-exp(49.97886.*SOCₙ-43.37888)-0.028221.*atan(22.52300.*SOCₙ-3.65328)-0.01308.*atan(28.34801.*SOCₙ-13.43960)
end

"""Open circuit potential in positive electrode"""
function Uₚ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    DoD = DoD(c₁, c₁ᵣ, I)
    4.2344-9.1296.*DoD.^6+25.8028.*DoD.^5-26.0238.*DoD.^4+11.1602.*DoD.^3-1.9671.*DoD.^2-0.2934.*DoD
end

# Solid potential in different domains of battery
# Refer equation 157.

"""Solid potential in negative electrode"""
function ϕ₁ₙ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I, x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    Uₙ(c₁, c₁ᵣ, I) + ϕ₂(x::Float64,c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux, I) + 2*R*T*asinh(jₙ(I)./(2*jₙ₀(c₁, c₁ᵣ, I, x, c₂ᵢₖ, q₂ᵢₖ)))/F
end

"""Solid potential in positive electrode"""
function ϕ₁ₚ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I, x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    Uₚ(c₁, c₁ᵣ, I) + ϕ₂(x::Float64,c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux, I) 
    + 2*R*T*asinh(jₚ(I)./(2*jₚ₀(c₁, c₁ᵣ, I, x, c₂ᵢₖ, q₂ᵢₖ)))/F
end

# Overall solid potential @x

"""Overall solid potential function"""
function ϕ₁(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I, x::Float64, c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    if x<=lₙ
        ϕ₁ₙ(c₁, c₁ᵣ, I, x, c₂ᵢₖ, q₂ᵢₖ)
    else
        ϕ₁ₚ(c₁, c₁ᵣ, I, x, c₂ᵢₖ, q₂ᵢₖ)
    end
end

