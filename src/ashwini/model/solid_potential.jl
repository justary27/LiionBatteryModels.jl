include("solid_conc.jl")
# Over-potential independent rate pre-factor.

"""Over-potential in negative electrode"""
function jₙ₀(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I, x::Float64, params)
    cₛₙ = cₛₙ(c₁, c₁ᵣ, I)
    kₙ*(c₁ₙₘₐₓ .- cₛₙ).^(0.5).*cₛₙ.^(0.5).*c₂(x, params).^(0.5)
end

"""Over-potential in positive electrode"""
function jₚ₀(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I, params)
    cₛₚ = cₛₚ(c₁, c₁ᵣ, I)
    kₚ*(c₁ₚₘₐₓ .- cₛₚ).^(0.5).*cₛₚ.^(0.5).*c₂(x, params).^(0.5)
end

"""State of Charge in negative electrode"""
function SOCₙ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    cₛₙ(c₁, c₁ᵣ, I)/c₁ₙₘₐₓ
end

"""State of Charge in positive electrode"""
function SOCₚ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    cₛₚ(c₁, c₁ᵣ, I)/c₁ₚₘₐₓ

end

# Open circuit potentials

"""Open circuit potential in negative electrode"""
function Uₙ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    socₙ = SOCₙ(c₁, c₁ᵣ, I)

    0.16 + 1.32*exp(-3*socₙ) + 10*exp(-2000*socₙ)
end

"""Open circuit potential in positive electrode"""
function Uₚ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    socₚ = SOCₚ(c₁, c₁ᵣ, I)

    4.1983 + 0.0565*tanh(-14.554*socₚ + 8.6094) - 0.0275*(-1.9011 + 1/(0.9984-socₚ).^0.4924) -0.1571*exp(-0.0474*socₚ.^8) + 0.8102*exp(-40*(socₚ-0.1339))
end


# Solid potential in different domains of battery

"""Solid potential in negative electrode"""
function ϕ₁ₙ(x, I, c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, params)
    Uₙ(c₁, c₁ᵣ, I) + 2*R*T*asinh(jₙ(x, dₙ₁, dₙ₂, I)./jₙ₀(c₁, c₁ᵣ, I, x, params))/F
end

"""Solid potential in positive electrode"""
function ϕ₁ₚ(x, I, c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, params)
    Uₚ(c₁, c₁ᵣ, I) + 2*R*T*asinh(jₚ(x, dₚ₁, dₚ₂, I)./jₚ₀(c₁, c₁ᵣ, I, x, params))/F
end

# Overall solid potential @x

"""Overall solid potential function"""
function ϕ₁(x, I, c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, params)
    if x<=lₙ
        ϕ₁ₙ(x, I, c₁, c₁ᵣ, params)
    elseif x>=lₙ + lₛ
        ϕ₁ₚ(x, I, c₁, c₁ᵣ, params)
    else
        0
    end
end