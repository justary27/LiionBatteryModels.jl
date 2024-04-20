include("solid_conc.jl")
# Over-potential independent rate pre-factor.

"""
Over-potential in negative electrode
"""
function jₙ₀(x, cₛₙ::Float64, params)
    kₙ*((c₁ₙₘₐₓ - cₛₙ)*cₛₙ*c₂(x, params))^0.5
end

function jₙ₀(cₛₙ, c₂)
    kₙ*((c₁ₙₘₐₓ - cₛₙ)*cₛₙ*c₂)^0.5
end

"""
Over-potential in positive electrode
"""
function jₚ₀(x, cₛₚ::Float64, params)
    kₚ*((c₁ₚₘₐₓ - cₛₚ)*cₛₚ*c₂(x, params))^0.5
end

function jₚ₀(cₛₚ, c₂)
    kₚ*((c₁ₚₘₐₓ - cₛₚ)*cₛₚ*c₂)^0.5
end

"""
State of Charge in negative electrode
"""
function SOCₙ(cₛₙ::Float64)
    cₛₙ/c₁ₙₘₐₓ
end

"""
State of Charge in positive electrode
"""
function SOCₚ(cₛₚ::Float64)
    cₛₚ/c₁ₚₘₐₓ
end

# Open circuit potentials

"""
Open circuit potential in negative electrode
"""
function Uₙ(cₛₙ)
    socₙ = SOCₙ(cₛₙ)

    # 0.16 + 1.32*exp(-3*socₙ) + 10*exp(-2000*socₙ)
    -0.16 + 1.32*exp(-3*socₙ)
end


"""
Open circuit potential in positive electrode
"""
function Uₚ(cₛₚ)
    socₚ = SOCₚ(cₛₚ)
    # 8.76 + -21.8*socₚ + 32.3*socₚ^2 + -16.1*socₚ^3
    7.87 + -17.9*socₚ + 27*socₚ^2 + -13.7*socₚ^3
    # -119 + 1259*socₚ  -5250*socₚ^2 + 11426*socₚ^3 + -13720*socₚ^4 + 8630*socₚ^5 + -2224*socₚ^6
    # 4.1983 + 0.0565*tanh(-14.554*socₚ + 8.6094) - 0.0275*(-1.9011 + 1/(0.9984-socₚ).^0.4924) -0.1571*exp(-0.0474*socₚ.^8) + 0.8102*exp(-40*(socₚ-0.1339))
end


# Solid potential in different domains of battery

"""
Solid potential in negative electrode
"""
function ϕ₁ₙ(x, I, cₛₙ::Float64, paramsₙ, paramsₛ, dₙ₁)
    Uₙ(cₛₙ) + ϕ₂ₙ(x, I, paramsₙ, paramsₛ, dₙ₁) + 2*R*T*asinh(jₙ(x, dₙ₁, dₙ₂(dₙ₁), I)./jₙ₀(x, cₛₙ, paramsₙ))/F
end

function ϕ₁ₙ(jₙ::Float64, cₛₙ::Float64, c₂::Float64, ϕ₂ₙ::Float64)
    Uₙ(cₛₙ) + 2*R*T*asinh(jₙ/jₙ₀(cₛₙ, c₂))/F + ϕ₂ₙ
end

"""
Solid potential in positive electrode
"""
function ϕ₁ₚ(x, I, cₛₚ::Float64, paramsₚ, paramsₛ, dₚ₁)
    Uₚ(cₛₚ) + ϕ₂ₚ(x, I, paramsₚ, paramsₛ, dₚ₁) + 2*R*T*asinh(jₚ(x, dₚ₁, dₚ₂(dₚ₁), I)./jₚ₀(x, cₛₚ, paramsₚ))/F
end

function ϕ₁ₚ(jₚ::Float64, cₛₚ::Float64, c₂::Float64, ϕ₂ₚ::Float64)
    Uₚ(cₛₚ) + 2*R*T*asinh(jₚ/jₚ₀(cₛₚ, c₂))/F + ϕ₂ₚ
end

# Overall solid potential @x

"""
Overall solid potential function
"""
function ϕ₁(x, I, cₛₖ::Float64, paramsₖ, paramsₛ, dₖ₁)
    if x<=lₙ
        ϕ₁ₙ(x, I, cₛₖ, paramsₖ, paramsₛ, dₖ₁)
    elseif x>=lₙ + lₛ
        ϕ₁ₚ(x, I, cₛₖ, paramsₖ, paramsₛ, dₖ₁)
    else
        0
    end
end