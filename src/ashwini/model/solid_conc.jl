include("reaction_rate.jl")
include("../../utils/rk4.jl")

# Volume averaged and radially averaged solid concentration

"""
Volume averaged solid concentration

### Fields
- c₁ₙ: Solid phase concentration in active material spheres in negative electrode.
- c₁ₚ: Solid phase concentration in active material spheres in positive electrode.
"""
struct cₛ
    c₁ₙ::Vector{Float64}
    c₁ₚ::Vector{Float64}
end

function cₛ(t::Vector{Float64}, tspan, I)
    c₁ₖ₀ = [c₁ₙ₀; c₁ₚ₀]
    c₁ₖ = rk4(dc₁dt, c₁ₖ₀, t, tspan, I)
    c₁ₙ = c₁ₖ[1, :]
    c₁ₚ = c₁ₖ[2, :]

    return cₛ(c₁ₙ, c₁ₚ)
end

function dc₁dt(t, c₁, args...)
    df = zeros(2, 1)
    I = args[1]

    dₙ = args[2]
    dₚ = args[3]

    dₙ₁ = dₙ[1]
    dₙ₂ = dₙ[2]

    dₚ₁ = dₚ[1]
    dₚ₂ = dₚ[2]

    if I isa Float64
        df[1] = -3*jₙ(x, dₙ₁, dₙ₂, I)/rₙ;
        df[2] = -3*jₚ(x, dₚ₁, dₚ₂, I)/rₚ;
    else
        i = args[2]
        df[1] = -3*jₙ(x, dₙ₁, dₙ₂, i)/rₙ;
        df[2] = -3*jₚ(x, dₚ₁, dₚ₂, i)/rₚ;
    end

    return df
end

"""
Radially averaged solid concentration

### Fields
- c₁ₙᵣ: Solid phase concentration in active material spheres in negative electrode.
- c₁ₚᵣ: Solid phase concentration in active material spheres in positive electrode.
"""
struct cₛᵣ
    c₁ₙᵣ::Vector{Float64}
    c₁ₚᵣ::Vector{Float64}
end

function cₛᵣ(t, tspan, I)
    c₁ₖᵣ₀ = Float64.([0.0; 0.0])
    c₁ₖᵣ = rk4(dc₁ᵣdt, c₁ₖᵣ₀, t, tspan, I)
    c₁ₙᵣ = c₁ₖᵣ[1, :]
    c₁ₚᵣ = c₁ₖᵣ[2, :]

    return cₛᵣ(c₁ₙᵣ, c₁ₚᵣ)
end

function dc₁ᵣdt(t, c₁ᵣ, args...)
    df = zeros(2, 1)
    I = args[1]

    dₙ = args[2]
    dₚ = args[3]

    dₙ₁ = dₙ[1]
    dₙ₂ = dₙ[2]

    dₚ₁ = dₚ[1]
    dₚ₂ = dₚ[2]



    if I isa Float64
        df[1] = -45*jₙ(x, dₙ₁, dₙ₂, I)/(2*rₙ^2) - 30*D₁ₙ*c₁ᵣ[1]/rₙ^2
        df[2] = -45*jₚ(x, dₚ₁, dₚ₂, I)/(2*rₚ^2) - 30*D₁ₚ*c₁ᵣ[2]/rₚ^2;
    else
        i = args[2]
        df[1] = -45*jₙ(x, dₙ₁, dₙ₂, i)/(2*rₙ^2) - 30*D₁ₙ*c₁ᵣ[1]/rₙ^2
        df[2] = -45*jₚ(x, dₚ₁, dₚ₂, i)/(2*rₚ^2) - 30*D₁ₚ*c₁ᵣ[2]/rₚ^2;
    end    

    return df
end

# Surface solid phase concentration
# Refer equation 149.

"""Surface solid phase concentration in negative electrode"""
function cₛₙ(x, dₙ₁, dₙ₂, c₁::cₛ, c₁ᵣ::cₛᵣ, I)
    c₁.c₁ₙ .- rₙ * jₙ(x, dₙ₁, dₙ₂, I) / (35 * D₁ₙ) .+ 8 * rₙ * c₁ᵣ.c₁ₙᵣ / 35
end

function cₛₙ(jₙ::Float64, c₁ₙ::Float64, c₁ₙᵣ::Float64)
    c₁ₙ - rₙ*jₙ/(35 * D₁ₙ) + 8*rₙ*c₁ₙᵣ/35
end

"""Surface solid phase concentration in positive electrode"""
function cₛₚ(x, dₚ₁, dₚ₂, c₁::cₛ, c₁ᵣ::cₛᵣ, I)
    c₁.c₁ₚ .- rₚ * jₚ(x, dₚ₁, dₚ₂, I) / (35 * D₁ₚ) .+ 8 * rₚ * c₁ᵣ.c₁ₚᵣ / 35
end

function cₛₚ(jₚ::Float64, c₁ₚ::Float64, c₁ₚᵣ::Float64)
    c₁ₚ - rₚ*jₚ/(35 * D₁ₚ) + 8*rₚ*c₁ₚᵣ/35
end

function c₁ₙ(δt, jₙ, c₁ₙₚ)
    c₁ₙₚ - δt * (3*jₙ/rₙ)
end

function c₁ₚ(δt, jₚ, c₁ₚₚ)
    c₁ₚₚ - δt * (3*jₚ/rₚ)
end

function c₁ₙᵣ(δt, jₙ, c₁ₙᵣₚ)
    c₁ₙᵣₚ - δt * (45*jₙ/(2*rₙ^2) - 30*D₁ₙ*c₁ₙᵣₚ/rₙ^2)
end

function c₁ₚᵣ(δt, jₚ, c₁ₚᵣₚ)
    c₁ₚᵣₚ - δt * (45*jₚ/(2*rₚ^2) - 30*D₁ₚ*c₁ₚᵣₚ/rₚ^2)
end