include("reaction_rate.jl")

# Volume averaged and radially averaged solid concentration
# Refer equations 125 & 153

"""
Volume averaged solid concentration

### Fields
- c₁ₙ: Solid phase concentration in active material spheres in negative electrode.
- c₁ₚ: Solid phase concentration in active material spheres in positive electrode.
"""
struct SolidConcentration
    c₁ₙ::Vector{Float64}
    c₁ₚ::Vector{Float64}
end

function SolidConcentration(t::Vector{Float64}, tspan, I)
    c₁ₖ₀ = [c₁ₙ₀; c₁ₚ₀]
    c₁ₖ = rk4(dc₁dt, c₁ₖ₀, t, tspan, I)
    c₁ₙ = c₁ₖ[1, :]
    c₁ₚ = c₁ₖ[2, :]

    return SolidConcentration(c₁ₙ, c₁ₚ)
end

function dc₁dt(t, c₁, args...)
    df = zeros(2, 1)
    I = args[1]


    if I isa Float64
        df[1] = -3*jₙ(I)/rₙ;
        df[2] = -3*jₚ(I)/rₚ;
    else
        i = args[2]
        df[1] = -3*jₙ(i)/rₙ;
        df[2] = -3*jₚ(i)/rₚ;
    end

    return df
end

"""
Radially averaged solid concentration

### Fields
- c₁ₙᵣ: Solid phase concentration in active material spheres in negative electrode.
- c₁ₚᵣ: Solid phase concentration in active material spheres in positive electrode.
"""
struct SolidRadialConcentration
    c₁ₙᵣ::Vector{Float64}
    c₁ₚᵣ::Vector{Float64}
end

function SolidRadialConcentration(t, tspan, I)
    c₁ₖᵣ₀ = Float64.([0.0; 0.0])
    c₁ₖᵣ = rk4(dc₁ᵣdt, c₁ₖᵣ₀, t, tspan, I)
    c₁ₙᵣ = c₁ₖᵣ[1, :]
    c₁ₚᵣ = c₁ₖᵣ[2, :]

    return SolidRadialConcentration(c₁ₙᵣ, c₁ₚᵣ)
end

function dc₁ᵣdt(t, c₁ᵣ, args...)
    df = zeros(2, 1)
    I = args[1]


    if I isa Float64
        df[1] = -45*jₙ(I)/(2*rₙ^2) - 30*D₁ₙ*c₁ᵣ[1]/rₙ^2
        df[2] = -45*jₚ(I)/(2*rₚ^2) - 30*D₁ₚ*c₁ᵣ[2]/rₚ^2;
    else
        i = args[2]
        df[1] = -45*jₙ(i)/(2*rₙ^2) - 30*D₁ₙ*c₁ᵣ[1]/rₙ^2
        df[2] = -45*jₚ(i)/(2*rₚ^2) - 30*D₁ₚ*c₁ᵣ[2]/rₚ^2;
    end    

    return df
end

# Surface solid phase concentration
# Refer equation 149.

"""Surface solid phase concentration in negative electrode"""
function cₛₙ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    c₁.c₁ₙ .- rₙ*jₙ(I)/(35*D₁ₙ) .+ 8*rₙ*c₁ᵣ.c₁ₙᵣ/35
end

"""Surface solid phase concentration in positive electrode"""
function cₛₚ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    c₁.c₁ₚ .- rₚ*jₚ(I)/(35*D₁ₚ) .+ 8*rₚ*c₁ᵣ.c₁ₚᵣ/35
end

