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

# Surface solid phase concentration
# Refer equation 149.

"""Surface solid phase concentration in negative electrode"""
function cₛₙ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    c₁.c₁ₙ - rₙ*jₙ(I)/(35*D₁ₙ) + 8*rₙ*c₁ᵣ.c₁ₙᵣ/35
end

"""Surface solid phase concentration in positive electrode"""
function cₛₚ(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I)
    c₁.c₁ₚ - rₚ*jₚ(I)/(35*D₁ₚ) + 8*rₚ*c₁ᵣ.c₁ₚᵣ/35
end

