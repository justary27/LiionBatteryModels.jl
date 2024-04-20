include("reaction_rate.jl")

# Surface solid phase concentration
# Refer equation 149.

"""Surface solid phase concentration in negative electrode"""
function cₛₙ(jₙ::Float64, c₁ₙ::Float64, c₁ₙᵣ::Float64)
    c₁ₙ - rₙ*jₙ/(35 * D₁ₙ) + 8*rₙ*c₁ₙᵣ/35
    # c₁ₙ - rₙ*jₙ/(35 * D₁ₙ) 
end

function cₛₚ(jₚ::Float64, c₁ₚ::Float64, c₁ₚᵣ::Float64)
    c₁ₚ - rₚ*jₚ/(35 * D₁ₚ) + 8*rₚ*c₁ₚᵣ/35
    # c₁ₚ - rₚ*jₚ/(35 * D₁ₚ)
end

"""
Average solid phase concentration in positive electode
"""
function c₁ₙ(δt, jₙ, c₁ₙₚ)
    c₁ₙₚ - δt * (3*jₙ/rₙ)
end

function c₁ₚ(δt, jₚ, c₁ₚₚ)
    c₁ₚₚ - δt * (3*jₚ/rₚ)
end

"""
Gradient of average solid phase concentration in positive electode
"""
function c₁ₙᵣ(δt, jₙ, c₁ₙᵣₚ)
    c₁ₙᵣₚ - δt * (45*jₙ/(2*rₙ^2) - 30*D₁ₙ*c₁ₙᵣₚ/rₙ^2)
end

function c₁ₚᵣ(δt, jₚ, c₁ₚᵣₚ)
    # 3*jₚ(exp(-δt) - 1)/4
    c₁ₚᵣₚ - δt * (45*jₚ/(2*rₚ^2) - 30*D₁ₚ*c₁ₚᵣₚ/rₚ^2)
end