include("electrolyte_potential.jl")

"""Local surface reaction rate in negative electrode"""
function jₙ(x, dₙ₁, dₙ₂, I)
    I*exp(dₙ₁*(x/lₙ - 1/dₙ₂))/(aₙ*F*lₙ)
end

"""Local surface reaction rate in positive electrode"""
function jₚ(x, dₚ₁, dₚ₂,I)
    -I*exp(dₚ₁*((L-x)/lₚ - 1/dₚ₂))/(aₚ*F*lₚ)
end

"""Average surface reaction rate in negative electrode """
function j̅ₙ(I)
    I/(aₙ*F*lₙ)
end

"""Average surface reaction rate in positive electrode """
function j̅ₚ(I)
    -I/(aₚ*F*lₚ)
end
