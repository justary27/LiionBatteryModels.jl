include("electrolyte_potential.jl")

# Local surface reaction rate
# Refer equations 16 & 17.

"""Local surface reaction rate in negative electrode"""
function jₙ(I)
    I/(aₙ*F*lₙ)
end

"""Local surface reaction rate in positive electrode"""
function jₚ(I)
    -I/(aₚ*F*lₚ)
end

