include("solid_potential.jl")

"""Cell voltage across the battery"""
function V(ϕ₁ₚₗ, ϕ₁ₙ₀)
    ϕ₁ₚₗ - ϕ₁ₙ₀
end
