include("params.jl")
"""Electrolyte diffusivity of material"""
function d₂(c₂)
    1e-4 * 10 ^ (-2.2e-4 * c₂ -4.43 * (54/(T-229-0.05*c₂)));
end

# Electrolyte diffusivity in different regions of battery

"""Electrolyte diffusivity in negative electrode"""
function d₂ₙ(c₂ₙ) 
    d₂(c₂ₙ) * (ϵ₂ₙ ^ brugₙ);
end

"""Electrolyte diffusivity in separator"""
function d₂ₛ(c₂ₚ) 
    d₂(c₂ₚ) * (ϵ₂ₛ ^ brugₛ);
end

"""Electrolyte diffusivity in positive electrode"""
function d₂ₚ(c₂ₚ) 
    d₂(c₂ₚ) * (ϵ₂ₚ ^ brugₚ);
end