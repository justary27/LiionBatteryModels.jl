"""Over-potential in negative electrode"""
function jₙ₀(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I, x::Float64, params)
    kₙ*(c₁ₙₘₐₓ .- cₛₙ(c₁, c₁ᵣ, I)).^(0.5).*cₛₙ(c₁, c₁ᵣ, I).^(0.5).*c₂(x, params).^(0.5)
end

"""Over-potential in positive electrode"""
function jₚ₀(c₁::SolidConcentration, c₁ᵣ::SolidRadialConcentration, I, x::Float64, params)
    kₚ*(c₁ₚₘₐₓ .- cₛₚ(c₁, c₁ᵣ, I)).^(0.5).*cₛₚ(c₁, c₁ᵣ, I).^(0.5).*c₂(x, params).^(0.5)
end
