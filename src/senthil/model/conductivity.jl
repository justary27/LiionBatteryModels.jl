include("electrolyte_conc.jl")


"""Electrloyte conductivity"""
function k₂(c₂)
    1.0793e-2 .+ 6.4761e-4*c₂ - 5.2245e-7*(c₂).^2 + 1.3605e-10*(c₂).^3 - 1.1724e-14*(c₂).^4;
end

"""Electrloyte conductivity in negative electrode"""
function k₂ₙ(c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux) 
    k₂(c₂̄ₙ(c₂ᵢₖ, q₂ᵢₖ)) * (ϵ₂ₙ^brug)
end

"""Electrloyte conductivity in separator"""
function k₂ₛ(c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)  
    k₂(c₂̄ₛ(c₂ᵢₖ, q₂ᵢₖ)) * (ϵ₂ₛ^brug)
end

"""Electrloyte conductivity in positive electrode"""
function k₂ₚ(c₂ᵢₖ::InterfacialConc, q₂ᵢₖ::InterfacialFlux)
    k₂(c₂̄ₚ(c₂ᵢₖ, q₂ᵢₖ)) * (ϵ₂ₚ^brug)
end

