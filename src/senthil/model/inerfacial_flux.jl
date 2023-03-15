include("../constants/constants.jl")
include("../../utils/rk4.jl")

"""
## Interfacial flux

### Fields
- q₂ᵢₙ: Interfacial flux in negative electrode.
- q₂ᵢₚ: Interfacial flux in positive electrode.
"""
struct InterfacialFlux
    q₂ᵢₙ::Vector{Float64}
    q₂ᵢₚ::Vector{Float64}
end

function InterfacialFlux(t::Vector{Float64}, tspan, I)
    q₂ᵢₖ₀ = Float64.([0.0; 0.0])
    q₂ᵢₖ = rk4(dqdt, q₂ᵢₖ₀, t, tspan, I)
    q₂ᵢₙ = q₂ᵢₖ[1, :]
    q₂ᵢₚ = q₂ᵢₖ[2, :]

    InterfacialFlux(q₂ᵢₙ, q₂ᵢₚ)
end

function dqdt(t, q, I)
    a = lₙ*ϵ₂ₙ*αᵢₙ + lₛ*lₙ*ϵ₂ₙ./(2*D₂ₛ) + lₙ^2*ϵ₂ₙ/(3*D₂ₙ);
    b = (lₙ*ϵ₂ₙ*αᵢₚ + (lₙ*lₛ*ϵ₂ₙ)/(2*D₂ₛ));
    d = lₚ*ϵ₂ₚ*αᵢₙ;
    e = (lₚ*ϵ₂ₚ*αᵢₚ - (lₚ^2*ϵ₂ₚ)/(3*D₂ₚ));

    df = zeros(2, 1)

    df[1] = (((-q[1] + (1-t₊)*(I/F)))*e - (q[2] + (1-t₊)*(I/F))*b)/(a*e - b*d)
    df[2] = (a*(q[2]+(1-t₊)*(I/F)) - (-q[1] + (1-t₊)*(I/F))*d)./(a*e - b*d)

    return df
end