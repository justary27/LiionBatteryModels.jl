# Universal Constants

"""Universal gas constant"""
const R = 8.314

"""Faraday's constant"""
const F = 96487

# Model Constants

"""System Temperature (in K)"""
const T = 298

"""Bruggeman factor"""
const brug = 1.5

"""Electrloyte transference number"""
const t₊ = 0.363

# Thickness of different regions & cell.

"""Thickness of negative electrode"""
const lₙ = 5e-5

"""Thickness of negative separator"""
const lₛ = 25e-6

"""Thickness of positive electrode"""
const lₚ = 57e-6

"""Total cell thickness"""
const L = lₙ + lₛ + lₚ

# Electrloyte volume fraction

"""Electrloyte volume fraction in negative electrode"""
const ϵ₂ₙ = 0.338

"""Electrloyte volume fraction in separator"""
const ϵ₂ₛ = 0.37

"""Electrloyte volume fraction in positive electrode"""
const ϵ₂ₚ = 0.1979

# Active material fraction

"""Active material fraction in negative electrode"""
const ϵ₁ₙ = 0.5871


"""Active material fraction in separator"""
const ϵ₁ₛ = 0

"""Active material fraction in positive electrode"""
const ϵ₁ₚ = 0.6539

# Electrolyte diffusivity

"""Electrolyte diffusivity of material"""
const D₂ = 7.5e-10

"""Effective electrolyte diffusivity in the negative electrode"""
const D₂ₙ = D₂ * (ϵ₂ₙ ^ brug)

"""Effective electrolyte diffusivity in the separator"""
const D₂ₛ = D₂ * (ϵ₂ₛ ^ brug)

"""Effective electrolyte diffusivity in the positive electrode"""
const D₂ₚ = D₂ * (ϵ₂ₚ ^ brug)

# Solid diffusivity

"""Solid diffusivity in negative electrode"""
const D₁ₙ = 3.9e-14

"""Solid diffusivity in positive electrode"""
const D₁ₚ = 1e-13

# Radii of active material spheres

"""Radii of active material spheres in negative electrode"""
const rₙ = 8e-6

"""Radii of active material sphere in positive electrode"""
const rₚ = 5e-6

# Specific surface area of active material

"""Specific surface area of active material in negative electrode"""
const aₙ = 3 * ϵ₁ₙ / rₙ

"""Specific surface area of active material in positive electrode"""
const aₚ = 3 * ϵ₁ₚ / rₚ

# Surface reaction rate constants

"Surface reaction rate constants in negative electrode"
const kₙ = 2e-6 / F

"Surface reaction rate constants in positive electrode"
const kₚ = 2e-6 / F

# Solid phase concentrations in active material spheres (in mol/m³)

"""Initial solid phase concentrations in negative electrode"""
const c₁ₙ₀ = 17e3

"""Maximum solid phase concentrations in negative electrode"""
const c₁ₙₘₐₓ = 31368

"""Initial solid phase concentrations in positive electrode"""
const c₁ₚ₀ = 22e3

"""Maximum solid phase concentrations in positive electrode"""
const c₁ₚₘₐₓ = 35555

# Electrloyte concentrations (in mol/m³)

"""Initial electrolyte concentration"""
const c₂₀ = 1000

# Parameter constants
# See equations 79 & 80 for reference.

"""Parameter constant at separator interface of the negative electrode"""
const αᵢₙ = - ((lₙ*lₛ*ϵ₂ₙ/(2*D₂ₛ)) + ((lₛ^2)*ϵ₂ₛ/(6*D₂ₛ)) + ((lₙ^2)*ϵ₂ₙ/(3*D₂ₙ))) / (lₙ*ϵ₂ₙ + lₛ*ϵ₂ₛ + lₚ*ϵ₂ₚ)

"""Parameter constant at separator interface of the positive electrode"""
const αᵢₚ = - ((lₙ*lₛ*ϵ₂ₙ/(2*D₂ₛ)) + ((lₛ^2)*ϵ₂ₛ/(3*D₂ₛ)) - ((lₚ^2)*ϵ₂ₚ/(3*D₂ₚ))) / (lₙ*ϵ₂ₙ + lₛ*ϵ₂ₛ + lₚ*ϵ₂ₚ)

const θ = R*T*(1-t₊)/F

const SOCₚₘᵢₙ=0.615617983

const SOCₚₘₐₓ=1

