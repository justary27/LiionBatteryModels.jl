# Universal Constants

"""Universal gas constant"""
const R = 8.314

"""Faraday's constant"""
const F = 96487

# Model Constants

"""System Temperature (in K)"""
const T = 298

"""Transferennce number of Li ion species dissolved in liquid"""
const t₊ = 0.363


# Transfer coefficients

"""Anodic transfer coefficient"""
const αₐ = 0.5

"""Cathodic transfer coefficient"""
const αᵪ = 0.5


# Bruggeman factor in different battery regions

"""Bruggeman factor in negative electrode"""
const brugₙ = 1.5

"""Bruggeman factor in separator"""
const brugₛ = 4.0

"""Bruggeman factor in positive electrode"""
const brugₚ = 1.5

# Volume fraction in different battery regions

"""Volume fraction in negative electrode"""
const ϵ₂ₙ = 0.33

"""Volume fraction in separator"""
const ϵ₂ₛ = 0.54

"""Volume fraction in positive electrode"""
const ϵ₂ₚ = 0.332
        
# Active material fraction in -ve & +ve electrode

"""Active material fraction in negative electrode"""
const ϵ₁ₙ = 0.5

"""Active material fraction in positive electrode"""
const ϵ₁ₚ = 0.49

"""Electronic conductivity in negative electode"""
const k₁ₙ = 100

"""Electronic conductivity in positive electode"""
const k₁ₚ = 3.8


# Thickness of different electrode regions (in m)

"""Thickness of negative electrode"""
const lₙ = 120e-6

"""Thickness of separator"""
const lₛ = 30e-6

"""Thickness of positive electrode"""
const lₚ = 150e-6

"""Total cell thickness"""
const L = lₙ + lₛ + lₚ

# Electrloyte concentrations (in mol/m³)

"""Initial electrolyte concentration"""
const c₂₀ = 1200.00

# Solid phase concentrations in active material spheres (in mol/m³)

"""Initial solid phase concentrations in negative electrode"""
const c₁ₙ₀ = 12e3

"""Initial solid phase concentrations in positive electrode"""
const c₁ₚ₀ = 9e3

# Solid phase concentrations in active material spheres (in mol/m³)

"""Maximum solid phase concentrations in negative electrode"""
const c₁ₙₘₐₓ = 26390.00

"""Maximum solid phase concentrations in positive electrode"""
const c₁ₚₘₐₓ = 22860.00

#  Characteristic radius of electrode particles (in m)
"""Radii of active material spheres in negative electrode"""
const rₙ = 12.5e-6

"""Radii of active material spheres in positive electrode"""
const rₚ = 8.5e-6

# Specific surface area of active material

"""Specific surface area of active material in negative electrode"""
const aₙ = 3 * ϵ₁ₙ / rₙ

"""Specific surface area of active material in positive electrode"""
const aₚ = 3 * ϵ₁ₚ / rₚ

# Miscellaneous
const θ = R*T*(1-t₊)/F