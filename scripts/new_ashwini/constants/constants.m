% Universal Constants

% Universal gas constant
R = 8.314;

% Faraday's constant
F = 96487;

% Model Constants

% System Temperature (in K)
T = 298;

% Transferennce number of Li ion species dissolved in liquid
t_plus = 0.363;

% Transfer coefficients

% Anodic transfer coefficient
alpha_a = 0.5;

% Cathodic transfer coefficient
alpha_c = 0.5;

% Bruggeman factor in different battery regions

% Bruggeman factor in negative electrode
brug_n = 1.5;

% Bruggeman factor in separator
brug_s = 4.0;

% Bruggeman factor in positive electrode
brug_p = 1.5;

% Volume fraction in different battery regions

% Volume fraction in negative electrode
epsilon_2n = 0.33;

% Volume fraction in separator
epsilon_2s = 0.54;

% Volume fraction in positive electrode
epsilon_2p = 0.332;

% Active material fraction in -ve & +ve electrode

% Active material fraction in negative electrode
epsilon_1n = 0.5;

% Active material fraction in positive electrode
epsilon_1p = 0.49;

% Electronic conductivity in negative electode
k_1n = 100;

% Electronic conductivity in positive electode
k_1p = 3.8;

% Thickness of different electrode regions (in m)

% Thickness of negative electrode
l_n = 120e-6;

% Thickness of separator
l_s = 30e-6;

% Thickness of positive electrode
l_p = 150e-6;

% Total cell thickness
L = l_n + l_s + l_p;

% Electrloyte concentrations (in mol/m³)

% Initial electrolyte concentration
c_20 = 1200.00;

% Solid phase concentrations in active material spheres (in mol/m³)

% Maximum solid phase concentrations in negative electrode
c_1n_max = 26390.00;

% Maximum solid phase concentrations in positive electrode
c_1p_max = 22860.00;

% Solid phase concentrations in active material spheres (in mol/m³)

% Initial solid phase concentrations in negative electrode
c_1n_0 = 0.74 * c_1n_max;

% Initial solid phase concentrations in positive electrode
c_1p_0 = 0.35 * c_1p_max;

% Solid diffusivity in negative electrode
D_1n = 3.9e-14;

% Solid diffusivity in positive electrode
D_1p = 7.51e-14;

% Characteristic radius of electrode particles (in m)

% Radii of active material spheres in negative electrode
r_n = 12.5e-6;

% Radii of active material spheres in positive electrode
r_p = 8.5e-6;

% Specific surface area of active material

% Specific surface area of active material in negative electrode
a_n = 3 * epsilon_1n / r_n;

% Specific surface area of active material in positive electrode
a_p = 3 * epsilon_1p / r_p;

k_n = 1.764e-11;

k_p = 3.626e-11;

% Miscellaneous
theta = R * T * (1 - t_plus) / F;