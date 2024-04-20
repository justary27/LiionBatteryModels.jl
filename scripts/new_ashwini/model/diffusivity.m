% Include params.m file
% params

% Electrolyte diffusivity of material
function d2 = d2(c2)
    d2 = 1e-4 * 10^(-2.2e-4*c2 - 4.43*(54/(T-229-0.005*c2)));
end

% Electrolyte diffusivity in different regions of battery

% Electrolyte diffusivity in negative electrode
function d2n = d2n(c2n) 
    d2n = d2(c2n) * (e2n ^ brugn);
end

% Electrolyte diffusivity in separator
function d2s = d2s(c2s) 
    d2s = d2(c2s) * (e2s ^ brugs);
end

% Electrolyte diffusivity in positive electrode
function d2p = d2p(c2p) 
    d2p = d2(c2p) * (e2p ^ brugp);
end
