% Include solid_potential.m file
% solid_potential

% Cell voltage across the battery
function V = cell_voltage(I, c1, c1r, paramsn, paramsp, dn1, dp1)
    V = solid_potential(L, I, c1, c1r, paramsp, dn1) - solid_potential(0, I, c1, c1r, paramsn, dp1);
end
