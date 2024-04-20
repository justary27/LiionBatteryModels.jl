% Include electrolyte_conc.m
run('electrolyte_conc.m')

% Effective electronic conductivity

% Effective electronic conductivity in negative electode
sigma_1_n = k_1_n * (epsilon_1_n ^ brug_n);

% Effective electronic conductivity in positive electode
sigma_1_p = k_1_p * (epsilon_1_p ^ brug_p);

% Ionic conductivity

% Electrloyte conductivity [Ref45 - 25]
function result = k_2(c_2)
    result = 1e-4*c_2*(-10.5 + 0.074*T - 6.69e-5*T^2 + 6.68e-4*c_2 - 1.78e-5*c_2*T + 2.8e-8*c_2*T^2 + 4.94e-7*c_2.^2 - 8.86e-10*T*c_2.^2).^2;
end

% Effective Electrloyte conductivity in negative electrode
function result = k_2_n(params)
    result = k_2(c_2_n_bar(params)) * (epsilon_2_n ^ brug_n);
end

function result = k_2_n(x, params)
    result = k_2(c_2_n(x, params)) * (epsilon_2_n ^ brug_n);
end

function result = k_2_n(c_2_n)
    result = k_2(c_2_n) * (epsilon_2_n ^ brug_n);
end

% Electrloyte conductivity in separator
function result = k_2_s(params)
    result = k_2(c_2_s_bar(params)) * (epsilon_2_s ^ brug_s);
end

function result = k_2_s(x, params)
    result = k_2(c_2_s(x, params)) * (epsilon_2_s ^ brug_s);
end

function result = k_2_s(c_2_s)
    result = k_2(c_2_s) * (epsilon_2_s ^ brug_s);
end

% Electrloyte conductivity in positive electrode
function result = k_2_p(params)
    result = k_2(c_2_p_bar(params)) * (epsilon_2_p ^ brug_p);
end

function result = k_2_p(x, params)
    result = k_2(c_2_p(x, params)) * (epsilon_2_p ^ brug_p);
end

function result = k_2_p(c_2_p)
    result = k_2(c_2_p) * (epsilon_2_p ^ brug_p);
end
