% Include constants.m
run('../constants/constants.m')

% Electrloyte concentration in negative electrode [68]
function result = c_2_n(x, params)
    n_0 = params(1); n_2 = params(2); n_3 = params(3);

    result = n_0 + n_2 * (l_n^2 - x^2) + n_3 * (l_n^3 - x^3);
end

% Electrloyte concentration in separator [70]
function result = c_2_s(x, params)
    s_0 = params(1); s_1 = params(2); s_2 = params(3);

    result = s_0 + s_1 * (x - l_n) + s_2 * (x - l_n)^2;
end

% Electrloyte concentration in positive electrode [69]
function result = c_2_p(x, params)
    p_0 = params(1); p_2 = params(2); p_3 = params(3);

    result = p_0 + p_2 * (l_p^2 - (L-x)^2) + p_3 * (l_p^3 - (L-x)^3);
end

% Overall electrolyte concentration @ x
function result = c_2(x, params)
    if x <= l_n
        result = c_2_n(x, params);
    elseif x > l_n && x <= l_n + l_s
        result = c_2_s(x, params);
    else
        result = c_2_p(x, params);
    end
end

% Average electrolyte concentration in negative electrode
function result = c_2_n_bar(params)
    result = integral(@(x)c_2(x, params), 0, l_n) / l_n;
end

% Average electrolyte concentration in separator
function result = c_2_s_bar(params)
    result = integral(@(x)c_2(x, params), l_n, l_s + l_n) / l_s;
end

% Average electrolyte concentration in positive electrode
function result = c_2_p_bar(params)
    result = integral(@(x)c_2(x, params), l_s + l_n, L) / l_p;
end

% Electrolyte concentration @ x=0 in negative electrode
function result = c_2_n_0(params)
    result = c_2(0, params);
end

% Interfacial electrolyte concentration in negative electrode
function result = c_2_i_n(params)
    result = c_2(l_n, params);
end

% Electrolyte concentration in the middle of battery
function result = c_2_m(params)
    result = c_2(l_n + l_s/2, params);
end

% Interfacial electrolyte concentration in positive electrode
function result = c_2_i_p(params)
    result = c_2(l_n + l_s, params);
end

% Electrolyte concentration @ x = L in positive electrode.
function result = c_2_p_L(params)
    result = c_2(L, params);
end
