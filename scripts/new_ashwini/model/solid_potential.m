% Over-potential independent rate pre-factor.

% Over-potential in negative electrode
function result = j_n_0(x, c_s_n, params)
    result = k_n*((c_1_n_max - c_s_n)*c_s_n*c_2(x, params))^0.5;
end

function result = j_n_0(c_s_n, c_2)
    result = k_n*((c_1_n_max - c_s_n)*c_s_n*c_2)^0.5;
end

% Over-potential in positive electrode
function result = j_p_0(x, c_s_p, params)
    result = k_p*((c_1_p_max - c_s_p)*c_s_p*c_2(x, params))^0.5;
end

function result = j_p_0(c_s_p, c_2)
    result = k_p*((c_1_p_max - c_s_p)*c_s_p*c_2)^0.5;
end

% State of Charge in negative electrode
function result = SOC_n(c_s_n)
    result = c_s_n/c_1_n_max;
end

% State of Charge in positive electrode
function result = SOC_p(c_s_p)
    result = c_s_p/c_1_p_max;
end

% Open circuit potential in positive electrode
function result = U_p(c_s_p)
    soc_p = SOC_p(c_s_p);

    result = 4.1983 + 0.0565*tanh(-14.554*soc_p + 8.6094) - 0.0275*(-1.9011 + 1/(0.9984-soc_p).^0.4924) -0.1571*exp(-0.0474*soc_p.^8) + 0.8102*exp(-40*(soc_p-0.1339));
end

% Solid potential in different domains of battery

% Solid potential in negative electrode
function result = phi1_n(x, I, c_s_n, params, d_n_1)
    result = U_n(c_s_n) + phi2(x, I, params, d_n_1) + 2*R*T*asinh(j_n(x, d_n_1, d_n_2(d_n_1), I)./j_n_0(x, c_s_n, params))/F;
end

function result = phi1_n(j_n, c_s_n, c_2, phi2_n)
    result = U_n(c_s_n) + 2*R*T*asinh(j_n/j_n_0(c_s_n, c_2))/F + phi2_n;
end

% Solid potential in positive electrode
function result = phi1_p(x, I, c_s_p, params, d_p_1)
    result = U_p(c_s_p) + phi2(x, I, params, d_p_1) + 2*R*T*asinh(j_p(x, d_p_1, d_p_2(d_p_1), I)./j_p_0(x, c_s_p, params))/F;
end

function result = phi1_p(j_p, c_s_p, c_2, phi2_p)
    result = U_p(c_s_p) + 2*R*T*asinh(j_p/j_p_0(c_s_p, c_2))/F + phi2_p;
end

% Overall solid potential @x

% Overall solid potential function
function result = phi1(x, I, c_s_k, params, d_k_1)
    if x<=l_n
        % Continue your code here
    end
end
