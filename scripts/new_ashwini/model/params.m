% Include cell_voltage.m
run('cell_voltage.m')

% The model parameter d_k_1

% TODO: Fix these functions

% d_n_1 function
function result = d_n_1(I, c_s_n_0, c_s_n_l, params)
    result = F*l_n*(I/(2*k_2_n(params)) - I/(2*sigma_1_n) - 2*theta*log(c_2_in(params)./c_2_n_0(params))/l_n - (U_n(c_s_n_l) - U_n(c_s_n_0))/l_n)/(2*R*T) + log(j_n_0(l_n, c_s_n_l, params)./j_n_0(0, c_s_n_0 , params));
end

% d_p_1 function
function result = d_p_1(I, c_s_p_0, c_s_p_l, params)
    result = F*l_p*(I/(2*k_2_p(params)) - I/(2*sigma_1_p) + 2*theta*log(c_2_ip(params)./c_2_p_l(params))/l_p - (U_p(c_s_p_l) - U_p(c_s_p_0))/l_p)/(2*R*T) + log(j_p_0(l_n + l_s, c_s_p_0, params)./j_p_0(L, c_s_p_l, params));
end

% The model parameter d_k_2

% d_n_2 function
function result = d_n_2(d_n_1)
    result = d_n_1/(log((exp(d_n_1) - 1)/d_n_1));
end

% d_p_2 function
function result = d_p_2(d_p_1)
    result = d_p_1/(log((exp(d_p_1) - 1)/d_p_1));
end
