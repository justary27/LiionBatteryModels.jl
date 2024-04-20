% Include conductivity.m
run('conductivity.m')

% Initial value of d_k_1

% d_n_1_0 function
function result = d_n_1_0(I)
    result = F*l_n*I*(1/k_2_n(c_2_0) - 1/sigma_1_n)/(4*R*T);
end

% d_p_1_0 function
function result = d_p_1_0(I)
    result = F*l_p*I*(1/k_2_p(c_2_0) - 1/sigma_1_p)/(4*R*T);
end
