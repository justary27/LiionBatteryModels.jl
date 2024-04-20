% Include electrolyte_potential.m
run('electrolyte_potential.m')

% Local surface reaction rate in negative electrode
function result = j_n(x, d_n_1, d_n_2, I)
    result = I*exp(d_n_1*(x/l_n - 1/d_n_2))/(a_n*F*l_n);
end

% Local surface reaction rate in positive electrode
function result = j_p(x, d_p_1, d_p_2, I)
    result = -I*exp(d_p_1*((L-x)/l_p - 1/d_p_2))/(a_p*F*l_p);
end

% Average surface reaction rate in negative electrode
function result = j_n_bar(I)
    result = I/(a_n*F*l_n);
end

% Average surface reaction rate in positive electrode
function result = j_p_bar(I)
    result = -I/(a_p*F*l_p);
end
