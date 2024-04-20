% Include reaction_rate.m and rk4.m
run('reaction_rate.m')
run('../../utils/rk4.m')


function result = c_s_n(j_n, c_1_n, c_1_n_r)
    result = c_1_n - r_n*j_n/(35 * D_1_n) + 8*r_n*c_1_n_r/35;
end

function result = c_s_p(j_p, c_1_p, c_1_p_r)
    result = c_1_p - r_p*j_p/(35 * D_1_p) + 8*r_p*c_1_p_r/35;
end

% Average solid phase concentration in positive electode
function result = c_1_n(delta_t, j_n, c_1_n_p)
    result = c_1_n_p - delta_t * (3*j_n/r_n);
end

function result = c_1_p(delta_t, j_p, c_1_p_p)
    result = c_1_p_p - delta_t * (3*j_p/r_p);
end

% Gradient of average solid phase concentration in positive electode
function result = c_1_n_r(delta_t, j_n, c_1_n_r_p)
    result = c_1_n_r_p - delta_t * (45*j_n/(2*r_n^2) - 30*D_1_n*c_1_n_r_p/r_n^2);
end

function result = c_1_p_r(delta_t, j_p, c_1_p_r_p)
    result = c_1_p_r_p - delta_t * (45*j_p/(2*r_p^2) - 30*D_1_p*c_1_p_r_p/r_p^2);
end
