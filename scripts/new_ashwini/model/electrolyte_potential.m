% Include initial_params.m
run('initial_params.m')

% initial_params

% Interfacial potential at negative electrode
function result = phi2in(I, params)
    result = 2*theta*log(c2in(params)./c2m(params)) + I*ls/(2*k2s(ln, params));
end

function result = phi2in(I, c2in, c2m)
    result = 2*theta*log(c2in./c2m) + I*ls/(2*k2s(c2in));
end

% Interfacial potential at positive electrode
function result = phi2ip(I, params)
    result = 2*theta*log(c2ip(params)./c2m(params)) - I*ls/(2*k2s(L-lp, params));
end

function result = phi2ip(I, c2ip, c2m)
    result = 2*theta*log(c2ip./c2m) - I*ls/(2*k2s(c2ip));
end

% Electrolyte potentials in different domains of battery

% TODO: Fix below 3 fns

% Electrolyte potential in negative electrode
function result = phi2n(x, I, params, varargin)
    dn1 = varargin{1};

    result = phi2in(I, params) + 2*theta*log(c2(x, params)./c2in(params)) + I*(ln-x)/k2n(x, params) - I*exp(dn1)*((ln-x) - ln*(1-exp(-dn1+dn1*x/ln))/dn1)/(k2n(x, params)*(exp(dn1) - 1));
end

function result = phi2n(x, I, c2, c2in, phi2in, dn1)
    result = phi2in + 2*theta*log(c2/c2in) + I*(ln-x)/k2n(c2) - I*exp(dn1)*((ln-x) - ln*(1-exp(-dn1+dn1*x/ln))/dn1)/(k2n(c2)*(exp(dn1) - 1));
end
