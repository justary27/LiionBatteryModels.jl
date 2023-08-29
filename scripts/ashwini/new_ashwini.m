% CHN 400A B.Tech Project

% Development of computationally-efficient physics-based 
% mathematical model for a lithium ion 
% battery management system

% Aryan Ranjan
% 20112026

% For the model refer:
% "A Closed Form Reduced Order Electrochemical Model 
% for Lithium-Ion Cells" by Ashwini Kumar Sharma

% Clearing terminal & memory
clc;
clear;


% Importing system constants & system functions
import Constants.*;

% Aliasing to use import as "c"
c = Constants;


% Model solving starts here...

% Discharge current in Ah
I = 13.5;

% Discharge duration and time step
tspan = 0 : 10 : 25000;
delT = 10;

% @ t = 0, we have dk1's as below
% Refer equations 99 & 100
dn10 = c.F*c.ln*I * (1/c.k2n(c.c20) - 1/c.k1n) /(4*c.R*c.T); % Check the concentration to put here
dp10 = c.F*c.lp*I * (1/c.k2p(c.c20) - 1/c.k1p) /(4*c.R*c.T);

% Refer equations 27 & 28
dn2 = @(dn1) dn1/(log((exp(dn1) - 1)/dn1));
dp2 = @(dp1) dp1/(log((exp(dp1) - 1)/dp1));

% Reaction rate profiles for x = 0 & x = L
% Refer equations 21 & 22
jn_0 = @(dn1) I*exp(-dn1/dn2(dn1))/(c.an*c.F*c.ln);
jp_l = @(dp1) I*exp(-dp1/dp2(dp1))/(c.ap*c.F*c.lp);

% Over-potential independent rate pre-factor
% Refer equation 56
jn0 = @(kn, csn, c2n0) kn * sqrt((c.csn_max - csn).*csn.*c2n0);
jp0 = @(kp, csp, c2pl) kp * sqrt((c.csp_max - csp).*csp.*c2pl);

% Electrolyte concentrations in various battery regions
% Refer equations 68, 70 & 69
c2n = @(n0, n2, n3, x) n0 + n2 * (c.ln^2 - x^2) + n3 * (c.ln^3 - x^3);
c2s = @(s0, s1, s2, x) s0 + s1 * (x - c.ln) + s2 * (x - c.ln)^2;
c2p = @(p0, p2, p3, x) p0 + p2 * (c.lp^2 - (c.L-x)^2) + p3 * (c.lp^3 - (c.L-x)^3);

% Electrolyte potential

% At electrolyte separator interfaces
% Refer equations 36 & 37
phi2in = @(c2in, c2mid, k2s) 2*c.theta*log(c2in./c2mid) + I*c.ls/(2*k2s);
phi2ip = @(c2ip, c2mid, k2s) 2*c.theta*log(c2ip./c2mid) - I*c.ls/(2*k2s);


% At ends of battery
% Refer equations 41 & 43
phi2n0 = @(phi2_in, c2n0, c2in, dn1) phi2_in + 2*c.theta*log(c2n0/c2in) ...
           + I*c.ln*(1-((dn1-1)*exp(dn1)+1)/(dn1*(exp(dn1)-1)))/c.k2n(c2n0);
phi2pl = @(phi2_ip, c2pl, c2ip, dp1) phi2_ip + 2*c.theta*log(c2pl/c2ip) ...
           - I*c.lp*(1-((dp1-1)*exp(dp1)+1)/(dp1*(exp(dp1)-1)))/c.k2p(c2pl);

% Solid potential at battery ends
% Refer equation 94
phi1n0 = @(csn, phi2n_0, c2n0, dn1, kn) c.Un(c.SoCn(csn)) ...
        + phi2n_0 + 2*c.R*c.T*asinh(jn_0(dn1)./jn0(kn, csn, c2n0))/c.F;
phi1pl = @(csp, phi2p_l, c2p0, dp1, kp) c.Up(c.SoCp(csp)) ...
        + phi2p_l + 2*c.R*c.T*asinh(jp_0(dp1)./jp0(kp, csp, c2p0))/c.F;


% Calculating volatage value at different time steps while discharging
dn1_curr = dn10;
dp1_curr = dp10;

% Different concentration time series
c2in = zeros(size(tspan, 2), 1); c2in(1) = c.c20;
c2ip = zeros(size(tspan, 2), 1); c2ip(1) = c.c20;

c2n_avg = zeros(size(tspan, 2), 1); c2n_avg(1) = c.c20;
c2s_avg = zeros(size(tspan, 2), 1); c2s_avg(1) = c.c20;
c2p_avg = zeros(size(tspan, 2), 1); c2p_avg(1) = c.c20;

c2n_0 = zeros(size(tspan, 2), 1); c2n_0(1) = c.c20;
c2p_l = zeros(size(tspan, 2), 1); c2p_l(1) = c.c20;
c2mid = zeros(size(tspan, 2), 1); c2mid(1) = c.c20;

% Different electrolyte potentials
phi2_in = zeros(size(tspan, 2), 1);
phi2_ip = zeros(size(tspan, 2), 1);

phi2n_0 = zeros(size(tspan, 2), 1);
phi2p_l = zeros(size(tspan, 2), 1);

for i = 2 : 3
    % Getting parameters
    params = param_solver( ...
                I, delT, c2n_avg(i-1), c2s_avg(i-1), c2p_avg(i-1), ...
                c2n_0(i-1), c2p_l(i-1), jn_0(dn1_curr), jp_l(dp1_curr) ...
             );
    
    n0 = params(1); n2 = params(3); n3 = params(4);
    p0 = params(5); p2 = params(7); p3 = params(8);
    s0 = params(9); s2 = params(10); s3 = params(11);

    % Calculating interfacial concentrations
    c2in(i) = c2n(n0, n2, n3, c.ln);
    c2ip(i) = c2p(p0, p2, p3, c.ls + c.ln);

    % Calculating concentration at battery ends and center
    c2n_0(i) = c2n(n0, n2, n3, 0);
    c2p_l(i) = c2p(p0, p2, p3, c.L);
    c2mid(i) = c2s(p0, p2, p3, c.ln + c.ls/2);

    % Interfacial Electrolyte potentials
    phi2_in(i) = phi2in(c2in(i), c2mid(i), k2s(c2s_avg(i)));
    phi2_ip(i) = phi2ip(c2ip(i), c2mid(i), k2s(c2s_avg(i)));

    % Eletrolyte potential at battery ends
    phi2n_0(i) = phi2n0(phi2_in(i), c2n_0(i), c2in(i), dn1_curr);
    phi2p_l(i) = phi2pl(phi2_ip(i), c2p_l(i), c2ip(i), dp1_curr);

    % Surface solid phase concentrations

    % Cell Voltage
    % Refer Equation 93
    % V(i) = phi1_l - phi1_0;

    dn1_curr = c.F*c.ln*(I/(2*c.k2n(c.c20)) - I/(2*c.k1n) - ...
               2*c.theta*log(c2in(i-1)/c2cc)/c.ln - c.Un())/(2*c.R*c.T) + ...
               log(jn0(kn, csn, c2n_0(i-1)));
end


% Plotting Voltage
% figure(1)
% plot(t,V);
% grid("on");
% 
% title("Variation of Voltage with time for a Lithium ion battery","FontSize",9);
% xlabel("Time (in s)");
% ylabel("Voltage (in V)");

function params = param_solver( ...
    I, delT, c2n_avg_prev, c2s_avg_prev, c2p_avg_prev, ...
    c2n_0_prev, c2p_l_prev, jn_0, jp_l ...
)

    import Constants.*;
    c = Constants;
    
    A = zeros(16, 16);
    B = zeros(16,1);
    
    % Eqn 73
    A(1, 2) = 1;
    A(2, 6) = 1;
    
    % Eqn 74
    A(3, 1) = -1; A(3, 9) = 1;

    % Eqn 75
    A(4, 1) = -1; A(4, 5) = 1; A(4, 10) = -c.ls; A(4, 11) = -c.ls^2;

    c.D2n(c2n_avg_prev)
    % Eqn 76
    A(5, 3) = -2*c.ln*c.D2n(c2n_avg_prev); 
    A(5, 4) = -3*c.ln^2*c.D2n(c2n_avg_prev);
    A(5, 10) = -c.D2s(c2s_avg_prev);
    A(5, 11) = -2*c.ln*c.D2s(c2s_avg_prev);


    % Eqn 77
    A(6, 7) = 2*c.lp*c.D2p(c2p_avg_prev); 
    A(6, 8) = 3*c.lp^2*c.D2p(c2p_avg_prev);
    A(6, 10) = -c.D2s(c2s_avg_prev);
    A(6, 11) = -2*(c.ln+c.ls)*c.D2s(c2s_avg_prev);

    % Eqn 78, 79, 80
    A(7, 1) = 1; A(7, 3) = 2*c.ln^2/3; A(7, 4) = 3*c.ln^3/4; A(7, 12) = -1;
    A(8, 5) = 1; A(8, 7) = 2*c.lp^2/3; A(8, 8) = 3*c.lp^3/4; A(8, 13) = -1;
    A(9, 9) = 1; A(9, 10) = c.ls/2; A(9, 11) = 2*c.ls^2/3; A(9, 14) = -1;

    % Eqn 81
    A(10, 12) = c.ln*c.e2n; A(10, 13) = c.lp*c.e2p; A(10, 14) = c.ls*c.e2s;

    % Eqn 82 & 83
    A(11, 3) = 2 * c.D2n(c2n_avg_prev); 
    A(11, 4) = 3 * c.ln * c.D2n(c2n_avg_prev);
    A(11, 12) = 1;
    A(11) = A(11)*delT/c.e2n;

    A(12, 7) = 2 * c.D2p(c2p_avg_prev); 
    A(12, 8) = 3 * c.lp * c.D2p(c2p_avg_prev);
    A(12, 13) = 1;
    A(12) = A(12)*delT/c.e2p;

    % Eqn 84 & 85
    A(13, 3) = 2*c.D2n(c2n_avg_prev); A(13, 15) = 1;
    A(13) = A(13)*delT/c.e2n;

    A(14, 7) = 2*c.D2p(c2p_avg_prev); A(14, 16) = 1;
    A(14) = A(14)*delT/c.e2p;

    % Eqn 86 & 87
    A(15, 1) = 1; A(15, 3) = c.ln^2; A(15, 4) = c.ln^3; A(15, 15) = -1;
    A(16, 5) = 1; A(16, 7) = c.lp^2; A(16, 8) = c.lp^3; A(16, 16) = -1;


    % Filling the B matrix to solve the eqn AX = B
    B(10) = c.c20 *(c.ln*c.e2n + c.ls*c.e2s + c.lp*c.e2p);
    B(11) = c.an*(1 - c.t_plus)*c.jn_avg(I)*delT/c.e2n + c2n_avg_prev;
    B(12) = c.an*(1 - c.t_plus)*c.jp_avg(I)*delT/c.e2p + c2p_avg_prev;
    B(13) = c.an*(1 - c.t_plus)*jn_0*delT/c.e2n + c2n_0_prev;
    B(14) = c.an*(1 - c.t_plus)*jp_l*delT/c.e2p + c2p_l_prev;
    
    A
    B
    params = A\B;
end
