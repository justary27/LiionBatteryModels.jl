% Include diffusivity.m file
% diffusivity

% param_solver
% The key function to get all the parameters required to solve 
% for c2 at a given timestep in the timeseries.
%
% Returns the 16 parameters: n0(t), n1(t), n2(t), n3(t),
% p0(t), p1(t), p2(t), p3(t), s0(t), s1(t), s2(t), <c2>n,
% <c2>p , <c2>s, c2|x=0 and c2|x=L in this order
%
% Params:
%     - I: The current at the timestep.
%     - Dt: The span of the timestep.
%     - c2np: The previous average concentration in negative electrode
%     - c2sp: The The previous average concentration in separator
%     - c2pp: The previous average concentration in positive electrode
%     - c2n0p: The previous concentration in negative electrode @ x=0
%     - c2plp: The previous average concentration in positive electrode @x=L
%     - jn_0: The reaction rate @x=0
%     - jp_l: The reaction rate @x=L

function X = param_solver(I, Dt, c2np, c2sp, c2pp, c2n0p, c2plp, jn_0, jp_l )

    A = zeros(16, 16);
    B = zeros(16, 1);
    
    % Eqn 73
    A(1, 2) = 1;
    A(2, 6) = 1;
    
    % Eqn 74
    A(3, 1) = -1; A(3, 9) = 1;

    % Eqn 75
    A(4, 1) = -1; A(4, 5) = 1; A(4, 10) = -ls; A(4, 11) = -ls^2;

    % Eqn 76
    A(5, 3) = -2*d2n(c2np)*ln;
    A(5, 4) = -3*d2n(c2np)*ln^2;
    A(5, 10) = -d2s(c2sp);
    A(5, 11) = -2*ln*d2s(c2sp);

    % Eqn 77
    A(6, 7) = 2*d2p(c2pp)*lp;
    A(6, 8) = 3*d2p(c2pp)*lp^2;
    A(6, 10) = -d2s(c2sp);
    A(6, 11) = -2*(ln+ls)*d2s(c2sp);

    % Eqn 78, 79, 80
    A(7, 1) = 1; A(7, 3) = (2*ln^2)/3; A(7, 4) = (3*ln^3)/4; A(7, 12) = -1;
    A(8, 5) = 1; A(8, 7) = (2*lp^2)/3; A(8, 8) = (3*lp^3)/4; A(8, 13) = -1;
    A(9, 9) = 1; A(9, 10) = ls/2; A(9, 11) = (2*ls^2)/3; A(9, 14) = -1;

    % Eqn 81
    A(10, 12) = ln*e2n; A(10, 13) = lp*e2p; A(10, 14) = ls*e2s;

    % Eqn 82 & 83
    A(11, 3) = 2 * d2n(c2np) * Dt/e2n;
    A(11, 4) = 3 * ln * d2n(c2np) * Dt/e2n;
    A(11, 12) = 1;

    A(12, 7) = 2 * d2p(c2pp) * Dt/e2p;
    A(12, 8) = 3 * lp * d2p(c2pp) * Dt/e2p;
    A(12, 13) = 1;

    % Eqn 84 & 85
    A(13, 3) = 2*d2n(c2np)*Dt/e2n; A(13, 15) = 1;
    A(14, 7) = 2*d2p(c2pp)*Dt/e2p; A(14, 16) = 1;

    % Eqn 86 & 87
    A(15, 1) = 1; A(15, 3) = ln^2; A(15, 4) = ln^3; A(15, 15) = -1;
    A(16, 5) = 1; A(16, 7) = lp^2; A(16, 8) = lp^3; A(16, 16) = -1;

    % Filling the B matrix to solve the eqn AX = B
    B(10) = c20 *(ln*e2n + ls*e2s + lp*e2p);
    B(11) = an*(1 - tplus)*jn(I)*Dt/e2n + c2np;
    B(12) = ap*(1 - tplus)*jp(I)*Dt/e2p + c2pp;
    B(13) = an*(1 - tplus)*jn_0*Dt/e2n + c2n0p;
    B(14) = ap*(1 - tplus)*jp_l*Dt/e2p + c2plp;

    X = A\B;
end
