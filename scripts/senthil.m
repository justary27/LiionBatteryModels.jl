% CHN 300 Case Study

% Topic
% Reduced order mathematical model for SOCE in 
% Lithium Ion Battery

% Group Members
% Aryan Ranjan
% Ashish Singh
% Pawan Lahoti

% For the model refer:
% "Reduced order model for a lithium ion cell with uniform 
% reaction rate approximation" by V. Senthil Kumar

% Clearing terminal and memory.
clc;
clear;

% Given Constants
csn_max=31368;
csp_max=35555;
R=8.314;
T=298;
brug=1.5;
F=96487;
I_t=13.5;



% tx=0:10:100;
% ix=[100,100,0,0,0,0,-100,-100,0,0,0];
% I_t=ix;
% plot(tx,ix)

ln=5e-5;
ls=25e-6;
lp=57e-6;

L=ln+ls+lp;

e1p=0.6539;
e1n=0.1979;
e2n=0.338;
e2p=0.1979;
e2s=0.37;

d2s=7.5e-10;
t_plus=0.363;
D2 = 7.5e-10;

rp=5e-6;
rn=8e-6;

c1n0=17e3;
c1p0=22e3;

D1n=3.9e-14;
D1p=1e-13;

kn=2e-6/F;
kp=kn;

% Calculated constants

%an
an=3*e1n/rn;
ap=3*e1p/rp;

%jk Eqn 16 17
jn=I_t/(an*F*ln);
jp=-I_t/(ap*F*lp);

D2n = D2*e2n.^brug;
D2p = D2*e2p.^brug;

% k2
k2=@(c2) 1.0793e-2 + 6.4761e-4.*c2 - 5.2245e-7.*(c2).^2 + 1.3605e-10.*(c2).^3 - 1.1724e-14.*(c2).^4;


% Model Solving

% Eqns 79, 80
alpha_in=ain(ln,ls,lp,e2n,e2p,e2s,d2s,D2n);
alpha_ip=aip(ln,ls,lp,e2n,e2p,e2s,d2s,D2p);

% Coeffs for 83, 84
e=(lp.*e2p.*alpha_ip-(lp.^2.*e2p)./(D2p.*3));
b=(ln.*e2n.*alpha_ip+(ln.*ls.*e2n)./(2.*d2s));
a=ln.*e2n.*alpha_in+ls.*ln.*e2n./(2.*d2s)+ln.^2.*e2n./(3.*D2n);
d=lp.*e2p.*alpha_in;

q0=[0,0];
tspan=0:1:1100;
% tspan=tx;
[t,Q] = ode45(@(t,q)odefun(t,q,t_plus,I_t,F,a,b,d,e),tspan,q0);

% % % Q

q2in = Q(:, 1);
q2ip = Q(:, 2);

% Eqn 77, 78
c20=1000;
c2ip=c20+alpha_in.*q2in+alpha_ip.*q2ip;
c2in=c2ip+ls.*(q2in+q2ip)./(2.*d2s);

% 60, 67, 76
c2x0t=c2in+q2in.*ln./(2.*D2n);
c2xlt=c2ip-q2ip.*lp./(2.*D2p);
c2mid=c2in-3.*ls.*q2in./(8.*d2s)-q2ip.*ls./(8.*d2s);

% 55, 63, 70;

% c2 negative
c2_neg=@(x)c2in+(ln.^2-x.^2).*q2in./(2.*ln.*D2n);

% c2 pos
c2_pos=@(x)c2ip-(lp.^2-(L-x).^2).*q2ip./(2.*lp.*D2p);

% c2 sep
c2_sep=@(x)c2in-((x-ln).*q2in)./d2s+((x-ln).^2.*(q2in-q2ip))./(2.*ls.*d2s);

% 105 & 106
c2_sep_avg=integral(c2_sep,ln,ln+ls,ArrayValued=true)./ls;

theta=R.*T.*(1-t_plus)./F;

k2s=k2(c2_sep_avg).*e2s.^brug;

phi2in = @(k2s) 2.*theta.*log(c2in./c2mid)+I_t.*ls./(2.*k2s);

phi2ip = @(k2s) 2.*theta.*log(c2ip./c2mid)-I_t.*ls./(2.*k2s);

% 104,113,116
phi2_sep= @(x) 2.*theta.*log(c2_sep(x)./c2mid)-I_t./k2s.*(x-(ln+ls/2));

c2_neg_avg=integral(c2_neg,0,ln,ArrayValued=true)./ln;
k2n=k2(c2_neg_avg).*e2n.^brug;

phi2_neg=@(x) phi2in(k2s)+2.*theta.*log(c2_neg(x)./c2in)+(I_t*(ln-x))./k2n-(I_t*(ln-x).^2)./(2.*k2n*ln);

c2_pos_avg=integral(c2_pos,ls+ln,L,ArrayValued=true)./lp;
k2p=k2(c2_pos_avg).*e2p.^brug;

phi2_pos=@(x) phi2ip(k2s)+2.*theta.*log(c2_pos(x)./c2ip)-(I_t*(x-ln-ls))./k2p+(I_t*(x-ln-ls).^2)./(2.*k2p*lp);

% see also 114 117 %%%

% solve the coupled ODE 125 and 153
j = [jn, jp];
r = [rn, rp];

c10_bar = [c1n0, c1p0];

[~, c1_bar] = ode45(@(t, c1_bar)odefunc(t, c1_bar, j, r), tspan, c10_bar);

c1r0_bar = [0, 0];
[~, c1r_bar] = ode45(@(t, c1r_bar)odef(t, c1r_bar, j, r, D1n, D1p), tspan, c1r0_bar);

% Eqn 149
c1n_bar = c1_bar(:,1);
c1p_bar = c1_bar(:,2);

c1rn_bar = c1r_bar(:,1);
c1rp_bar = c1r_bar(:,2);

csn = c1n_bar - rn*jn/(35*D1n)+8*rn*c1rn_bar/35;
csp = c1p_bar - rp*jp/(35*D1p)+8*rp*c1rp_bar/35;

% Eqn 156
jn0 = @(x) kn.*(csn_max-csn).^(0.5).*csn.^(0.5).*c2_neg(x).^(0.5);
jp0 = @(x) kp.*(csp_max-csp).^(0.5).*csp.^(0.5).*c2_pos(x).^(0.5);

% Eqn 157

SOCn=csn./csn_max;
SOCp=csp./csp_max;

SOCpmin=0.615617983;
SOCpmax=1;

DoD = (SOCp-SOCpmin)./(SOCpmax-SOCpmin);

Un=0.13966+0.68920.*exp(-49.20361.*SOCn)+0.41903.*exp(-254.40067.*SOCn)-exp(49.97886.*SOCn-43.37888)-0.028221.*atan(22.52300.*SOCn-3.65328)-0.01308.*atan(28.34801.*SOCn-13.43960);
Up=4.2344-9.1296.*DoD.^6+25.8028.*DoD.^5-26.0238.*DoD.^4+11.1602.*DoD.^3-1.9671.*DoD.^2-0.2934.*DoD;

phi1_neg= @(x) Un+phi2_neg(x)+2*R*T.*asinh(jn./(2.*jn0(x)))/F;
phi1_pos= @(x) Up+phi2_pos(x)+2*R*T.*asinh(jp./(2.*jp0(x)))/F;

% eqn 158
V= phi1_pos(L) - phi1_neg(0);

% Plotting Voltage for constant current I_T=13.5
figure(1)
plot(t,V);
grid("on");

title("Variation of Voltage with time for a Lithium ion battery","FontSize",9);
xlabel("Time (in s)");
ylabel("Voltage (in V)");

% %%%%%%%%%%%%%% End of program %%%%%%%%%%%%%%

% Functions Used

function f=ain(ln,ls,lp,e2n,e2p,e2s,d2s,D2n)
f=-((ln.*ls.*e2n)./(2.*d2s)+((ls.*ls).*e2s./(6.*d2s))+ln.^2.*e2n./(3.*D2n))./(ln.*e2n+ls.*e2s+lp.*e2p);

end

function f=aip(ln,ls,lp,e2n,e2p,e2s,d2s,D2p)
f=-((ln.*ls.*e2n)./(2.*d2s)+((ls.*ls).*e2s./(3.*d2s))-lp.^2.*e2p./(3.*D2p))./(ln.*e2n+ls.*e2s+lp.*e2p);
end

function dqdt=odefun(t,q,t_plus,I_t,F,a,b,d,e)
    dqdt=zeros(2,1);
    dqdt(1)=(((-q(1)+(1-t_plus).*(I_t./F))).*e-(q(2)+(1-t_plus).*(I_t./F)).*b)./(a.*e-b.*d); %%
    dqdt(2)=(a.*(q(2)+(1-t_plus).*(I_t./F))-(-q(1)+(1-t_plus).*(I_t./F)).*d)./(a.*e-b.*d);
end

function dc1dt_bar=odefunc(t, c1_bar, j, r)
    dc1dt_bar = zeros(2,1);
    dc1dt_bar(1) = -3.*j(1)./r(1);
    dc1dt_bar(2) = -3.*j(2)./r(2);
end

function dc1rdt_bar = odef(t, c1r_bar, j, r, D1n, D1p)
    dc1rdt_bar = zeros(2,1);
    dc1rdt_bar(1) = -45.*j(1)./(2.*r(1).^2) - 30*D1n*c1r_bar(1)/r(1).^2;
    dc1rdt_bar(2) = -45.*j(2)./(2.*r(2).^2) - 30*D1p*c1r_bar(2)/r(2).^2;
end


