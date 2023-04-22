clc
clear all
% options = optimoptions('fsolve'); 
% options.MaxIterations = 1000;
% options.MaxFunctionEvaluations = 5000;
R=8.314;
T=298;
brug=1.5;
F=96487;
I_t=13.5;
tspan=0:1:3400;

ln=120e-6;
ls=30e-6;
lp=150e-6;
L=ln+lp+ls;

rp=8.5e-6;
rn=12.5e-6;
e1n=0.5;
e1p=0.49;

e2n=0.33;
e2s=0.54;
e2p=0.332;

an=3*e1n/rn;
ap=3*e1p/rp;

D2 = 7.5e-10;
D2n = D2*e2n.^brug;
D2p = D2*e2p.^brug;
D2s = 7.5e-10;

D1n=3.9e-14;
D1p=1e-13;

jn=I_t/(an*F*ln);
jp=-I_t/(ap*F*lp);
t_plus=0.363;

delT=0.01;

thetha=R*T/F*(1-t_plus);

k2=@(c2) 1e-4*c2*(-10.5+0.074*T)^2;

csn_max=26390;
csp_max=22860;

kn=2e-6/F;
kp=kn;

%unknowns
sign=100;
sig1n=sign * ((e1n)^brug);
sigp=3.8;
sig1p=sigp*((e1p)^brug);
k2p=0.0229;
k2n= 0.0511;
initial=[I_t,I_t];

% eqn 99 100 for t=0 only
dn1=((F*ln)./(2*R*T))*(I_t./(2*k2n) - I_t./(2*sig1n));
dp1=((F*lp)./(2*R*T))*(I_t./(2*k2p) - I_t./(2*sig1p));

%%%%%%%%%%% LOOOP START
for i=2:10

%eq 27 and 28
dn2=dn1/log((-1+exp(dn1))/(dn1));
dp2=dp1/log((-1+exp(dp1))/(dp1));

% eqn 21 22
jn_xt=@(x) (I_t/(an*F*ln))*exp(dn1*(x/ln - 1/dn2));
jp_xt=@(x) (I_t/(ap*F*lp))*exp(dp1*((L-x)/lp - 1/dp2));

%eqn 59 and 66
fdn1=@(dn1) ((dn1-1)*exp(dn1)+1)/(dn1*(-1+exp(dn1)));
fdp1=@(dp1) ((dp1-1)*exp(dp1)+1)/(dp1*(-1+exp(dp1)));

% Sub Routines for solving the simentaneous eqn 
% we need to find coeffs of --> eqn 68 69 70

interval=2e5;
c2_avg_p=zeros(interval+2,1);
c2_avg_n=zeros(interval+2,1);
c2_avg_s=zeros(interval+2,1);

c2_x0=zeros(interval+2,1);
c2_xL=zeros(interval+2,1);

%%% diff eqation solve of eqn 82-85
initial_guess=ones(16,1)+1;
lst=fsolve( @(param)set_eqn(param,i,c2_avg_n,c2_avg_p,c2_avg_s,c2_x0,c2_xL,D2n,D2p,D2s,lp,ln,ls,L,e2n,e2s,e2p,an,ap,t_plus,jn,jp,delT,jn_xt,jp_xt),initial_guess);
c2_avg_n(i)=lst(1);
c2_avg_p(i)=lst(2);
c2_avg_s(i)=lst(3);
c2_x0(i)=lst(4);
c2_xL(i)=lst(5);

s0=lst(14);
s1=lst(15);
s2=lst(16);

c2_s=@(x) s0+s1*(x-ln)+s2*(x-ln)^2;
c2_sep_avg=integral(c2_s,ln,ln+ls,ArrayValued=true)./ls; %%
k2s=k2(c2_sep_avg); %%

p0=lst(10);
p1=lst(11);
p2=lst(12);
p3=lst(13);
c2_p=@(x) p0+p1*(lp-L+x)+p2*(lp^2-(L-x)^2)+p3*(lp^3-(L-x)^3);
% c2_pos_avg=integral(c2_p,ls+ln,L,ArrayValued=true)./lp;
% k2p=k2(c2_pos_avg);

n0=lst(6);
n1=lst(7);
n2=lst(8);
n3=lst(9);
c2_n=@(x) n0+n1*(ln-x)+n2*(ln-x^2)+n3*(ln^3-x^3);
% c2_neg_avg=integral(c2_n,0,ln,ArrayValued=true)./ln;
% k2n=k2(c2_neg_avg);

c2ip=c2_s(ln+ls);
c2in=c2_s(ln);
c2mid=c2_s(L/2);

% eqn 36 37

phi2_in = 2*thetha*log(c2in/c2mid)+I_t*ls/(2*k2s);
phi2_ip = 2*thetha*log(c2in/c2mid)-I_t*ls/(2*k2s);
% c2_p(L)

%eqn 41 43
phi2_0=phi2_in +2*thetha*log(c2_n(0)/c2in) + (I_t*ln/k2n)*(1-((dn1-1)*exp(dn1)+1)/(dn1*(-1+exp(dn1))));
phi2_L=phi2_ip +2*thetha*log(c2_p(L)/c2ip) - (I_t*lp/k2p)*(1-((dp1-1)*exp(dp1)+1)/(dp1*(-1+exp(dp1))));

%Solving coupled Odes's 89 and 90
% c1n_bar and c1rp_bar will help us to calc  Uk --> Un , Up
% here four points are x=0, x=ln, x=ls+ln, x=L;
c1n_bar_ig=1000;
c1ln_bar_ig=1000;
c1lp_bar_ig=1000;
c1p_bar_ig=1000;

c1n_bar=c1n_bar_ig+delT*(-3*jn_xt(0))/rn;
c1ln_bar=c1ln_bar_ig+delT*(-3*jn_xt(ln))/rn;
c1lp_bar=c1lp_bar_ig+delT*(-3*jp_xt(ls+ln))/rp;
c1p_bar=c1p_bar_ig+delT*(-3*jp_xt(L))/rp;

c1rn_bar_ig=0;
c1rln_bar_ig=0;
c1rlp_bar_ig=0;
c1rp_bar_ig=0;

c1rn_bar=c1rn_bar_ig+delT*(-45*jn_xt(0)-30*D1n*c1rn_bar_ig)/rn^2;
c1rln_bar=c1rln_bar_ig+delT*(-45*jn_xt(ln)-30*D1n*c1rln_bar_ig)/rn^2;
c1rlp_bar=c1rlp_bar_ig+delT*(-45*jp_xt(ls+ln)-30*D1p*c1rlp_bar_ig)/rp^2;
c1rp_bar=c1rp_bar_ig+delT*(-45*jp_xt(L)-30*D1p*c1rp_bar_ig)/rp^2;

% Un(x=0), Un(x=ln)-->Uln , Up(ls+ln)-->Ulp, Up(ls+ln)

csn = c1n_bar - rn*jn_xt(0)/(35*D1n)+8*rn*c1rn_bar/35;
SOCn=csn./csn_max;
Un=0.13966+0.68920.*exp(-49.20361.*SOCn)+0.41903.*exp(-254.40067.*SOCn)-exp(49.97886.*SOCn-43.37888)-0.028221.*atan(22.52300.*SOCn-3.65328)-0.01308.*atan(28.34801.*SOCn-13.43960);

csln = c1ln_bar - rn*jn_xt(ln)/(35*D1n)+8*rn*c1rln_bar/35;
SOCln=csln./csn_max;
Uln=0.13966+0.68920.*exp(-49.20361.*SOCln)+0.41903.*exp(-254.40067.*SOCln)-exp(49.97886.*SOCln-43.37888)-0.028221.*atan(22.52300.*SOCln-3.65328)-0.01308.*atan(28.34801.*SOCln-13.43960);

cslp = c1lp_bar - rp*jp_xt(ls+ln)/(35*D1p)+8*rp*c1rlp_bar/35;
SOClp=cslp./csp_max;
SOClpmin=0.615617983; % taken from previous paper;
SOClpmax=1;
DoDlp = (SOClp-SOClpmin)./(SOClpmax-SOClpmin);
Ulp=4.2344-9.1296.*DoDlp.^6+25.8028.*DoDlp.^5-26.0238.*DoDlp.^4+11.1602.*DoDlp.^3-1.9671.*DoDlp.^2-0.2934.*DoDlp;

csp = c1p_bar - rp*jp_xt(L)/(35*D1p)+8*rp*c1rp_bar/35;
SOCp=csp./csp_max;
SOCpmin=0.615617983; % taken from previous paper;
SOCpmax=1;
DoD = (SOCp-SOCpmin)./(SOCpmax-SOCpmin);
Up=4.2344-9.1296.*DoD.^6+25.8028.*DoD.^5-26.0238.*DoD.^4+11.1602.*DoD.^3-1.9671.*DoD.^2-0.2934.*DoD;

%eqn 97
jp0 = @(x)kp*sqrt((csp_max-csp)*csp*c2_p(x));
jn0 = @(x)kn*sqrt((csn_max-csn)*csn*c2_n(x));
Voltage(i)=Up-Un+(2*R*T/F)*asinh(jp_xt(L)/jp0(L))-(2*R*T/F)*asinh(jn_xt(0)/jn0(0))-(I_t*ls)/k2s - (I_t*lp*(1-fdp1(dp1)))/k2p -(I_t*ln*(1-fdn1(dn1)))/k2n +2*thetha*log(c2_p(L)/c2_n(0));

% after calc of V we will have to calc d1k eqn 61 67
c2cc_n=c2_n(0);
c2cc_p=c2_p(L);
dn1_nt=((F*ln)/(2*R*T))*((I_t/(2*k2n)) - (I_t/(2*sig1n)) - (2*thetha*log(c2in/c2cc_n)) - (Uln-Un)/ln + log((jn0(ln))/(jn0(0))) );
dn1_pt=((F*lp)/(2*R*T))*((I_t/(2*k2p)) - (I_t/(2*sig1p)) - (2*thetha*log(c2ip/c2cc_p)) - (Ulp-Up)/lp + log((jp0(ln+ls))/(jp0(L))) );

% we have to update dk1 again by eqn 61 and 67
dn1=dn1_nt;
dp1=dn1_pt;
end
plot(Voltage)
fplot(@(x) exp(-(x+5)^2))
%we have to substitute x and y with c2_avg_n(i) and c2_avg_p(i).

function F=set_eqn(param,i,c2_avg_n,c2_avg_p,c2_avg_s,c2_x0,c2_xL,D2n,D2p,D2s,lp,ln,ls,L,e2n,e2s,e2p,an,ap,t_plus,jn,jp,delT,jn_xt,jp_xt)
x=param(1);
y=param(2);
z=param(3);
p=param(4);
q=param(5);
n0=param(6);
n1=param(7);
n2=param(8);
n3=param(9);
p0=param(10);
p1=param(11);
p2=param(12);
p3=param(13);
s0=param(14);
s1=param(15);
s2=param(16);
F(1)= n1;
F(2)= p1;
F(3)= s0-n0;
F(4)= s2*ls^2+s1*ls+n0-p0;
F(5)= D2n*(-2*n2*ln-3*n3*ln^2)+D2s*(-2*s2*ln-s1);
F(6)= D2p*(2*p2*lp+3*p3*lp^2)-D2s*(2*s2*(ln+ls)+s1);
F(7)= x-n0-(2*n2*ln^2)/3-(3*n3*ln^3)/4;
F(8)= y-p0-(2*p2*lp^2)/3-(3*p3*lp^3)/4;
F(9)= z-s0-(2*s2*ls^2)/3-(s1*ls)/2;
F(10)=ln*e2n*x+ls*e2s*z+lp*e2p*y-(ln*e2n+ls*e2s+p*e2p)*1000;%just took c20=1
F(11)=(-D2n*(2*n2+3*n3*ln)+an*(1-t_plus)*jn)/e2n-(x-c2_avg_n(i-1))/delT;
F(12)=(-D2p*(2*p2+3*p3*lp)+ap*(1-t_plus)*jp)/e2p-(y-c2_avg_p(i-1))/delT;
F(13)=(-2*D2n*n2+an*(1-t_plus)*jn_xt(0))/e2n -(p-c2_x0(i-1))/delT;
F(14)=(-2*D2p*p2+ap*(1-t_plus)*jp_xt(L))/e2p -(q-c2_xL(i-1))/delT;
F(15)=-p+n2*ln^2+n3*ln^3+n0;
F(16)=-q+p2*lp^2+p3*lp^3+p0;
end

function dc1n_dt=odefuncn(t,c1_bar,jn,rn,D1n)
    dc1n_dt=zeros(2,1);
    dc1n_dt(1)=-(3*jn)/rn;
    dc1n_dt(2)=-(45*jn)/(2*(rn)^2)-(30*D1n*c1_bar(2))/rn^2;
end

function dc1p_dt=odefuncp(t,c1_bar,jp,rp,D1p)
    dc1p_dt=zeros(2,1);
    dc1p_dt(1)=-(3*jp)/rp;
    dc1p_dt(2)=-(45*jp)/(2*(rp)^2)-(30*D1p*c1_bar(2))/rp^2;
end










