%% Tensao induzida pelo ima
alpha = 117;
Ml = 14e-3;
Mh = 2.3e-3;
Mshgt = 2.5e-3;
Br = 1.15;
r_perc = -0.09;
ur = 1.05;

ut = (Mshgt*ur/Mh)/(1 + ur*(Mshgt - Mh)/Mh);
Br_t = Br - Br*(T - 20)*r_perc;
Br_ts = Br_t/(1 + ut*(1 + Mshgt - Mh)/Mh);
phi_r = A_ima*Br_ts;
f_lkg = 1;
Pmo = Mh/(ut*A_ima);
phi_m = phi_r/(1 + f_lkg*Pmo*Rpol);
phi_l = f_lkg*phi_m;
Ag = pi*RIE*CP*.5/P;
Bg_l = phi_l/Ag;
phi_l_fund = Bg_l*2*RIE*CP/P;
phi = Neep*phi_m;
V = w*phi/sqrt(2);