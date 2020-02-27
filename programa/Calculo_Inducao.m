%% Entradas
global t h1 h2 alpha P g l RIE Bs y wm urec lm Br re Kc CP hs1 bp1 Ns
global D f Bpk_dente Bpk_coroa Bg g_lin Xd Xq Eq Kw1

uo = 4*pi*1e-7;


%% Fluxo Entreferro

% Alternativa (nao usa aproximacao, nao considera saturacao)
g_lin = Kc*g;
Ag = alpha*pi*RIE*l/P;
Rg = g_lin/(uo*Ag);
Pg = 1/Rg;
phi_y = Bs*y*l*2;
Am = wm*l;
phi_r = Br*Am;
Pm = urec*uo*Am/lm;
phi_g = (phi_r - phi_y)/(1 + Pm*Rg);
Bg = phi_g/Ag;

% Ag = alpha*RIE*l*pi/P;
% Am = wm*l;
Abdg = t*l;

Cphi = Am/Ag;
betaI = urec*Kc*g_lin*Cphi/lm;
n = lm*(h1 + h2)/(4*D*urec*wm);

lambda = (1 + 2*n + (1/betaI))/(2*(Am/Abdg)*(Br/Bs) - 4);
k1 = 4*sin(alpha*pi*.5)/pi;
k1ad = alpha + sin(alpha*pi)/pi;
k1aq = alpha + 0.1 + (sin(0.1*pi) - sin(alpha*pi))/pi;
kalphad = sin(alpha*pi*0.5)/(alpha*pi*.5);
gd = Kc*g_lin/(k1ad - k1*kalphad/(1 + betaI*(1 + 2*n + 4*lambda)));
gq = Kc*g_lin/k1aq;

Bg = Cphi*Br/(1 + betaI*(1 + 2*n + 4*lambda));
Bm = Br*(1 + betaI*(2*n + 4*lambda))/(1 + betaI*(1 + 2*n + 4*lambda));
BM1 = k1*Bg;

phi_g = Bg*Ag;
phiM1 = phi_g*k1*2/(pi*alpha);
Dl = phiM1*P/BM1;
phi_fun = Bg*2*RIE*CP/P;
phi_concatenado = Neep*phi_fun;

% Tensao Induzida
% Eq = 2*pi*Kw1*N1*phiM1*f/sqrt(2);
Eq = w*phi_concatenado/sqrt(2);

% Reatancias Xd e Xq ( Kw1N1 = Kw1*N1)
Kw1N1 = Eq*sqrt(2)*.5/(f*pi*phiM1);
Xd = 6*uo*Dl*f*power(Kw1N1,2)/(P*P*gd);
Xq = 6*uo*Dl*f*power(Kw1N1,2)/(P*P*gq);

%% Perdas na coroa
Vol_s = pi*CP*(re*re - RIE*RIE);  % Volume do estator
Vol_coroa_s = Vol_s - pi*CP*power(RIE + hs1,2) + 2*pi*CP*RIE*RIE; % Volume da coroa do estator
Bpk_coroa = phi_g/(2*P*CP*(RIE - re - hs1)); % Inducao de pico na coroa

%% Perdas no dente
% Vol_dente = pi*CP*(power(RIE + hs1,2) - RIE*RIE) - 4*CP*Area_ranhura;
Bpk_dente = phi_g*P*2/(CP*bp1*Ns); % Inducao de pico na coroa






























% falta Bmg e Ki




% slot_pitch = pi*D1_in/Ns;
% stator_tooth_width = pi*(D1_in + 2*h14 + b12)/Ns - b12;
% B_stator_teeth = Bmg*slot_pitch/(stator_tooth_width*ki);
% 
% stator_height_tooth = 0.5*b11 + h11 + h12 + 0.5*b12 + h14;
% stator_height_yolk = 0.5*(D1_out - D1_in) - h1t;
