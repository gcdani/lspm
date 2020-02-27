global f alpha RIE l P wm t urec g_lin lm h1 h2 D Br_mag Bs_mag Kw1
global Kc Eq Xd Xq

alpha = .7820;
g_lin = 1e-3;
t = 1.5e-3;
D = 4e-3;
h1 = 15.904e-3;
h2 = 8.885e-3;
wm = 27.7e-3;
lm = 8.1e-3;
Br_mag = 1.05;
Bs_mag = 1.88;
urec = 1.05;
P = 3;
f = 360;
Ns = 30;
Kw1 = 0.64;
l = 80e-3;
uo = 4*pi*1e-7;


Ag = alpha*RIE*l*pi/P;
Am = wm*l;
Abdg = t*l;

Cphi = Am/Ag;
betaI = urec*Kc*g_lin*Cphi/lm;
n = lm*(h1 + h2)/(4*D*urec*wm);

lambda = (1 + 2*n + (1/betaI))/(2*(Am/Abdg)*(Br_mag/Bs_mag) - 4);
k1 = 4*sin(alpha*pi*.5)/pi;
k1ad = alpha + sin(alpha*pi)/pi;
k1aq = alpha + 0.1 + (sin(0.1*pi) - sin(alpha*pi))/pi;
kalphad = sin(alpha*pi*0.5)/(alpha*pi*.5);
gd = Kc*g_lin/(k1ad - k1*kalphad/(1 + betaI*(1 + 2*n + 4*lambda)));
gq = Kc*g_lin/k1aq;

Bg = Cphi*Br_mag/(1 + betaI*(1 + 2*n + 4*lambda));
Bm = Br_mag*(1 + betaI*(2*n + 4*lambda))/(1 + betaI*(1 + 2*n + 4*lambda));
BM1 = k1*Bg;

phig = Bg*Ag;
phiM1 = phig*k1*2/(pi*alpha);
Dl = phiM1*P/BM1;

Eq = 2*pi*Kw1*Ns*phiM1*f/sqrt(2);
Xd = 6*uo*Dl*f*Kw1*Kw1*Ns*Ns/(P*P*gd);
Xq = 6*uo*Dl*f*Kw1*Kw1*Ns*Ns/(P*P*gq);