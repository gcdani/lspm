RTC = 10e-3;
RCMp = 40e-3;
RCMa = 37e-3;

%% Calculo das Indutancias
global Kc Kw1 CP Sb Nr g btte bttr hs2 bs2 hpeq bs P Nep Ns Nea RIE hs1a hp1 bp1 hs1 ref Xxm uo w; 

Kc = 1 - bs*bs*Ns/(pi*(5*g + bs)*RIE); Kc = 1/Kc; % Fator de carter

% Caluclo de Kw1
pole_pitch = pi*RIE/(2*P);
q1 = Ns/(2*P);

% If full-pitch coil
wc = pole_pitch;

Kp1 = sin(pi*.5*wc/pole_pitch);
Kd1 = sin(pi*.5)/(q1*sin(pi*.5/q1));
Kw1 = Kd1*Kp1;

% Retornando as Indutancias
Lmp = Neep*Neep*uo*pi*RIE*CP*Xxm*.5/(Kc*g); % Indutancia magnetizacao assincrona
Lma = Neea*Neea*uo*pi*RIE*CP*Xxm*.5/(Kc*g); % Indutancia magnetizacao assincrona

% Nao existe no manual
Xmp = Lmp*2*pi*f;
Xma = Lma*2*pi*f*power(Neep/Neea,2);



Lslotpp = sum(Nep.^2)*uo*CP*hp1*4*P/bs; % Indutancia de dispersao de pescoco de ranhura
Xslotpp = w*Lslotpp;
Lslotrp = sum(Nep.^2)*uo*CP*hs1*2*P/bp1; % Indutancia de dispersao de corpo de ranhura
Xslotrp = w*Lslotrp;

hs1a = hs1*.5;
Lslotpa = sum(Nea.^2)*uo*CP*hp1*4*P/bs; % Indutancia de dispersao de pescoco de ranhura
Xslotpa = w*Lslotpa*power(Neep/Neea,2);
Lslotra = sum(Nea.^2)*uo*CP*hs1a*2*P/bp1; % Indutancia de dispersao de corpo de ranhura
Xslotra = w*Lslotra*power(Neep/Neea,2);

Lcabp = power(Neep,2)*RCMp*RCMp*uo*.5/RTC; % Indutancia de cabeca de bobina do estator
Xcabp = w*Lcabp;

Lcaba = power(Neea,2)*RCMa*RCMa*uo*.5/RTC; % Indutancia de cabeca de bobina do estator
Xcaba = w*Lcaba*power(Neep/Neea,2);

Lzzp = Ntotp*Ntotp*uo*btte*CP*.5/(g*Ns); % Indutancia zig-zag no estator
Xzzp = w*Lzzp;
Lzza = Ntota*Ntota*uo*btte*CP*.5/(g*Ns); % Indutancia zig-zag no estator
Xzza = w*Lzza*power(Neep/Neea,2);

Lzzr = uo*bttr*CP*.5*power(Neep,2)/(g*Nr); % Indutancia zig-zag no rotor
Xzzr = w*Lzzr;

Lslotr = uo*hs2*CP*2/(bs2*Nr); % Indutancia de ranhura do rotor
Xslotr = w*Lslotr*power(Neep,2);

Xeq = Xslotrp + Xslotpp + Xcabp + Xzzp + Xzzr + Xslotr;

Xrr = 0;

for idx = 1:20

    IBT = Vlin*(Neep + Neea)/(Xeq + Rp + R2 + Xrr); % IBT = 2.145e3;

    IBM = IBT*2/Nr;
    delta_theta = 2*pi/Nr;
    
    soma = 0;
    
    for i = 0:(.5*Nr - 1)
        soma = soma + sin(delta_theta*i + .5*delta_theta);
    end
    
    ABT = Sb*soma;
    JB = IBT/ABT;
    fluxo_sat = Bsat*hpeq*CP;
    Lrr = fluxo_sat*2/(IBM*Nr); % Indutancia de pescoco do rotor
    Xrr = w*Lrr*Neep*Neep;

end

%% Nao tem no manual
Vq = Vlin*Neep/(1i*Neea);
TT = 60; 
Zn = 869.3; % ????
Cap = 61e-6;
R2ini = R2t;
Rot_sinc = 60*f; % Rotacao sincrona
Xcap = (1/(w*Cap*1i))*power(Neep/Neea,2);
Rpd = 26.97*Zn; % ????
Rpq = 28.14*Zn; % ????

R1d = Rp*(ref + TT)/(ref + 25);
R1q = Ra*power(Neep/Neea,2)*(ref + TT)/(ref + 25);
X1d = (Xslotpp + Xslotrp + Xcabp + Xzzp)*1i;
X1q = (Xslotpa + Xslotra + Xcaba + Xzza)*1i;

X2d = (Xslotr + Xzzr + Xrr)*1i;
X2q = X2d;

Xmd = Xmp*1i;
Xmq = Xma*1i;
Zmd = 1/((1/Xmd) + (1/Rpd));
Zmq = 1/((1/Xmq) + (1/Rpq));
Z1d = R1d + X1d;
Z1q = R1q + X1q + Xcap;

int = 0;

for Rot_assinc = 0:Rot_sinc
    
    int = int + 1;
    s = 1 - Rot_assinc/Rot_sinc;
    Rot_assinc(int) = Rot_sinc*(1 - s);
    R2F(int) = R2ini/s;
    R2T(int) = R2ini/(2 - s);
    Z2d(int) = R2T(int) + X2d;
    Z2q(int) = R2T(int) + X2q;
    RR(int) = .5*(R2F(int) - R2T(int));

end

%% Precisa de mais dados para calcular Ld e Lq
% alpha = 117;
% Ml = 14e-3;
% Mh = 2.3e-3;
% Mshgt = 2.5e-3;
% Br = 1.15;
% r_perc = -0.09;
% ur = 1.05;
% 
% Lstk = 7e-3; % Faltou...
% 
% K1ad = 4*sin(alpha*pi/360)/pi;
% K1 = alpha/180 + sin(alpha*pi/180)/pi;
% K_alpha_d = sin(alpha*pi/360)/(alpha*pi/360);
% A_ima = Ml*Lstk;
% Mperm = uo*A_ima/Mh;
% AirA = pi*RIE*CP*.5/P;
% Rpol = Kc/(AirA*uo);
% Md = K1ad - (K1*K_alpha_d)/(1 + Mperm*Rpol);
% Lmo = 4*Ml*Mh*Neep*uo/(pi*Kc*P*P);
% Lmd = Lmo*Md; % Inducao sincrona de eixo direto
% 
% Mq = alpha - sin(alpha*pi)/pi;
% Lmq = Lmo*Mq; % Inducao sincrona em quadratura