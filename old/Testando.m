%% Calculo Geometricos
Vlin = 198;
f = 50;
uo = 4*pi*1e-7;
pho_cu = 1.72e-8;
pho_rotor = 2.64e-8;
g = .35e-3;
Bsat = 1.1;
Ns = 24;
Nr = 28;
RIE = 30.04e-3;
CP = 43.5e-3;
bs = 1.36e-3;
hp1 = .5e-3;
hs1 = 11.6e-3;
bp1 = 6e-3;
btte = 6.5e-3;
bttr = 6.5e-3;
hp2 = .936e-3;
hs2 = 5.8e-3;
bs2 = 3.6e-3;
Sb = 21.5e-6;
Sa = 80e-6;
ra = 45e-3;
hpeq = 0.936e-3;
MIP = 1;
MIA = 1;
NC = 1;
CI = 1;
P = 1;
HRp = 8e-3;
HRa = 5e-3;
HC = 3e-3;
phi_p = 0.5;
phi_a = 0.515;
pho_cond = 8890;
ref = 237.5;
T = 80;
w = 2*pi*f;
Xxm = 0.25;
hs1a = hs1*.5;
Nep = [123 168 95 55 24 0];
Nea = [63 73 81 45 0 0];

Ntotp = 2*P*sum(Nep); % Numero total de ranhuras principal
Ntota = 2*P*sum(Nea); % Numero total de ranhuras auxiliar

for i = 1:length(Nep)
    
    AR(i) = 360*(i-1)/Ns + 360/(2*Ns); % Angulo das ranhuras
    AR2rad(i) = AR(i)*pi/180;
    LCabp(i) = 2*pi*(RIE + HRp)*cos(AR2rad(i))*Nep(i)*NC;
    LCaba(i) = 2*pi*(RIE + HRa)*cos(AR2rad(i))*Nea(i)*NC;
    
end

somap = 0;
somaa = 0;

for i = 1:length(Nep)
    somap = somap + Nep(i)*cos(AR2rad(i))*2;
    somaa = somaa + Nea(i)*cos(AR2rad(i))*2;
end

Neep = somap; % Espiras efetivas principal
Neea = somaa; % Espiras efetivas auxiliar

% Comprimento da cabeca de bobina + colarino
LCtotp = sum(LCabp)*MIP + Ntotp*2*HC; 
LCtota = sum(LCaba)*MIA + Ntota*2*HC; 

Ltotp = LCtotp + Ntotp*CP; % Comprimento do fio total principal
Ltota = LCtota + Ntota*CP; % Comprimento do fio total auxiliar

A_fiop = power(phi_p,2)*pi*.25;
A_fioa = power(phi_a,2)*pi*.25;
Vtotp = Ltotp*A_fiop/1e3; % Volume total de cobre principal
Vtota = Ltota*A_fioa/1e3; % Volume total de cobre auxiliar
Peso_p = pho_cond*Vtotp*2*NC; % Peso total principal

% Nao deu igual
Peso_a = pho_cond*Vtota*2*NC; % Peso total auxiliar

Rb = pho_rotor*CP/Sb; % Resistencia das barras
Rr = pho_rotor*pi*ra/Sa; % Resistencia dos aneis
Rrp = Rr*.5;

soma = 0;

for i = 1:Nr*.5
    soma = soma + power(sin(2*pi*i/Nr),2);
end

Rbp2 = Rb/soma; % Resistencia em paralelo de cada barra considerando o angulo entre elas
R2 = (2*Rbp2 + 2*Rrp)*power(Neep,2); % Resistencia total da gaiola referida pelo estator

R2t = R2*((ref + T)/(ref + 25)); % Corrigida pela temperatura

% Diferente
% Rp = Ltotp*pho_enrol/(A_fiop*power(NC*CI,2)); % Resistencia do enrolamento principal
Rp = Ltotp*2/(57*A_fiop*power(NC*CI,2));
Ra = Ltota*2/(57*A_fioa*power(NC*CI,2));

RTC = 10e-3;
RCMp = 40e-3;
RCMa = 37e-3;

%% Calculo das Indutancias

Kc = 1 - bs*bs*Ns/(pi*(5*g + bs)*RIE); Kc = 1/Kc; % Fator de carter
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
% K1ad = 4*sin(alpha*pi/360)/pi;
% K1 = alpha/180 + sin(alpha*pi/180)/pi;
% K_alpha_d = sin(alpha*pi/360)/(alpha*pi/360);
% A_ima = Ml*Lstk;
% Mperm = uo*A_ima/Mh;
% AirA = pi*RIE*CP*.5/P;
% Rpol = Kc/(AirA*uo);
% Md = K1ad - (K1*K_alpha_d)/(1 + Perm*Rpol);
% Lmo = 4*Sint*Neep*uo/(pi*Kc*P*P);
% Lmd = Lmo*Md; % Inducao sincrona de eixo direto
% 
% Mq = alpha - sin(alpha*pi)/pi;
% Lmq = Lmo*Mq; % Inducao sincrona em quadratura

% Entradas Temporarias
E = 160.3917;   % ???

Ld = 0.8869;    % Falta                     
Lq = 2.0268;    % Falta
% Xls = 19.7167;  % Falta Como Calcular???
Xls = X1d;
Xmd = Xmd/1i;
Xmq = Xmq/1i;
X2d = X2d/1i;
X2q = X2q/1i;
X1d = X1d/1i;

% R1 = 48.9793; 
% R2 = 48.9793; 
% Rp = 27.7141;
% Ra = 13.5156;
% Cap = 9.8448e-6;
Xd = w*Ld;  
Xq = w*Lq;


P_ferro = 3.0345;
delta_carga = 0:(90/2999):90;

% Fim

beta = Neep/Neea;
Xcap = 2*pi*f*Cap;

a1 = [1i -1i; 1 1];
a2 = [1 1i;1 -1i];
a3 = [beta beta; 1i -1i];
a4 = [1/beta 1/beta; 1i -1i];

V_baixa = Vlin/beta;
E_baixa = E/beta;
Xd_baixa = Xd/power(beta,2);
Xq_baixa = Xq/power(beta,2);
R = Rp/power(beta,2);

% Circuito de Calculos
tol = 1e-5;

% Z2 = R2 + X1d + 1i*.5*Xmd*((.5*Ra + X2d)/(1i*Xmd + .5*Ra + X2d)) + ...
%             1i*.5*Xmq*((.5*Ra + X2q)/(1i*Xmq + .5*Ra + X2q));
  
Z2 = Rp + 1i*X1d + 1i*.5*Xmd*((.5*R2 + 1i*X2d)/(1i*Xmd + .5*R2 + 1i*X2d)) + ...
            1i*.5*Xmq*((.5*R2 + 1i*X2q)/(1i*Xmq + .5*R2 + 1i*X2q));

eixo_d = 1i*.5*Xmd*((.5*R2 + 1i*X2d)/(1i*Xmd + .5*R2 + 1i*X2d));
eixo_q = 1i*.5*Xmq*((.5*R2 + 1i*X2q)/(1i*Xmq + .5*R2 + 1i*X2q));
         
Z2_baixa = Z2/power(beta,2);
ZCap = 1/(1i*Xcap);
Zmd_menos = 1/(1/(1i*Xmd)+2/(R2 +1i*2*X2d));
Zmq_menos = 1/(1/(1i*Xmq)+2/(R2 +1i*2*X2q));
Zdf_menos = Zmd_menos -R/3 + 1i*Xls;
Zdb_menos = Zmd_menos - R + 1i*Xls;
Zqf_menos = Zmq_menos -R/3 + 1i*Xls;
Zqb_menos = Zmq_menos -R + 1i*Xls;
Zeq_tras1 = (Zdb_menos - Zqb_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);
Zeq_tras2 =(Zdf_menos + Zqf_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);
Re_neg = 0.5*real(-Zmd_menos -Zmq_menos + ((Zmd_menos - Zmq_menos)^2)/(Zqf_menos + Zdf_menos));


for int = 1:length(delta_carga)

    flag = 1;
    teste = 0;
    
    while(teste==0)
        
        V1n = V_baixa*cos(delta_carga(int)*pi/180)*1i - V_baixa*sin(delta_carga(int)*pi/180);
        
        A1 = [V_baixa*cos(delta_carga(int)*pi/180) - E_baixa;
            V_baixa.*sin(delta_carga(int)*pi/180)];
        B1 = [Xd_baixa, -1i*R;
            -R, -1i*Xq_baixa];
        
        Im = B1\A1;
        I1n = sum(Im);
        Z1n = V1n/I1n;
        
        a1n = 1i + 1i*ZCap/Z1n;
        a2n = -1i - 1i*ZCap/Z2_baixa;
        
        V2n = V1n*((beta - a1n)/(a2n - beta));
        Vs = V1n + V2n;
        
        flag = abs(Vs)*beta/Vlin;
        V_baixa = V_baixa/flag;
        
        if (flag < 1 + tol && flag > 1 - tol)
            teste = 1;
        end
    
    end
    
    I2n = V2n/Z2_baixa;
    V_alpha_beta = (1/sqrt(2)).*a1*[V1n; V2n];
    I_alpha_beta = (1/sqrt(2)).*a1*[I1n; I2n];
    V_1_2 = (1/sqrt(2)).*a2*V_alpha_beta;
    I_1_2 = (1/sqrt(2)).*a2*I_alpha_beta;
    
    V_m_a = a3*[V1n; V2n];
    I_m_a = a4*[I1n; I2n];
    
    P_in(int) = 2*real(V1n*conj(I1n) + V2n*conj(I2n));
    
    T_rel(int) = 2*P*(abs(Xd_baixa) - abs(Xq_baixa))*real(Im(1))*imag(Im(2))/w;
    T_al(int) = 2*P*imag(Im(2))*E_baixa/w;
    T_frente(int) = T_rel(int) + T_al(int);
    % P_cup = R1d*I_m_a(1);
    % P_cua = R2d*I_m_a(2);
    P_cu_frente(int) = 2*R*abs(power(I1n,2));
    P_cu_tras(int) = 2*R*abs(power(I2n,2));
    % P_cu_rotor = 2*V2n*conj(I2n) - P_cu_tras;
    % P_air = P_in - P_cup - P_cua;
    % P_out = P_air - P_cu_rotor - P_ferro;
    
    %Potencia pra tras.
    If_menos(int) = -V2n*Zeq_tras1;
    Ib_menos(int) = V2n*Zeq_tras2;
    
    T_f_menos(int) = 2*(R/(-3))*power(abs(If_menos(int)),2)/w;
    T_b_menos(int) = 2*Re_neg*power(abs(Ib_menos(int)),2)/w; 
    T_tras_tot(int) = T_f_menos(int) + T_b_menos(int);
    P_tras_tot = w*T_tras_tot(int);
    
    % T_tras(int) = T_tras_tot(int) - P_cu_tras(int)/w;
    % P_tras = w*T_tras(int);
    P_tras = real(2*conj(I2n)*V2n) - P_cu_tras(int);
    T_tras(int) = P_tras/w;
     
    %Potencia no condutor do rotor.   
    V_air_tras = V2n - (Rp + 1i*Xls)*(I2n)/power(beta,2);
    V_air_q = V_air_tras*eixo_q/(eixo_d + eixo_q);
    V_air_d = V_air_tras*eixo_d/(eixo_d + eixo_q);
    I_d_negativo = V_air_d*power(beta,2)/(.5*(.5*R2 + 1i*X2d));
    I_q_negativo = V_air_q*power(beta,2)/(.5*(.5*R2 + 1i*X2q));
    Wcur1 = power(abs(I_d_negativo),2)*R2/power(beta,2) ;
    Wcur2 = power(abs(I_q_negativo),2)*R2/power(beta,2);
    Wcur(int) = Wcur1 +  Wcur2;
    T_Wcur = Wcur(int)/w;
           
    P_out(int) = P_in(int) - P_cu_frente(int) - P_cu_tras(int) - P_ferro;
    T_out(int) = P_out(int)/w;
    P_tras_frente(int) = - P_out(int) + T_frente(int)*w;
    T_tras_frente = P_tras_frente(int)/w;
    
    T_total(int) = 10.2*(T_tras(int) + T_frente(int)); % 10.2 multiplicado para dar kgf.cm
    
    Eff(int) = (P_out(int)/P_in(int))*100;
    
    if Eff(int) > 100
        Eff(int) = 1;
    elseif Eff(int) < 0
        Eff(int) = 0;
    end

end

