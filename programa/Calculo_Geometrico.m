%% Calculo Geometricos
% Vlin = 198;
% f = 50;
% uo = 4*pi*1e-7;
% pho_cu = 1.72e-8;
% pho_rotor = 2.64e-8;
% g = .35e-3;
% Bsat = 1.1;
% Ns = 24;
% Nr = 28;
% RIE = 30.04e-3;
% CP = 43.5e-3;
% bs = 1.36e-3;
% hp1 = .5e-3;
% hs1 = 11.6e-3;
% bp1 = 6e-3;
% btte = 6.5e-3;
% bttr = 6.5e-3;
% hp2 = .936e-3;
% hs2 = 5.8e-3;
% bs2 = 3.6e-3;
% Sb = 21.5e-6;
% Sa = 80e-6;
% ra = 45e-3;
% hpeq = 0.936e-3;
% MIP = 1;
% MIA = 1;
% NC = 1;
% CI = 1;
% P = 1;
% HRp = 8e-3;
% HRa = 5e-3;
% HC = 3e-3;
% phi_p = 0.5;
% phi_a = 0.515;
% pho_cond = 8890;
% ref = 237.5;
% T = 80;
% w = 2*pi*f;
% Xxm = 0.25;
% hs1a = hs1*.5;
% Nep = [123 168 95 55 24 0];
% Nea = [63 73 81 45 0 0];

global stator beta CP HC P Nep Ns Nea RIE HRa HRp phi_p phi_a pho_cond NC MIA MIP pho_rotor pho_cu Xxm ref;

uo = 4*pi*1e-7;
pho_cu = 1.72e-8;
pho_rotor = 2.64e-8;
MIA = 1;
MIP = 1;
Xxm = 0.25;
ref = 237.5;

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
beta = Neep/Neea;

stator.main_wdg = Neep;
stator.aux_wdg = Neea;

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