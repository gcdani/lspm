%% Programa LSPM - Sincrono
%%% Programa para calculo dos valores de eficiencia,
%%% Angulo de carga do motor estudado
%%% Dados de entrada extraidos do SPEED - 29/Jul/2016


global w raiz f Cap Vlin P beta
global delta Pelec Tps Eff P_cu_tras P_cu_frente P_out T_rel Wcur T_total Tns T_pma;


delta = 0:90;

% P = 1;                  % Par de polos
% f = 50;                 % Frequencia da rede
% beta = 1.7836889285222159;   % Relacao de transformacao entre main e aux.

a_XL = 1;     % Co

%% Entradas
Pwatt_ferro = 3.0345;   % Perda no ferro (falta calcular)
WWF = 0.01;             % Usado para calcular a perda mecanica
Vlin = 220;             % Tensao de entrada (rms)
Cap = 4.75e-6;          % Capacitancia
Ld = a_XL*969.2281e-3;       % Indutancia de eixo direto
Lq = a_XL*5688.9918e-3;      % Indutancia de eixo em quadratura
R1 = 27.7141;           % Resistencia do estator
R2 = 13.5156;           % Resistencia do rotor
E = 154.0069;           % Tensao induzida em aberto
R2d = 49.2121;
X1 = 21.2265;    
R2q = R2d;
X2d = 12.6949;    
X2 = X2d;
Lmd = 914.2586e-3;
Lmq = 5634.0223e-3;

%% Inicio dos calculos
Xmd = 2*pi*f*Lmd;
Xmq = 2*pi*f*Lmq;
V_baixa = Vlin/beta;
E_baixa = E/beta;
R = R1/beta^2;
Xd = 2*pi*f*Ld;
Xq = 2*pi*f*Lq;
Xcap = 2*pi*Cap*f;
ZCap = 1/(1i*Xcap);
Xd_baixa = Xd/beta^2;
Xq_baixa = Xq/beta^2;
% w = 2*pi*f;



Zd = 0.5*1i*Xmd*(R2d/2 + 1i*X2d)/(1i*Xmd + R2d/2 + 1i*X2d);
Zq = 0.5*1i*Xmq*(R2q/2 + 1i*X2)/(1i*Xmq + R2q/2 + 1i*X2);
Z2 = R1 + 1i*X1 + Zd + Zq;


Z2_baixa = Z2/beta^2;

Zmd_menos = 1/(1/(1i*Xmd)+2/(R2d +1i*2*X2d));
Zmq_menos = 1/(1/(1i*Xmq)+2/(R2q +1i*2*X2));
Zdf_menos = Zmd_menos - R/3 + 1i*X1;
Zdb_menos = Zmd_menos - R + 1i*X1;
Zqf_menos = Zmq_menos - R/3 + 1i*X1;
Zqb_menos = Zmq_menos - R + 1i*X1;
Zeq_tras1 = (Zdb_menos - Zqb_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);
Zeq_tras2 =(Zdf_menos + Zqf_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);
Re_neg = 0.5*real(-Zmd_menos -Zmq_menos + ((Zmd_menos - Zmq_menos)^2)/(Zqf_menos + Zdf_menos));

%%
% Matrizes de Transformadas

a1 = [1 1; -1i 1i];
a2 = [1 1; -1i 1i];
a3 = [beta beta;1i -1i];
a4 = [1/beta 1/beta;1i -1i];
a5 = [-1i*beta, 1i*beta;1, 1];
a6 = [-1i/beta, 1i/beta;1, 1];

%%
% Loop

tol = 1e-6;

% Inicializando
len_d = length(delta);
V1 = zeros(len_d,1); V2 = zeros(len_d,1); 
I1 = zeros(len_d,1); I2 = zeros(len_d,1);
Z1 = zeros(len_d,1); Eff = zeros(len_d,1);
Tns = zeros(len_d,1); T_pma = zeros(len_d,1);
T_rel = zeros(len_d,1); Pps = zeros(len_d,1);
P_ferro = zeros(len_d,1); PF = zeros(len_d,1);
Pelec = zeros(len_d,1); WCu1 = zeros(len_d,1);
WCu2 = zeros(len_d,1); WCuR = zeros(len_d,1);
Pshaft = zeros(len_d,1); Tps = zeros(len_d,1);
T_total = zeros(len_d,1); I_sim_nat_baixa = zeros(2,len_d);
V_sim_nat_baixa = zeros(2,len_d); I_sim_baixa = zeros(2,len_d);
V_sim_baixa = zeros(2,len_d); V_alphabeta_baixa = zeros(2,len_d);
I_alphabeta_baixa = zeros(2,len_d); T_total2 = zeros(len_d,1);


for k = 1:length(delta)
    
    flag = 1;
    
    % Calculo das tensoes de seq. positiva e negativa, da impendancia
    % e corrente de seq. positiva.
    while (flag > 0)
        
        V1d = (V_baixa)*sin(delta(k)*pi/180);
        V1q = 1i*(V_baixa)*cos(delta(k)*pi/180);
        V1(k) = -V1d + V1q;                 % Tensao de seq. Positiva
        A = [V1q-1i*E_baixa ; V1d];
        B = [1i*Xd_baixa R; -R -1i*Xq_baixa];
        In = B\A;
        I1(k) = In(1) + In(2);              % Corrente de seq. Positiva
        Z1(k) = V1(k)/I1(k);                % Impedancia de seq. Positiva
        
        a1n = 1 + ZCap/Z1(k);
        a2n = 1 + ZCap/Z2_baixa;
        
        V2(k) = V1(k)*(beta - 1i*a1n)/(beta + 1i*a2n);    % Tensao de seq. -
        
        Vx = -1i*(V1(k) - V2(k));       % Tensao principal
        divisor = abs(Vx)/(Vlin/beta);
        V_baixa = V_baixa/divisor;
        
        P_ferro(k) = Pwatt_ferro;
        
        if (abs(divisor) < 1 + tol)  && (abs(divisor) > 1 - tol)
            flag = 0;
        end
        
    end
    
    % Corrente de seq. Negativa
    I2(k) = V2(k)/(Z2_baixa);               
        
    I_sim_nat_baixa(:,k) = [I1(k);I2(k)];
    V_sim_nat_baixa(:,k) = [V1(k);V2(k)];
    
    V_alphabeta_baixa(:,k) = (1/sqrt(2))*a1*[V1(k);V2(k)];
    I_alphabeta_baixa(:,k) = (1/sqrt(2))*a1*[I1(k);I2(k)];
    
    V_sim_baixa(:,k) = (1/sqrt(2))*a2*V_alphabeta_baixa(:,k);
    I_sim_baixa(:,k) = (1/sqrt(2))*a2*I_alphabeta_baixa(:,k);

    %%
    % Tens?o e corrente principal e auxiliar, e o fator de pot?ncia
    V_m_a = a5*[V1(k); V2(k)];
    I_m_a = a6*[I1(k); I2(k)];
    Is = sum(I_m_a);                            % Corrente de entrada
    PF(k) = cos(phase(V_m_a(1)) - phase(Is));   % Fator de potencia
    Pelec(k) = abs(V_m_a(1)*conj(Is))*PF(k);    % Potencia de entrada

    %Potencia pra frente.
    T_rel(k) = 2*(P/w)*(Xd_baixa-Xq_baixa)*real(I_sim_nat_baixa(1,k))*imag(I_sim_nat_baixa(1,k));
    T_pma(k) = 2*(P/w)*E_baixa*imag(I_sim_nat_baixa(1,k));
    Tps(k) = T_pma(k) + T_rel(k);        % Torque sequencia positiva
    Pps(k) = Tps(k)*w;
    
    % Calculo das perdas
    WCu_main = R1*power(abs(I_m_a(1)),2);
    WCu_aux = R2*power(abs(I_m_a(2)),2);
    WCu1(k) = 2*R*power(abs(I1(k)),2);
    WCu2(k) = 2*R*power(abs(I2(k)),2);
    
    %Torque pra tras.
    T_f_menos = -R*power(abs(-V2(k)*Zeq_tras1),2)/(3*w);
    T_b_menos = Re_neg*power(abs(V2(k)*Zeq_tras2),2)/w;
    Tns(k) = T_f_menos + T_b_menos;
    
    % Perdas no Rotor
    WCuR(k) = abs(Tns(k)*w*2); % o 2 vem de s = 2.
    
    % Calculo das potencias
    Pgap = Pelec(k) - WCu_main - WCu_aux - WCuR(k)*.5; % O 0.5 se deve porque o slip = 1
    Pshaft(k) = Pgap - WCu2(k) - (Pwatt_ferro + WWF); % No SPEED se usa WCuR em vez de WCu2
    
    
    if Pshaft(k) >=0
        Eff(k) = 100*Pshaft(k)/Pelec(k);
    else
        Eff(k) = 100*(1 - (Pshaft(k) - Pelec(k))/Pshaft(k));
    end
    
    if Eff(k) > 100
        Eff(k) = 100;
    elseif Eff(k) < 0
        Eff(k) = 0;
    end
   
    
    T_total2 (k) = 10.2*(Tps(k) + Tns(k)); % 10.2 multiplicado para dar kgf.cm
    T_total (k) = 10.2*Pshaft(k)/w;
        
end

plot(T_total, Eff, 'b:.');
hold on
LSPM_Embraco
plot(Torque_med, Eff_med, 'r');
plot(Speed(:,1), Speed(:,2),'g');
grid on
title('Eficiencia (%)');
xlabel('Torque (kgf.cm)');
ylabel('Eficiencia (%)');


%%%%%%%%%%%%%%%%%%%%  CALCULO DAS PERDAS NO FERRO %%%%%%%%%%%%%%%%%%%%%%
%%%%% W = (Cfch)(f)(B_pk)^Cfa + (CfCe)(2)((pi)(f)(B_pk))^2

kw1 = 0.85;
Tph = 930;
freq = 50;
wm = 2*pi*freq;
Eq1 = 148.2921;
Rad1 = 29.8e-3;
Gap = 0.35e-3;
D = 2*(Rad1 + Gap);
L = 43e-3;

fluxo = Eq1*sqrt(2)/(wm*Tph*kw1);

CfCh = 0.013817;
Cfa = 1.8;
Cfb = 0;
CfCe = 1.6136e-5;

B_pkt = 0.6593;

B_pky = fluxo/(D*L);

Why = (CfCh)*(freq)*(B_pky)^Cfa;
Wey = (CfCe)*(2)*((pi)*(freq)*(B_pky))^2;

Wht = (CfCh)*(freq)*(B_pkt)^Cfa;
Wet = (CfCe)*(2)*((pi)*(freq)*(B_pkt))^2;

%% Insercao MVFL
%%
Wferro = Wht + Wet ;

