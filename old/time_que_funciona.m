%% Programa LSPM - Sincrono
%%% Programa para calculo dos valores de eficiencia,
%%% Angulo de carga do motor estudado
%%% Dados de entrada extraidos do SPEED - 29/Jul/2016


clear; clc;

delta = 0:90;

p = 1;                  % Par de polos
f = 50;                 % Frequencia da rede
beta = 1.7836889285222159;   % Relacao de transformacao entre main e aux.

a_XL = 1;     % usado no SPEED para multiplicar todas as indutancias

%% Entradas
Pwatt_ferro = 3.0345;   % Perda no ferro (falta calcular)
WWF = 0.01;             % Usado para calcular a perda mecanica
Vlin = 220;             % Tensao de entrada (rms)
Cap = 4.75e-6;          % Capacitancia
Ld = 969.2281e-3;       % Indutancia de eixo direto
Lq = 5688.9918e-3;      % Indutancia de eixo em quadratura
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
Xmd = a_XL*2*pi*f*Lmd;
Xmq = a_XL*2*pi*f*Lmq;
V_baixa = Vlin/beta;
E_baixa = E/beta;
R = R1/beta^2;
Xd = 2*pi*f*Ld;
Xq = 2*pi*f*Lq;
ZCap = 1/(1i*2*pi*Cap*f);
ws = 2*pi*f;

Xcap = 2*pi*Cap*f;


negativo_eixo_d = 0.5*1i*Xmd*(R2d/2 + 1i*X2d)/(1i*Xmd + R2d/2 + 1i*X2d);
negativo_eixo_q = 0.5*1i*Xmq*(R2q/2 + 1i*X2)/(1i*Xmq + R2q/2 + 1i*X2);
Z2_alta = R1 + 1i*X1 + 1.0*(negativo_eixo_d + negativo_eixo_q);
Xd_baixa = Xd/beta^2;
Xq_baixa = Xq/beta^2;

Z2_baixa = Z2_alta/beta^2;

Zmd_menos = 1/(1/(1i*Xmd)+2/(R2d +1i*2*X2d));
Zmq_menos = 1/(1/(1i*Xmq)+2/(R2q +1i*2*X2));
Zdf_menos = Zmd_menos -R/3 + 1i*X1;
Zdb_menos = Zmd_menos - R + 1i*X1;
Zqf_menos = Zmq_menos -R/3 + 1i*X1;
Zqb_menos = Zmq_menos -R + 1i*X1;
Zeq_tras1 = (Zdb_menos - Zqb_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);
Zeq_tras2 =(Zdf_menos + Zqf_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);
Re_neg = 0.5*real(-Zmd_menos -Zmq_menos + ((Zmd_menos - Zmq_menos)^2)/(Zqf_menos + Zdf_menos));

% a1 = [1 1; -1i 1i];
% a2 = [1 1i; 1 -1i];
a1 = [1 1; -1i 1i];
a2 = [1 1; -1i 1i];
a3 = [beta beta;1i -1i];
a4 = [1/beta 1/beta;1i -1i];
a5 = [-1i*beta, 1i*beta;1, 1];
a6 = [-1i/beta, 1i/beta;1, 1];


tol = 0.00001;

%%I_1n = zeros(1,lenght(delta)); %Kleyton sugest?o

V1_lin = Vlin;




for k = 1:length(delta)
    
    sentinela = 1;
    
    % Para componentes naturais.
    while (sentinela > 0)
        
        % V_1n(:,k) = (V_baixa)*cos(delta(k)*pi/180)*1i -(V_baixa)*sin(delta(k)*pi/180);
        V1d = (V_baixa)*sin(delta(k)*pi/180);
        V1q = 1i*(V_baixa)*cos(delta(k)*pi/180);
        V_1n(:,k) = -V1d + V1q;
        A_1n = [V1q-1i*E_baixa ; V1d];
        B_1n = [1i*Xd_baixa R; -R -1i*Xq_baixa];
        In = B_1n\A_1n;
        I_1n(:,k) = In(1) + In(2);
        Z_1n(:,k) = V_1n(:,k)/I_1n(:,k);
        
        a1n = 1 + ZCap/Z_1n(:,k);
        a2n = 1 + ZCap/Z2_baixa;
        
        V_2n(:,k) = V_1n(:,k)*(beta - 1i*a1n)/(beta + 1i*a2n);
        
        % Vx = V_1n(:,k) + V_2n(:,k);
        Vx = -1i*(V_1n(:,k) - V_2n(:,k));
        divisor = abs(Vx)/(Vlin/beta);
        V_baixa = V_baixa/divisor;
        
        P_ferro(:,k) = Pwatt_ferro;
        
        if (abs(divisor) < 1 + tol)  && (abs(divisor) > 1 - tol)
            sentinela = 0;
        end
        
    end
    
    I_2n(:,k) = V_2n(:,k)/(Z2_baixa);
        
    I_sim_nat_baixa(:,k) = [I_1n(:,k);I_2n(:,k)];
    V_sim_nat_baixa(:,k) = [V_1n(:,k);V_2n(:,k)];
    
    % Para componentes normais eixos "para frente" e "para tras".
    % V_alphabeta_baixa(:,k) = [1i -1i;1 1]*[V_1n(:,k);V_2n(:,k)];
    % I_alphabeta_baixa(:,k) = [1i -1i;1 1]*[I_1n(:,k);I_2n(:,k)];
    V_alphabeta_baixa(:,k) = (1/sqrt(2))*a1*[V_1n(:,k);V_2n(:,k)];
    I_alphabeta_baixa(:,k) = (1/sqrt(2))*a1*[I_1n(:,k);I_2n(:,k)];
    
    V_sim_baixa(:,k) = (1/sqrt(2))*a2*V_alphabeta_baixa(:,k);
    I_sim_baixa(:,k) = (1/sqrt(2))*a2*I_alphabeta_baixa(:,k);
    
    % Tensao nos enrolamentos auxiliar e principal.
%     V_m_a = a3*V_sim_nat_baixa(:,k);
%     I_m_a = a4*I_sim_nat_baixa(:,k);
    V_m_a = a5*[V_1n(:,k); V_2n(:,k)];
    I_m_a = a6*[I_1n(:,k); I_2n(:,k)];
    Is = sum(I_m_a);
    PF(k) = cos(phase(Is));
    
    %%%%Torques & Pot?ncias
    
    %Potencia de entrada.
    I_in = -1i*(I_1n - I_2n)/beta;
    V_in = -1i*beta*(V_1n - V_2n);
    P_entrada(:,k) = 2*real(conj(I_1n(:,k))*V_1n(:,k)) + 2*real(conj(I_2n(:,k))*V_2n(:,k));
%     PF(k) = cos(phase(V_in) - phase(I_in));
    % P_entrada(:,k) = Vm_rms*(Is)*PF;
%     P_in2(k) = abs(V_m_a(1))*abs(sum(I_m_a))*PF(k);
%     S_in(k) = abs(V_m_a(1))*abs(sum(I_m_a));
    %Potencia pra frente.
    %    T_rel(:,k) = 2*(p/ws)*(Xd_baixa-Xq_baixa)*real(I_sim_nat_baixa(1,k))*imag(I_sim_nat_baixa(1,k));
    %    T_alin(:,k) = 2*(p/ws)*E_baixa*imag(I_sim_nat_baixa(1,k));
    T_rel(:,k) = (p/ws)*(Xd-Xq)*real(In(1))*imag(In(2))/beta;
    T_alin(:,k) = p*imag(In(2))*E/ws;
    T_frente(:,k) = (T_alin(:,k) + T_rel(:,k));
    P_frente(:,k) = T_frente(:,k)*ws;
    
    %Potencia nos enrolamentos.
    P_cu_estator_frente(:,k) = (R*(abs(I_1n(:,k)))^2);
    T_cu_estator_frente(:,k) = P_cu_estator_frente(:,k)/ws;
    P_cu_estator_tras(:,k) = (R*(abs(I_2n(:,k)))^2);
    T_cu_estator_tras(:,k) = P_cu_estator_tras(:,k)/ws;
    
    %Potencia pra tras.
    If_menos(k) = -V_2n(:,k)*Zeq_tras1;
    Ib_menos(k) = V_2n(:,k)*Zeq_tras2;
    
    T_f_menos(k) = (1/(ws))*(R/(-3))*abs(If_menos(k))^2;
    T_b_menos(k) = (1/(ws))*Re_neg*abs(Ib_menos(k))^2;
    
    T_tras_tot(k) = T_f_menos(k) + T_b_menos(k);
    P_tras_tot(k) = ws*T_tras_tot(k);
    T_tras(k) = T_tras_tot(k) - T_cu_estator_tras(:,k);
    P_tras(k) = ws*T_tras(k);
    P_tras(k) = real(conj(I_2n(:,k))*V_2n(:,k)) - P_cu_estator_tras(:,k);
    T_tras(k) = P_tras(k)/ws;
    P1(k) = real(conj(I_1n(:,k))*V_1n(:,k)) - P_cu_estator_frente(:,k);
    T1(k) = P1(k)/ws;
    
    %Potencia no condutor do rotor.
    V_entreferro_tras(:,k) = V_2n(:,k)- ((R1 + 1i*X1)/(beta^2))*(I_2n(:,k));
    V_entreferro_q(:,k) = V_entreferro_tras(:,k)*(negativo_eixo_q/(negativo_eixo_d + negativo_eixo_q));
    V_entreferro_d(:,k) = V_entreferro_tras(:,k)*(negativo_eixo_d/(negativo_eixo_d + negativo_eixo_q));
    I_d_negativo(:,k) = (V_entreferro_d(:,k))/((0.5*(R2d/2 + 1i*X2d ))/beta^2);
    I_q_negativo(:,k) = (V_entreferro_q(:,k))/((0.5*(R2q/2 + 1i*X2))/beta^2);
    Wcur1(:,k) = R2d/(beta^2)*(abs(I_d_negativo(:,k)))^2 ;
    Wcur2(:,k) = R2q/(beta^2)*(abs(I_q_negativo(:,k)))^2;
    Wcur(:,k) = Wcur1(:,k) +  Wcur2(:,k);
    T_Wcur(:,k) = Wcur(:,k)/ws;
    
    P_entreferro(:,k) = P_entrada(:,k) - P_cu_estator_frente(:,k) - P_cu_estator_tras(:,k) - P_ferro(:,k);
    T_entreferro(:,k) = P_entreferro(:,k)/ws;
    P_tras_frente(:,k) = -P_entreferro(:,k) + P_frente(:,k);
    T_tras_frente(:,k) = P_tras_frente(:,k)/(ws);
    
    P1(k) = real(V_1n(:,k)*conj(I_1n(:,k))) - R*abs(power(I_1n(:,k),2));
    P2(k) = real(V_2n(:,k)*conj(I_2n(:,k))) - R*abs(power(I_2n(:,k),2));
    Pout(k) = P_frente(:,k) - P2(k);
    Tmotor(k) = Pout(k)/ws;
    
    WCu_main = R1*power(abs(I_m_a(1)),2);
    WCu_aux = R2*power(abs(I_m_a(2)),2);
    WCu1 = 2*R*power(abs(I_1n(:,k)),2);
    WCu2 = 2*R*power(abs(I_2n(:,k)),2);
    WCuR = abs(T_tras(k)*ws*2); % o 2 vem de s = 2.
    
    Pgap = P_entrada(:,k) - WCu_main - WCu_aux - WCuR*.5; % O 0.5 se deve porque o slip = 1
    Pshaft(k) = Pgap - WCu2 - (Pwatt_ferro + WWF); % No SPEED se usa WCuR em vez de WCu2
    
    if Pshaft(k) < 0
        Pshaft(k) = 0;
    end
    
    
    %Eff(k) = 100*P_entreferro(:,k)/P_entrada(:,k); %100 multiplicado para dar em %
    Eff(k) = 100*Pshaft(k)/P_entrada(:,k);
    Eff_2(k) = 100*Pout(k)/P_entrada(:,k);
    
    if Eff(k) > 100
        Eff(k) = 100;
    elseif Eff(k) < 0
        Eff(k) = 0;
    end
    
    if Eff_2(k) > 100
        Eff_2(k) = 100;
    elseif Eff_2(k) < 0
        Eff_2(k) = 0;
    end
    
    T_frente (k) = T_rel(:,k) + T_alin(:,k);
    
    T_total (k) = 10.2*(-T_tras(:,k) + T_frente(:,k)); % 10.2 multiplicado para dar kgf.cm
    T_total (k) = 10.2*Pshaft(k)/ws;
    
    Tns(k) = p*real(Ld*real(I_2n(k))*imag(I_2n(k)) - Lq*(-imag(I_2n(k)))*real(I_2n(k)));
    
    
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

% %grafico_potent = plot(delta,P_entreferro(:,k));
% subplot(3,3,1);
% grafico_potent = plot(delta, P_entrada, 'b:.');
% %grafico_potent = plot(delta,Eff);
% grid on
% %title('Potencia de entrada (Pot. elec no SPEED)');
% xlabel('Delta (graus)');
% ylabel('Pot. entrada(W)');
%
% subplot(3,3,2);
% grafico_potent = plot(delta,P_cu_estator_frente, 'b:.');
% grid on
% %title('Potencia Cu estator para frente');
% xlabel('Delta (graus)');
% ylabel('Pcu est. frente (W)');
%
% subplot(3,3,3);
% grafico_potent = plot(delta, P_cu_estator_tras, 'b:.');
% grid on
% %title('Potencia Cu estator para tras');
% xlabel('Delta (graus)');
% ylabel('Pcu est. tras (W)');
%
% subplot(3,3,4);
% grafico_potent = plot(delta,P_entreferro, 'b:.');
% grid on
% %title('Pot. entreferro');
% xlabel('Delta (graus)');
% ylabel('Pot. entref. (W)');
%
% subplot(3,3,5);
% grafico_potent = plot(delta,T_alin, 'b:.');
% grid on
% %title('Potencia no condutor do rotor');
% xlabel('Delta (graus)');
% ylabel('Tor. alin. (N.m)');
%
% subplot(3,3,6);
% grafico_potent = plot(delta,T_rel, 'b:.');
% grid on
% %title('Eficiencia');
% xlabel('Delta (graus)');
% ylabel('Torq. rel (N.m)');
%
% subplot(3,3,7);
% grafico_potent = plot(delta,T_total, 'b:.');
% grid on
% %title('Eficiencia');
% xlabel('Delta (graus)');
% ylabel('Torq. total (kgf.cm)');
%
% subplot(3,3,8);
% grafico_potent = plot(delta, Eff, 'b:.');
% grid on
% %title('Eficiencia');
% xlabel('Delta (graus)');
% ylabel('Efici?ncia (%)')
%
% subplot(3,3,9);
% grafico_potent = plot(T_total, Eff, 'b:.');
% grid on
% %title('Efici?ncia (%)');
% xlabel('Torque (kgf.cm)');
% ylabel('Efici?ncia (%)');

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

