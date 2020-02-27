clear; clear all; clc;
%%% Programa para calculo dos valores de eficiencia,
%%% Angulo de carga do motor estudado
%%% Dados de entrada extraidos do SPEED - 29/Jul/2016

%% Parte  Sincrona
% global w raiz f Cap Vlin P
% global delta P_in T_frente Eff P_cu_tras P_cu_frente P_out T_rel Wcur T_total T_tras;

delta = 0:90;

%%%% Valor hipotetico para angulo de carga para comparar resultados do teste.
P = 1;
f = 50;                        %input('Entre com a frequencia da rede: ')

%%%% Valores do circuito equivalente tirados do SPEED.

%%%% Fatores de ajustes do SPEED
a_XL = 0.91;     % usado no SPEED para multiplicar todas as indut?ncias

%%%% Valores do circuito equivalente tirados do SPEED.
beta = 1.7836889285222159;
Pwatt_ferro = 3.0345;
Vm_rms = 220;                         %input('Entre com o valor de tensao da fonte de alimenta?ao: ');
V_baixa = Vm_rms/beta;
Cap = 11.45e-6;                        %input('Entre com o valor do capacitor:');
Ld = a_XL*807.1198e-3;                %input('Entre com o valor da indutancia de eixo direto:');
Lq = a_XL*1331.7623e-3;               %input('Entre com o valor da indutancia de eixo em quadratura:');
R1 = 27.7141;                         %input('Entre com o valor da resistencia do estator para o lado de alta'):
R2 = 13.5156;
E = 154.0069;                         %input('Entre com o valor da tensao de enrolamento direto em circuito aberto'):
Rrd = 49.2121;          % R2q = R2d do SPEED
Xls = a_XL*21.0763;     % X1 do SPEED
Rrq = Rrd;
Xlrd = a_XL*12.6949;    % X2 do SPEED
Xlrq = Xlrd;
Lmd = a_XL*757.0975e-3;
Lmq = a_XL*1281.7401e-3;
Xmd = 2*pi*f*Lmd;
Xmq = 2*pi*f*Lmq;
E_baixa = E/beta;
R = R1/beta^2;
Xd = 2*pi*f*Ld;
Xq = 2*pi*f*Lq;
ZCap = 1/(1i*2*pi*Cap*f);
w = 2*pi*f;

Xcap = 2*pi*Cap*f;

% beta = Neep/Neea;
a1 = [1 1; -1i 1i];
a2 = [1 1; 1i -1i];
a3 = [beta beta;1i -1i];
a4 = [1/beta 1/beta;1i -1i];

E_baixa = E/beta;
R = R1/beta^2;

negativo_eixo_d = 0.5*1i*Xmd*(Rrd/2 + 1i*Xlrd)/(1i*Xmd + Rrd/2 + 1i*Xlrd);
negativo_eixo_q = 0.5*1i*Xmq*(Rrq/2 + 1i*Xlrq)/(1i*Xmq + Rrq/2 + 1i*Xlrq);
Z2 = R1 + 1i*Xls + (negativo_eixo_d + negativo_eixo_q);
Xd_baixa = Xd/beta^2;
Xq_baixa = Xq/beta^2;

Z2_baixa = Z2/beta^2;

Zmd_menos = 1/(1/(1i*Xmd)+2/(Rrd +1i*2*Xlrd));
Zmq_menos = 1/(1/(1i*Xmq)+2/(Rrq +1i*2*Xlrq));
Zdf_menos = Zmd_menos -R/3 + 1i*Xls;
Zdb_menos = Zmd_menos - R + 1i*Xls;
Zqf_menos = Zmq_menos -R/3 + 1i*Xls;
Zqb_menos = Zmq_menos -R + 1i*Xls;
Zeq_tras1 = (Zdb_menos - Zqb_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);
Zeq_tras2 =(Zdf_menos + Zqf_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);
Re_neg = 0.5*real(-Zmd_menos -Zmq_menos + ((Zmd_menos - Zmq_menos)^2)/(Zqf_menos + Zdf_menos));

Pwatt_ferro = 3.0345;
P_ferro = Pwatt_ferro;

tol = 1e-4;

for int = 1:length(delta)
    
    flag = 1;
    
    % Para componentes naturais.
    while (flag > 0)
        
        V1n = (V_baixa)*cos(delta(int)*pi/180)*1i - (V_baixa)*sin(delta(int)*pi/180);
        A_1n = [(V_baixa)*cos(delta(int)*pi/180)-E_baixa ; (V_baixa)*sin(delta(int)*pi/180)];
        B_1n = [1i*Xd_baixa R; -R -1i*Xq_baixa];
        In = B_1n\A_1n;
        I1n = 1i*(In(1)) + (In(2));
        Z1n = V1n/I1n;
        
        a1n = 1 + ZCap/Z1n;
        a2n = 1 + ZCap/Z2_baixa;
        
        V2n = V1n*(beta-1i*a1n)/(beta+1i*a2n);
        
        Vx = V1n + V2n;
        divisor = abs(Vx)/(Vm_rms/beta);
        V_baixa = V_baixa/divisor;
        
        P_ferro(:,int) = Pwatt_ferro;
        
        if (abs(divisor) < 1 + tol)  && (abs(divisor) > 1 - tol)
            flag = 0;
        end
        
    end
    
    I2n = V2n/Z2_baixa;
    
    %     I_sim_nat_baixa(:,int) = [I1n;I2n];
    %     V_sim_nat_baixa(:,int) = [V1n;V2n];
    
    % Para componentes normais eixos "para frente" e "para tras".
    
    V_alpha_beta = (1/sqrt(2)).*a1*[V1n; V2n];
    I_alpha_beta = (1/sqrt(2)).*a1*[I1n; I2n];
    V_1_2 = (1/sqrt(2)).*a2*V_alpha_beta;
    I_1_2 = (1/sqrt(2)).*a2*I_alpha_beta;
    V_m_a = a3*[V1n; V2n];
    I_m_a = a4*[I1n; I2n];
    PF(int) = cos(phase(V_m_a(1)) - phase(sum(I_m_a)));
    P_in2(int) = abs(V_m_a(1))*abs(sum(I_m_a))*PF(int);
    S_in(int) = abs(V_m_a(1))*abs(sum(I_m_a));
    
    %     V_sim_baixa(:,int) = (1/sqrt(2))*[1 1i; 1 -1i]*V_alpha_beta;
    %     I_sim_baixa(:,int) = (1/sqrt(2))*[1 1i; 1 -1i]*I_alpha_beta;
    %
    %     %Tensao nos enrolamentos auxiliar e principal.
    %     V_ma_aux(:,int) = [beta beta;1i -1i]*V_sim_nat_baixa(:,int);
    %     I_ma_aux(:,int) = [1/beta 1/beta;1i -1i]*I_sim_nat_baixa(:,int);
    
    %%%%Torques & Pot?ncias
    
    %Potencia de entrada.
    P_in(int) = 2*real(conj(I1n)*V1n) + 2*real(conj(I2n)*V2n);
    
    %Potencia pra frente.
    T_rel(int) = P*(Xd_baixa - Xq_baixa)*real(In(1))*imag(In(2))/w; % Olhar as correntes
    T_alin(int) = P*E_baixa*imag(In(2))/w; % Olhar corrente
    T_frente(int) = T_alin(int) + T_rel(int);
    P_frente = T_frente(int)*w;
    
    %Potencia nos enrolamentos.
    P_cu_estator_frente(int) = R*power(abs(I1n),2);
    P_cu_estator_tras(int) = R*power(abs(I2n),2);
    T_cu_estator_tras(int) = P_cu_estator_tras(int)/w;
    
    % Potencia para tras
    P_tras(int) = real(2*conj(I2n)*V2n) - P_cu_estator_tras(int);
    T_tras(int) = P_tras(int)/w;
    
    %Potencia pra tras.
    %     If_menos(int) = -V2n*Zeq_tras1;
    %     Ib_menos(int) = V2n*Zeq_tras2;
    %     T_f_menos(int) = 2*(1/(w))*(R/(-3))*abs(If_menos(int))^2;
    %     T_b_menos(int) = 2*(1/(w))*Re_neg*abs(Ib_menos(int))^2;
    %     T_tras_tot(int) = T_f_menos(int) + T_b_menos(int);
    %     P_tras_tot(int) = w*T_tras_tot(int);
    %     T_tras(int) = T_tras_tot(int) - T_cu_estator_tras(int);
    %     P_tras(int) = w*T_tras(int);
    
    
    %Potencia no condutor do rotor.
    V_entreferro_tras(int) = V2n - ((R1 + 1i*Xls)/(beta^2))*(I2n);
    V_entreferro_q(int) = V_entreferro_tras(int)*(negativo_eixo_q/(negativo_eixo_d + negativo_eixo_q));
    V_entreferro_d(int) = V_entreferro_tras(int)*(negativo_eixo_d/(negativo_eixo_d + negativo_eixo_q));
    I_d_negativo(int) = (V_entreferro_d(int))/((0.5*(Rrd/2 + 1i*Xlrd ))/beta^2);
    I_q_negativo(int) = (V_entreferro_q(int))/((0.5*(Rrq/2 + 1i*Xlrq))/beta^2);
    Wcur1(int) = Rrd/(beta^2)*(abs(I_d_negativo(int)))^2 ;
    Wcur2(int) = Rrq/(beta^2)*(abs(I_q_negativo(int)))^2;
    Wcur(int) = Wcur1(int) +  Wcur2(int);
    T_Wcur(int) = Wcur(int)/w;
    
    P_entreferro(int) = P_in(:,int) - P_cu_estator_frente(:,int) - P_cu_estator_tras(:,int) - P_ferro(1);
    % T_entreferro(int) = P_entreferro(:,int)/w;
    
    %P_out(int) = w*(T_frente(int) + T_tras(int));
    %P_entreferro(int) = P_in(int) - P_out(int);
    T_entreferro(int) = P_entreferro(:,int)/w;
    
    %     P_tras_frente(int) = -P_entreferro(:,int) + P_frente;
    %     T_tras_frente(int) = P_tras_frente(:,int)/(w);
    
    Eff(int) = 100*P_entreferro(int)/P_in(int); %100 multiplicado para dar em %
    
    if Eff(int) > 100
        Eff(int) = 100;
    elseif Eff(int) < 0
        Eff(int) = 0;
    end
    
    T_total(int) = 10.2*(T_tras(int) + T_frente(int)); % 10.2 multiplicado para dar kgf.cm
    
end

plot(T_total, Eff, 'b:.');
grid on
hold on
LSPM_Embraco
plot(Torque_med, Eff_med, 'r');
plot(Speed(:,1), Speed(:,2),'g');
xlabel('Torque (kgf.cm)');
ylabel('Eficiencia (%)');

% %grafico_potent = plot(delta,P_entreferro(:,k));
% subplot(3,3,1);
% plot(delta, P_in, 'b:.');
% %grafico_potent = plot(delta,Eff);
% grid on
% title('Potencia de entrada (Pot. elec no SPEED)');
% xlabel('Delta (graus)');
% ylabel('Potencia (W)');
% 
% subplot(3,3,2);
% plot(delta, T_total, 'b:.');
% grid on
% title('Torque total');
% xlabel('Delta (graus)');
% ylabel('Potencia(W)')
% 
% subplot(3,3,3);
% plot(delta, Eff, 'b:.');
% grid on
% title('Efficiencia');
% xlabel('Delta (graus)');
% ylabel('Efficiencia')
% 
% subplot(3,3,4);
% plot(delta,P_cu_estator_frente, 'b:.');
% grid on
% title('Potencia Cu estator para frente');
% xlabel('Delta (graus)');
% ylabel('Potencia (W)');
% 
% subplot(3,3,5);
% plot(delta, P_cu_estator_tras, 'b:.');
% grid on
% title('Potencia Cu estator para tras');
% xlabel('Delta (graus)');
% ylabel('Potencia (W)');
% 
% subplot(3,3,6);
% plot(delta,Wcur, 'b:.');
% grid on
% title('Potencia no condutor do rotor');
% xlabel('Delta (graus)');
% ylabel('Potencia (W)');
% 
% subplot(3,3,7);
% plot(delta,T_alin, 'b:.');
% grid on
% title('Torque Alinhamento');
% xlabel('Delta (graus)');
% ylabel('Torque');
% 
% subplot(3,3,8);
% plot(delta,T_rel, 'b:.');
% grid on
% title('Torque Relutancia');
% xlabel('Delta (graus)');
% ylabel('Torque');
% 
% subplot(3,3,9);
% plot(T_total, Eff, 'b:.');
% grid on
% title('Efici?ncia (%) x Torque');
% xlabel('Torque (kgf.cm)');
% ylabel('Efici?ncia(%)');

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

