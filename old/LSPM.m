clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                          LSPM                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nome: Guilherme Corrêa Danielski                                                               %%
%% Matrícula: 12100578                                                                            %%
%% Engenharia Elétrica                                                                            %%
%% Versão: 1.0v                                                                                   %%
%% Tempo de Implementação: 5h                                                                     %%
%% Engenharia Elétrica - Sétima Fase                                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input
f = 60;             % frequency (Hz)
Eo = 105.6;         % No-load induced voltage (V)
Vm = 230;           % Complex main supply voltage (V)
C_run = 100e-6;     % Run capacitor value (F)
Xc = 2*pi*f*C_run;  % Capacitive impedance (Ohms)
Nm = 1;             % Number of turns on main stator winding
Na = 1;             % Number of turns on auxiliary stator winding
Ra = 27.7;          % Auxiliary stator winding resistance (Ohms)
Xla = 21.1;         % Auxiliary stator leakage reactance (Ohms)
Xd = 90.3;          % Synchronous reactance for d axis (Ohms)
Xq = 160.5;         % Synchronous reactance for q axis (Ohms)
Xmq = 56.2;         % Magnetization reactance for q axis (Ohms)
Xmd = 49.6;         % Magnetization reactance for d axis (Ohms)
Rrd = 11;           % Rotor resistance for d axis (Ohms)
Rrq = 25;           % Rotor resistance for q axis (Ohms)
Xlrd = 49.6;        % Rotor leakage reactance for d axis (Ohms)
Xlrq = 56.2;        % Rotor leakage reactance for q axis (Ohms)
m = 1;              % Phases number
P = 8;              % Poles number
csi = pi/2;         % Shift electrical angle between stator windings (rad)

%% Parameters
a = Na/Nm;                  % Turns ratio (auxiliary/main)
beta = 1/a;                 % Turns ration (main/auxiliary)
h = (1+cos(csi))/sin(csi);  % Eq. (11)
Xls = power(beta,2)*Ra;     % Citation [1]
Xlm = Xls;                  % Citation [1]
Rs = power(beta,2)*Xla;     % Equivalent stator winding resistance (Ohms)
Rm = Rs;                    % Main stator winding resistance (Ohms)
Ns = 120*f/P;               % Motor syncronous speed (RPM)
w = 4*pi*f/P;       % Synchronous speed (rad/s)

%% Initializing the loop

int = 0;

% Initializing the variables
pre_allocate = length(0:0.001:1);

T_pulse_2sf = zeros(pre_allocate,1); T_pulse_4_2sf = zeros(pre_allocate,1); T_pulse_2f = zeros(pre_allocate,1);
T_pulse_2_2sf = zeros(pre_allocate,1); T_pulse_sf = zeros(pre_allocate,1); T_pulse_2_sf = zeros(pre_allocate,1);
Tm = zeros(pre_allocate,1); RPM = zeros(pre_allocate,1); T_f_pos = zeros(pre_allocate,1);
T_b_pos = zeros(pre_allocate,1); T_f_neg = zeros(pre_allocate,1); T_b_neg = zeros(pre_allocate,1); 
Tavg = zeros(pre_allocate,1);
 
%% Loop
for s = 0:0.001:1           % Slip

    int = int + 1;
    
    Z_pos = Rs + 1i*Xls + 0.5*(1i*Xmd*(Rrd/s + 1i*Xlrd)/(Rrd/s + 1i*(Xmd + Xlrd)) ...
        + 1i*Xmq*(Rrd/s + 1i*Xlrd)/(Rrd/s + 1i*(Xmq + Xlrq)));                                  % Eq. (12)
    Z_neg = Rs + 1i*Xls + 0.5*(1i*Xmd*(Rrd/(2 - s) + 1i*Xlrd)/(Rrd/(2 - s) + ... 
        1i*(Xmd + Xlrd)) + 1i*Xmq*(Rrd/(2 - s) + 1i*Xlrd)/(Rrd/(2 - s) + 1i*(Xmq + Xlrq)));     % Eq. (13)

    % Eq. (11)
    a1 = 1 - (1i*Xc/Z_pos)*(1 + 1i/(beta*tan(csi)));                
    a2 = 1 - (1i*Xc/Z_neg)*(1 - 1i/(beta*tan(csi)));                



    V_pos = Vm*(sqrt(2)*sin(csi)/beta)*(h*beta + 1i*a2)/(a1+a2);    % Eq. (9)
    V_neg = Vm*(sqrt(2)*sin(csi)/beta)*(h*beta - 1i*a1)/(a1+a2);    % Eq. (10)
    
    % Eq. (2)
    Vd_pos = V_pos;                                                 
    Vq_pos = 1i*V_pos;                                              
    
    % Eq. (3)
    Vd_neg = V_neg;                                                 
    Vq_neg = 1i*V_neg;                                              

    % Eq. (30)
    Zmd_pos = inv(1/(1i*Xmd) + s/(Rrd + 1i*s*Xlrd));
    Zmq_pos = inv(1/(1i*Xmq) + s/(Rrq + 1i*s*Xlrq));
    Zmd_neg = inv(1/(1i*Xmd) + (2 - s)/(Rrd + 1i*(2 - s)*Xlrd));
    Zmq_neg = inv(1/(1i*Xmq) + (2 - s)/(Rrq + 1i*(2 - s)*Xlrq));

    % Eq. (22)
    Xd_pos = -1i*(Zmd_pos + 1i*Xls);
    Xd_neg = -1i*(Zmd_neg + 1i*Xls);
    Xq_pos = -1i*(Zmq_pos + 1i*Xls);
    Xq_neg = -1i*(Zmq_neg + 1i*Xls);

    Zd_pos = 1i*Xd_pos;
    Zd_neg = 1i*Xd_neg;
    Zq_pos = 1i*Xq_pos;
    Zq_neg = 1i*Xq_neg;

    % Eq. (23)
    Zd_f_pos = Zmd_pos + Rs + 1i*Xls;
    Zq_f_pos = Zmq_pos + Rs + 1i*Xls;

    Zd_b_pos = Zmd_pos + (Rs/(2*s - 1)) + 1i*Xls;
    Zq_b_pos = Zmq_pos + (Rs/(2*s - 1)) + 1i*Xls;

    % Eq. (24)
    Zd_f_neg = Zmd_neg + (Rs/(2*s - 3)) + 1i*Xls;
    Zq_f_neg = Zmq_neg + (Rs/(2*s - 3)) + 1i*Xls;

    Zd_b_neg = Zmd_neg - Rs + 1i*Xls;
    Zq_b_neg = Zmq_neg - Rs + 1i*Xls;

    %% Calculo do Torque medio
    I_f_pos = V_pos*(Zd_b_pos + Zq_b_pos)/(Zd_f_pos*Zq_f_pos + Zd_b_pos*Zq_f_pos);              % Eq. (18)
    I_b_pos = -V_pos*(Zd_f_pos - Zq_f_pos)/(Zd_f_pos*Zq_f_pos + Zd_b_pos*Zq_f_pos);             % Eq. (19)
    I_f_neg = -V_neg*(Zd_b_neg - Zq_b_neg)/(Zd_f_neg*Zq_f_neg + Zd_b_neg*Zq_f_neg);             % Eq. (20)
    I_b_neg = V_neg*(Zd_f_neg - Zq_f_neg)/(Zd_f_neg*Zq_f_neg + Zd_b_neg*Zq_f_neg);              % Eq. (21)

    % Eq. (29)
    Re_pos = 0.5*real(Zmd_pos + Zmq_pos - power(Zmd_pos - Zmq_pos,2)/(Zd_b_pos + Zq_b_pos));    
    Re_neg = 0.5*real(- Zmd_neg - Zmq_neg + power(Zmd_neg - Zmq_neg,2)/(Zd_f_neg + Zq_f_neg));

    T_f_pos(int) = m*P*0.5*Re_pos*power(I_f_pos,2)/w;                       % Eq. (25)
    T_b_pos(int) = m*P*0.5*Rs*power(I_b_pos,2)/(w*(2*s - 1));               % Eq. (26)
    T_f_neg(int) = m*P*0.5*Rs*power(I_f_neg,2)/(w*(2*s - 3));               % Eq. (27)
    T_b_neg(int) = m*P*0.5*Re_neg*power(I_b_neg,2)/w;                       % Eq. (28)

    Tavg(int) = T_f_pos(int) + T_b_pos(int) + T_f_neg(int) + T_b_neg(int);  % Eq. (31)

    %% Pulsating Torque

    D_pos = power(Rs,2) + (1 - 2*s)*Xd_pos*Xq_pos + 1i*s*Rs*(Xd_pos + Xq_pos);          % Eq. (7)
    D_neg = power(Rs,2) + (2*s - 3)*Xd_neg*Xq_neg + 1i*(2 - s)*Rs*(Xd_neg + Xq_neg);    % Eq. (8)

    % Eq. (5)
    Id_pos = 1i*V_pos*(-Rs + 1i*(2*s - 1)*Xd_pos)/D_pos;                                
    Iq_pos = -V_pos*(-Rs + 1i*(2*s - 1)*Xq_pos)/D_pos;  
    
    % Eq. (6)
    Id_neg = V_neg*(Rs + 1i*(3*s - 2)*Xd_neg)/D_neg;
    Iq_neg = 1i*V_neg*(Rs + 1i*(3*s - 2)*Xq_pos)/D_neg;
    
    % Eq. (32)
    Idm = -power(1 - s,2)*(Xq - Xc)*Eo/(power(Rs,2) + Xd*(Xq - Xc)*power(1 - s,2));
    Iqm = -(1 - s)*Rs*Eo/(power(Rs,2) + Xd*(Xq - Xc)*power(1 - s,2));

    % Eq. (33)
    psi_dm = (Xd*Idm + Eo)/w;
    psi_qm = (Xq - Xc)*Iqm/w;
    
    % Eq. (4)
    psi_d_pos = Xd_pos*1i*s*Id_pos/w;
    psi_q_pos = Xq_pos*1i*s*Iq_pos/w;
    psi_d_neg = Xd_neg*1i*s*Id_neg/w;
    psi_q_neg = Xq_neg*1i*s*Iq_neg/w;
    
    % Reluctance pulsating torque
    T_pulse_2sf(int) = 0.25*m*P*abs(psi_q_pos*Id_pos - psi_d_pos*Iq_pos);                           % Eq. (35)
    T_pulse_4_2sf(int) = 0.25*m*P*abs(psi_q_neg*Id_neg - psi_d_neg*Iq_neg);                         % Eq. (36)
    
    % Unbalanced stator pulsating torque
    T_pulse_2f(int) = 0.25*m*P*abs(psi_q_pos*Id_neg - psi_d_pos*Iq_neg);                            % Eq. (37)
    T_pulse_2_2sf(int) = 0.25*m*P*abs(psi_q_neg*Id_pos - psi_d_neg*Iq_pos);                         % Eq. (38)
    
    % PM excitantion pulsating torque
    T_pulse_sf(int) = 0.5*m*P*abs(psi_q_pos*Idm + psi_qm*Id_pos - psi_d_pos*Iqm - psi_dm*Iq_pos);   % Eq. (39)
    T_pulse_2_sf(int) = 0.5*m*P*abs(psi_q_neg*Idm + psi_qm*Id_neg - psi_d_neg*Iqm - psi_dm*Iq_neg); % Eq. (40)

    %% Magnetic Breaking Torque
    Tm(int) = 0.5*P*sin(csi)*(beta*psi_dm*Iqm - psi_qm*Idm/beta);   % Eq. (34)

    %% Speed
    RPM(int) = (1-s)*Ns; % Motor assyncronous speed.
end

%% Plotting

% COMENTÁRIO: Os torques tem componentes complexas. Isso é normal?

% Average torque and magnetic breaking torque x Motor speed
figure(1);
plot(RPM,Tavg,RPM,Tm);
xlabel('Motor Speed (RPM)');
ylabel('Motor Torque [N.m]');
legend('Tavg','Tm');
title('Average torque and magnetic breaking torque x Motor speed');

% Average cage torque components x Motor speed
figure(2);
plot(RPM,T_f_pos,RPM,T_f_neg,RPM,T_b_pos,RPM,T_b_neg);
xlabel('Motor Speed (RPM)');
ylabel('Motor Torque [N.m]');
legend('T_f_pos','T_f_neg','T_b_pos','T_b_neg');
title('Average cage torque components x Motor speed');

% Pulsating torques x Motor speed
figure(3);
plot(RPM,T_pulse_2sf,RPM,T_pulse_4_2sf,RPM,T_pulse_2f,RPM,T_pulse_2_2sf,RPM,T_pulse_sf,RPM,T_pulse_2_sf);
xlabel('Motor Speed (RPM)');
ylabel('Motor Torque [N.m]');
legend('T_pulse_2sf','T_pulse_4_2sf','T_pulse_2f','T_pulse_2_2sf','T_pulse_sf','T_pulse_2_sf');
title('Pulsating torques x Motor speed');

% Sum of all torques x Motor speed
figure(4);
plot(RPM,Tavg+Tm+T_pulse_2sf+T_pulse_4_2sf+T_pulse_2f+T_pulse_2_2sf+T_pulse_sf+T_pulse_2_sf);
xlabel('Motor Speed (RPM)');
ylabel('Motor Torque [N.m]');
title('Sum of all torques x Motor speed');