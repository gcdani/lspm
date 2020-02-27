%% Parte  Sincrona

%% Observacoes
% Capacitor de partida e de run

%% Valores Globais
 global w raiz f Cap Vlin P
 global delta_carga P_in T_frente Eff P_cu_tras P_cu_frente P_out T_rel Wcur T_total T_tras T_al;

%% Valores de entrada
a_XL = 0.91;    % Fatores de ajustes do SPEED

% Entradas Temporarias
E = 154.0069;   % Tensao induzida
f = 50;         % Frequencia da rede
P = 1;          % Par de polos
w = 2*pi*f;
Ld = a_XL*807.1198e-3;    % Indutancia sincrona de eixo direto                     
Lq = a_XL*1331.7623e-3;   % Indutancia sincrona de eixo em quadratura
Rp = 27.7141;   % Resistencia do enrolamento principal
Ra = 13.5156;   % Resistencia do enrolamento auxiliar
R2q = 49.2121;  % Resistencia sincrona de eixo direto
R2d = 49.2121;  % Resistencia sincrona de eixo em quadratura
R2 = R2d;       % Resistencia total da gaiola
X1d = a_XL*21.0763; % Reatancia 
Xls = X1d;
X2q = a_XL*12.6949;
X2d = X2q;
Lmd = a_XL*757.0975e-3;
Lmq = a_XL*1281.7401e-3;
Cap = 4.75e-6;  % Capacitor de partida
Vlin = 220;     % Tensao de entrada


Xmd = w*Lmd;
Xmq = w*Lmq;


% Xmd = Xmd/1i;
% Xmq = Xmq/1i;
% X2d = X2d/1i;
% X2q = X2q/1i;
% X1d = X1d/1i;
% Xls = X1d;

Xd = w*Ld;  
Xq = w*Lq;

%% Inicializacao de variaveis
P_ferro = 3.0345;           % Perda no ferro
delta_carga = 0:(90/25):90; % Angulo de carga

%% Inicio Programa

% beta = Neep/Neea;         % Relacao enrolamento por auxiliar
beta = 1.7836889285222159;  % Temporario
Xcap = 2*pi*f*Cap;          % Reatancia devido ao capacitor

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
tol = 1e-5;     % Tolerancia


Z2 = Rp + 1i*X1d + 1i*.5*Xmd*((.5*R2 + 1i*X2d)/(1i*Xmd + .5*R2 + 1i*X2d)) + ...
            1i*.5*Xmq*((.5*R2 + 1i*X2q)/(1i*Xmq + .5*R2 + 1i*X2q));

eixo_d = 1i*.5*Xmd*((.5*R2 + 1i*X2d)/(1i*Xmd + .5*R2 + 1i*X2d));
eixo_q = 1i*.5*Xmq*((.5*R2 + 1i*X2q)/(1i*Xmq + .5*R2 + 1i*X2q));
         
Z2_baixa = Z2/power(beta,2);
ZCap = 1/(1i*Xcap);

%% Calculo das sequencias positivas e negativas
for int = 1:length(delta_carga)

    flag = 1;
    teste = 0;
    
    while(teste==0)
        
        V1n = V_baixa*cos(delta_carga(int)*pi/180)*1i - V_baixa*sin(delta_carga(int)*pi/180);   % Tensao de sequencia positiva
        
        A1 = [V_baixa*cos(delta_carga(int)*pi/180) - E_baixa;
            V_baixa.*sin(delta_carga(int)*pi/180)];
        B1 = [1i*Xd_baixa, -R;
            -R, -1i*Xq_baixa];
        
        Im = B1\A1;         % Id1 e Iq1
        I1n = sum(Im);      % Corrente de sequencia positiva
        Z1n = V1n/I1n;      % Impedancia de sequencia positiva

        a1n = 1 + ZCap/Z1n;
        a2n = 1 + ZCap/Z2_baixa;

        V2n = V1n*(beta-1i*a1n)/(beta+1i*a2n);  % Tensao de sequencia negativa
        Vs = V1n + V2n;     % Tensao de entrada
        
        flag = abs(Vs)*beta/Vlin;
        V_baixa = V_baixa/flag;
        
        if (flag < 1 + tol && flag > 1 - tol)
            teste = 1;
        end
    
    end
    
    I2n = V2n/Z2_baixa;     % Corrente de sequencia negativa
    
    V_alpha_beta = (1/sqrt(2)).*a1*[V1n; V2n];
    I_alpha_beta = (1/sqrt(2)).*a1*[I1n; I2n];
    V_1_2 = (1/sqrt(2)).*a2*V_alpha_beta;
    I_1_2 = (1/sqrt(2)).*a2*I_alpha_beta;        
    V_m_a = a3*[V1n; V2n];
    I_m_a = a4*[I1n; I2n];
    
    %% Calculo das Saidas
    
    % Potencia de entrada
    P_in(int) = 2*(real(V1n*conj(I1n)) + real(V2n*conj(I2n)));
    
    % Torque para frente
    T_rel(int) = P*(abs(Xd_baixa) - abs(Xq_baixa))*real(Im(1))*imag(Im(2))/w;
    T_al(int) = P*imag(Im(2))*E_baixa/w;
    T_frente(int) = T_rel(int) + T_al(int);
    
    % Perda no condutor
    P_cu_frente(int) = R*abs(power(I1n,2));
    P_cu_tras(int) = R*abs(power(I2n,2));
    
    % Torque para tras
    P_tras = real(conj(I2n)*V2n) - P_cu_tras(int)/w;
    T_tras(int) = P_tras/w;
     
    %Potencia no condutor do rotor.   
    Wcur(int) = P*2*T_tras(int);
    
    % Potencia de saida
    P_out(int) = P_in(int) - P_cu_frente(int) - P_cu_tras(int) - P_ferro;
    T_out = P_out(int)/w;
    
    % Potencia do motor
    T_total(int) = 10.2*(T_tras(int) + T_frente(int)); % 10.2 multiplicado para dar kgf.cm
    
    % Eficiencia
    Eff(int) = (P_out(int)/P_in(int))*100;
    
    if Eff(int) > 100
        Eff(int) = 100;
    elseif Eff(int) < 0
        Eff(int) = 0;
    end

end

cd(raiz)
cd interface
LSPM_Plot_Manager_Sincrono;