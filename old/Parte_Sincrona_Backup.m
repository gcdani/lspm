%% Parte  Sincrona
 global w raiz f Cap Vlin P
 global delta_carga P_in T_frente Eff P_cu_tras P_cu_frente P_out T_rel Wcur T_total T_tras T_al;

%%%% Fatores de ajustes do SPEED
a_XL = 0.91;

% Entradas Temporarias
E = 154.0069;   % ???
f = 50;
P = 1;
w = 2*pi*f;
Ld = a_XL*807.1198e-3;    % Falta                     
Lq = a_XL*1331.7623e-3;    % Falta
Rp = 27.7141;
Ra = 13.5156;
R2q = 49.2121; 
R2d = 49.2121; 
R2 = R2d;
X1d = a_XL*21.0763;
Xls = X1d;
X2q = a_XL*12.6949;
X2d = X2q;
Lmd = a_XL*757.0975e-3;
Lmq = a_XL*1281.7401e-3;
Cap = 4.75e-6;
Vlin = 220;


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


P_ferro = 3.0345;
delta_carga = 0:(90/25):90;

% Fim

% beta = Neep/Neea;
beta = 1.7836889285222159;
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
        
%         a1n = 1i + 1i*ZCap/Z1n;
%         a2n = -1i - 1i*ZCap/Z2_baixa;

        a1n = 1 + ZCap/Z1n;
        a2n = 1 + ZCap/Z2_baixa;

        %V2n = V1n*((beta - a1n)/(a2n - beta));
        V2n = V1n*(beta-1i*a1n)/(beta+1i*a2n);
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
    
    Vd = abs(V_alpha_beta(1)) + abs(V_alpha_beta(2))*exp(1i*phase(sum(V_alpha_beta)));
    Vq = -1i*Vd;
    
    
    H = Rp + 1i*2*Xd;
    G = Xq;
    
    V_m_a = a3*[V1n; V2n];
    I_m_a = a4*[I1n; I2n];
    
    P_in(int) = 2*(real(V1n*conj(I1n)) + real(V2n*conj(I2n)));
    
    T_rel(int) = P*(abs(Xd_baixa) - abs(Xq_baixa))*real(Im(1))*imag(Im(2))/w;
    T_al(int) = P*imag(Im(2))*E_baixa/w;
    T_frente(int) = T_rel(int) + T_al(int);
    % P_cup = R1d*I_m_a(1);
    % P_cua = R2d*I_m_a(2);
    P_cu_frente(int) = 2*R*abs(power(I1n,2));
    P_cu_tras(int) = R*abs(power(I2n,2));
    % P_cu_rotor = 2*V2n*conj(I2n) - P_cu_tras;
    % P_air = P_in - P_cup - P_cua;
    % P_out = P_air - P_cu_rotor - P_ferro;
    
    %Potencia pra tras.
    If_menos = -V2n*Zeq_tras1;
    Ib_menos = V2n*Zeq_tras2;
    
    T_f_menos = 2*(R/(-3))*power(abs(If_menos),2)/w;
    T_b_menos = 2*Re_neg*power(abs(Ib_menos),2)/w; 
    T_tras_tot = T_f_menos + T_b_menos;
    P_tras_tot = w*T_tras_tot;
    
    % T_tras = T_tras_tot - P_cu_tras(int)/w;
    % P_tras = w*T_tras;
    P_tras = real(conj(I2n)*V2n) - P_cu_tras(int)/w;
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
    T_out = P_out(int)/w;
    P_tras_frente(int) = - P_out(int) + T_frente(int)*w;
    T_tras_frente = P_tras_frente(int)/w;
    
    T_total(int) = (T_tras(int) + T_frente(int)); % 10.2 multiplicado para dar kgf.cm
    
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