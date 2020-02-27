%% Parte  Sincrona
global w raiz f Cap Vlin P
global delta_carga P_in T_frente Eff P_cu_tras P_cu_frente P_out T_rel Wcur T_total T_tras;

% Entradas Temporarias
E = 160.3917;   % ???

Ld = 0.8869;    % Falta                     
Lq = 2.0268;    % Falta
% Xls = 19.7167;  % Falta Como Calcular???
Xmd = Xmd/1i;
Xmq = Xmq/1i;
X2d = X2d/1i;
X2q = X2q/1i;
X1d = X1d/1i;
Xls = X1d;

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
    
    P_in(int) = 2*(real(V1n*conj(I1n)) + real(V2n*conj(I2n)));
    
    T_rel(int) = 2*P*(abs(Xd_baixa) - abs(Xq_baixa))*real(Im(1))*imag(Im(2))/w;
    T_al(int) = 2*P*imag(Im(2))*E_baixa/w;
    T_frente(int) = T_rel(int) + T_al(int);
    P_cup = Rp*real(I_m_a(1)*conj(I_m_a(1)));
    P_cua = Ra*real(I_m_a(2)*conj(I_m_a(2)));
    P_cu_frente(int) = 2*R*power(abs(I1n),2);
    P_cu_tras(int) = 2*R*power(abs(I2n),2);
    P_tras = 2*conj(I2n)*V2n - P_cu_tras(int);
    T_tras(int) = P_tras/w;
    Wcur(int) = 2*real(V2n*conj(I2n)) - P_cu_tras(int);
    P_air(int) = P_in(int) - P_cup - P_cua;
    P_out(int) = P_air(int) - Wcur(int) - P_ferro;
    
    T_total(int) = 10.2*(T_tras(int) + T_frente(int));
    
    Eff(int) = (P_out(int)/P_in(int))*100;
    
    if Eff(int) > 100
        Eff(int) = 1;
    elseif Eff(int) < 0
        Eff(int) = 0;
    end

end

cd(raiz)
cd interface
LSPM_Plot_Manager_Sincrono;