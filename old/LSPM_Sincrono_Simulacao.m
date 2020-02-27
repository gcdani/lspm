%%% Programa para calculo dos valores de eficiencia, 
%%% Angulo de carga do motor estudado.
%%% Dados de entrada extraidos de PC - 45 do WORD

global freqRede Cap Vm_rms s beta E Ld Lq p R1 R2 Rrd Rrq Xd Xq Xlrd Xlrq Xls Xmd Xmq;
global delta P_entrada T_frente Eff P_cu_estator_tras P_cu_estator_frente P_entreferro T_rel Wcur;

delta = 0:(90/2999):90;


V_baixa = Vm_rms/beta;
              
E_baixa = E/beta;

R = R1/beta^2;

ZCap = 1/(1i*2*pi*Cap*freqRede);   

ws = 2*pi*freqRede;    

negativo_eixo_d = 0.5*1i*Xmd*(Rrd/2 + 1i*Xlrd)/(1i*Xmd + Rrd/2 + 1i*Xlrd); 

negativo_eixo_q = 0.5*1i*Xmq*(Rrq/2 + 1i*Xlrq)/(1i*Xmq + Rrq/2 + 1i*Xlrq);

Z2_alta = R1 + 1i*Xls + 0.5*(negativo_eixo_d + negativo_eixo_q);

%Xd_baixa = 78.387; %Xd/beta^2;

%Xq_baixa = 320;%Xq/beta^2;

Z2_baixa = 12.3 + 11.17i;

Z2_baixa = 15.326 + 10.73i; %Z2_alta/beta^2;


Zmd_menos = 1/(1/(1i*Xmd)+2/(Rrd +1i*2*Xlrd));

Zmq_menos = 1/(1/(1i*Xmq)+2/(Rrq +1i*2*Xlrq));

Zdf_menos = Zmd_menos -R/3 + 1i*Xls;

Zdb_menos = Zmd_menos - R + 1i*Xls;

Zqf_menos = Zmq_menos -R/3 + 1i*Xls;

Zqb_menos = Zmq_menos -R + 1i*Xls;

Zeq_tras1 = (Zdb_menos - Zqb_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);

Zeq_tras2 =(Zdf_menos + Zqf_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);

Re_neg = 0.5*real(-Zmd_menos -Zmq_menos + ((Zmd_menos - Zmq_menos)^2)/(Zqf_menos + Zdf_menos));

Xd_baixa = 83;

Xq_baixa = 530;

tol = 1e-6;

a2n = -1i -1i*ZCap/Z2_baixa;

I_1n = zeros(1,length(delta)); %Kleyton sugestão

for k = 1:length(delta)
    
    sentinela = 1;
    
    % Para componentes naturais.
    while (sentinela > 0)
        
    V_1n(:,k) = (V_baixa)*cos(delta(k)*pi/180)*1i -(V_baixa)*sin(delta(k)*pi/180); 
    A_1n = [(V_baixa)*cos(delta(k)*pi/180)-E_baixa ; (V_baixa)*sin(delta(k)*pi/180)];
    B_1n = [-1i*Xd_baixa -1i*R; -R Xq_baixa];
    In = inv(B_1n)*A_1n;
    I_1n(:,k) = In(1) + In(2);
    Z_1n(:,k) = V_1n(:,k)/I_1n(:,k);
    
    a1n = 1i + 1i*ZCap/Z_1n(:,k);
    
    V_2n(:,k) = V_1n(:,k)*(beta-a1n)/(a2n-beta);
    
    Vx = V_1n(:,k) + V_2n(:,k);
    divisor = abs(Vx)/(Vm_rms/beta);
    V_baixa = V_baixa/divisor;
    
    P_ferro(:,k) = 3.0345;
    
    if (abs(divisor) < 1 + tol)  && (abs(divisor) > 1 - tol)
        sentinela = 0;    
    end
    
end
       
    I_2n(:,k) = V_2n(:,k)/Z2_baixa;     
    
    I_sim_nat_baixa(:,k) = [I_1n(:,k);I_2n(:,k)];
    V_sim_nat_baixa(:,k) = [V_1n(:,k);V_2n(:,k)];
    
    % Para componentes normais eixos "para frente" e "para tras".
    
    V_alphabeta_baixa(:,k) = [1i -1i;1 1]*[V_1n(:,k);V_2n(:,k)]; 
    I_alphabeta_baixa(:,k) = [1i -1i;1 1]*[I_1n(:,k);I_2n(:,k)]; 
    
    V_sim_baixa(:,k) = (1/sqrt(2))*[1 1i; 1 -1i]*V_alphabeta_baixa(:,k);
    I_sim_baixa(:,k) = (1/sqrt(2))*[1 1i; 1 -1i]*I_alphabeta_baixa(:,k);
    
    %Tensao nos enrolamentos auxiliar e principal.
    V_ma_aux(:,k) = [beta beta;1i -1i]*V_sim_nat_baixa(:,k);
    I_ma_aux(:,k) = [1/beta 1/beta;1i -1i]*I_sim_nat_baixa(:,k);
   
    %%%%Torques & Potências
    
    %Potencia de entrada.
    P_entrada(:,k) = 2*real(conj(I_1n(:,k))*V_1n(:,k)) + 2*real(conj(I_2n(:,k))*V_2n(:,k));     
    
    %Potencia pra frente.
    T_rel(:,k) = 2*(p/ws)*(Xd_baixa-Xq_baixa)*real(I_sim_nat_baixa(1,k))*imag(I_sim_nat_baixa(1,k));
    T_alin(:,k) = 2*(p/ws)*E_baixa*imag(I_sim_nat_baixa(1,k));
    T_frente(:,k) = (T_alin(:,k) + T_rel(:,k));
    P_frente(:,k) = T_frente(:,k)*ws;
    
    %Potencia nos enrolamentos.
    P_cu_estator_frente(:,k) = 2*(R*(abs(I_1n(:,k)))^2);
    T_cu_estator_frente(:,k) = P_cu_estator_frente(:,k)/ws;
    P_cu_estator_tras(:,k) = 2*(R*(abs(I_2n(:,k)))^2);
    T_cu_estator_tras(:,k) = P_cu_estator_tras(:,k)/ws;
        
    %Potencia pra tras.
    If_menos(k) = -V_2n(:,k)*Zeq_tras1;
    Ib_menos(k) = V_2n(:,k)*Zeq_tras2;
    
    T_f_menos(k) = 2*(1/(ws))*(R/(-3))*abs(If_menos(k))^2;
    T_b_menos(k) = 2*(1/(ws))*Re_neg*abs(Ib_menos(k))^2; 
    
    T_tras_tot(k) = T_f_menos(k) + T_b_menos(k);
    P_tras_tot(k) = ws*T_tras_tot(k);
    T_tras(k) = T_tras_tot(k) - T_cu_estator_tras(:,k);
    P_tras(k) = ws*T_tras(k);
    P_tras(k) = real(2*conj(I_2n(:,k))*V_2n(:,k)) - P_cu_estator_tras(:,k);
    T_tras(k) = P_tras(k)/ws;
     
    %Potencia no condutor do rotor.   
    V_entreferro_tras(:,k) = V_2n(:,k)- ((R1 + 1i*Xls)/(beta^2))*(I_2n(:,k));
    V_entreferro_q(:,k) = V_entreferro_tras(:,k)*(negativo_eixo_q/(negativo_eixo_d + negativo_eixo_q));
    V_entreferro_d(:,k) = V_entreferro_tras(:,k)*(negativo_eixo_d/(negativo_eixo_d + negativo_eixo_q));
    I_d_negativo(:,k) = (V_entreferro_d(:,k))/((0.5*(Rrd/2 + 1i*Xlrd ))/beta^2);
    I_q_negativo(:,k) = (V_entreferro_q(:,k))/((0.5*(Rrq/2 + 1i*Xlrq))/beta^2);
    Wcur1(:,k) = Rrd/(beta^2)*(abs(I_d_negativo(:,k)))^2 ;
    Wcur2(:,k) = Rrq/(beta^2)*(abs(I_q_negativo(:,k)))^2;
    Wcur(:,k) = Wcur1(:,k) +  Wcur2(:,k);
    T_Wcur(:,k) = Wcur(:,k)/ws;
           
    P_entreferro(:,k) = 2*real(conj(I_1n(:,k))*V_1n(:,k)) + 2*real(conj(I_2n(:,k))*V_2n(:,k)) - P_cu_estator_frente(:,k) - P_cu_estator_tras(:,k) - P_ferro(:,k);
    T_entreferro(:,k) = P_entreferro(:,k)/ws;
    P_tras_frente(:,k) = -P_entreferro(:,k) +P_frente(:,k);
    T_tras_frente(:,k) = P_tras_frente(:,k)/(ws);
   
    Eff(k) = P_entreferro(:,k)/P_entrada(:,k);
   
    T_frente (k) = T_rel(:,k) + T_alin(:,k);
                     
end

cd ..
cd interface
LSPM_Plot_Manager;

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

