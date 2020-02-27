%%% Programa para calculo dos valores de eficiencia,
%%% Angulo de carga do motor estudado
%%% Dados de entrada extraidos do SPEED - 29/Jul/2016

clear;
clc;
disp('Come?ou !!!');

format long;

delta = 0:90;     %input('Entre com o valor do angulo de carga: ')
%delta = 45;

%%%% Valor hipotetico para angulo de carga para comparar resultados do teste.

p = 1;
freqRede = 50;                        %input('Entre com a frequencia da rede: ')
beta = 1.7836889285222159;             %input('Entre com a rela?ao de transformacao entre os enrolamentos
%principal e auxiliar: ')
s = 2;

%%%% Fatores de ajustes do SPEED
a_XL = 0.91;     % usado no SPEED para multiplicar todas as indut?ncias

%%%% Valores do circuito equivalente tirados do SPEED.
P_ferro = 3.0345;
Vm_rms = 220;                         %input('Entre com o valor de tensao da fonte de alimenta?ao: ');
V_baixa = Vm_rms/beta;
Cap = 4.75e-6;                        %input('Entre com o valor do capacitor:');
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
Xmd = 2*pi*freqRede*Lmd;
Xmq = 2*pi*freqRede*Lmq;
E_baixa = E/beta;
R = R1/beta^2;
Xd = 2*pi*freqRede*Ld;
Xq = 2*pi*freqRede*Lq;
ZCap = 1/(1i*2*pi*Cap*freqRede);
ws = 2*pi*freqRede;

Xcap = 2*pi*Cap*freqRede;

negativo_eixo_d = 0.5*1i*Xmd*(Rrd/2 + 1i*Xlrd)/(1i*Xmd + Rrd/2 + 1i*Xlrd);
negativo_eixo_q = 0.5*1i*Xmq*(Rrq/2 + 1i*Xlrq)/(1i*Xmq + Rrq/2 + 1i*Xlrq);
Z2_alta = R1 + 1i*Xls + 1.0*(negativo_eixo_d + negativo_eixo_q);
Xd_baixa = Xd/beta^2;
Xq_baixa = Xq/beta^2;

Z2_baixa = Z2_alta/beta^2;

Zmd_menos = 1/(1/(1i*Xmd)+2/(Rrd +1i*2*Xlrd));
Zmq_menos = 1/(1/(1i*Xmq)+2/(Rrq +1i*2*Xlrq));
Zdf_menos = Zmd_menos -R/3 + 1i*Xls;
Zdb_menos = Zmd_menos - R + 1i*Xls;
Zqf_menos = Zmq_menos -R/3 + 1i*Xls;
Zqb_menos = Zmq_menos -R + 1i*Xls;
Zeq_tras1 = (Zdb_menos - Zqb_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);
Zeq_tras2 =(Zdf_menos + Zqf_menos)/(Zdf_menos*Zqb_menos + Zdb_menos*Zqf_menos);
Re_neg = 0.5*real(-Zmd_menos -Zmq_menos + ((Zmd_menos - Zmq_menos)^2)/(Zqf_menos + Zdf_menos));

a1 = [1 1; -1i 1i];
a2 = [1 1i; 1 -1i];
a3 = [beta beta;1i -1i];
a4 = [1/beta 1/beta;1i -1i];


tol = 0.00001;

%%I_1n = zeros(1,lenght(delta)); %Kleyton sugest?o

for k = 1:length(delta)
    
    V1_lin = V_baixa;
    Vd_1 = -V1_lin*sin(delta(k)*pi/180);
    Vq_1 = V1_lin*cos(delta(k)*pi/180);
    V1_lin = Vd_1 + 1i*Vq_1;
    
    for int = 1:10000
        Vd_1 = real(V1_lin);
        Vq_1 = imag(V1_lin);
        Iq_1 = -Vd_1/(1i*Xq_baixa);
        Id_1 = (E_baixa - Vq_1)/(1i*Xd_baixa);
        I1_lin = Id_1 + 1i*Iq_1;
        Z1 = V1_lin/I1_lin;
        a1n = 1 - (1i*Xcap/Z1);
        a2n = 1 - (1i*Xcap/Z2_baixa);
        V1_teste = (beta + 1i*a2n)*sqrt(2)*V_baixa/(beta*(a1n + a2n));
        
%         Rc = 2*power(E_baixa,2)/P_ferro;
%         Ic = sqrt(P_ferro/Rc);
%         Ied = (R*Rc*Rc*Vd_1 + (R + Rc)*Rc*Xq*Vq_1 - (R + Rc)*Xq*E_baixa)/(R*R*Rc*Rc + power(R + Rc,2)*Xd*Xq);
%         Ieq = (R*Rc*Rc*Vq_1 - (R + Rc)*Rc*Xq*Vd_1 - power(R + Rc,2)*R*Rc*E_baixa)/(R*R*Rc*Rc + power(R + Rc,2)*Xd*Xq);
%         Id_1 = Ied + Ic;
%         Iq_1 = Ieq + Ic;
%         I1 = Id_1 + 1i*Iq_1;
%         Z1 = V1_lin/I1;
%         a1n = 1 - (1i*ZCap/Z1);
%         a2n = 1 - (1i*ZCap/Z2_baixa);
%         V1_teste = (beta + 1i*a2n)*sqrt(2)*Vm_rms/(beta*(a1n + a2n));
        
        
        if abs(V1_teste - V1_lin) > tol
            V1_lin = V1_teste;
        else
            break
        end
    end
    
    T_rel(k) = (p/ws)*(Xd_baixa-Xq_baixa)*Id_1*Iq_1;
    T_alin(k) = p*Iq_1*E_baixa/ws;
    T_frente(k) = (T_alin(k) + T_rel(k));
    
end