%% Alternativa para calcular Id1 e Iq1
V_1n(:,k) = -V1d + V1q;
h = Xd_baixa*Xq_baixa + (R2^2);
B = [R2/h, -1i*Xq_baixa/h; -1i*Xd_baixa/h, R2/h];
Inn = B*[V1d; V1q - 1i*E_baixa];

%% Alternativa torque total - Menos preciso
T_total = 2*(real(Z_1n)*power(abs(I_1n),2) - real(Z2_baixa)*power(abs(I_2n),2))/ws;
T_frente = 2*real(Z_1n)*power(abs(I_1n),2)/ws;
T_tras = 2*real(Z2_baixa)*power(abs(I_2n),2)/ws;

%% Potencia de entrada
Pelec(:,k) = 2*real(conj(I1(k))*V1(k)) + 2*real(conj(I2(k))*V2(k));

%% Torques
T_rel(:,k) = (p/ws)*(Xd-Xq)*real(In(1))*imag(In(2))/beta;
T_alin(:,k) = p*imag(In(2))*E/ws;

%% Potencia tras
T_tras(k) = T_tras_tot(k) - T_cu_estator_tras(:,k);
P_tras(k) = ws*T_tras(k);
P_tras(k) = real(conj(I2(:,k))*V2(:,k)) - P_cu_estator_tras(:,k);
T_tras(k) = P_tras(k)/ws;

%% Potencia Entreferro
P_entreferro(k) = Pelec(k) - P_cu_estator_frente(k) - P_cu_estator_tras(k) - P_ferro(k);
T_entreferro(k) = P_entreferro(k)/ws;
P_tras_frente(k) = -P_entreferro(k) + P_ps(k);
T_tras_frente(k) = P_tras_frente(k)/(ws);

%% Perdas no Rotor
V_entreferro_tras(:,k) = V2(:,k)- ((R1 + 1i*X1)/(beta^2))*(I2(:,k));
V_entreferro_q(:,k) = V_entreferro_tras(:,k)*(Zq/(Zd + Zq));
V_entreferro_d(:,k) = V_entreferro_tras(:,k)*(Zd/(Zd + Zq));
I_d_negativo(:,k) = (V_entreferro_d(:,k))/((0.5*(R2d/2 + 1i*X2d ))/beta^2);
I_q_negativo(:,k) = (V_entreferro_q(:,k))/((0.5*(R2q/2 + 1i*X2))/beta^2);
Wcur1(:,k) = R2d/(beta^2)*(abs(I_d_negativo(:,k)))^2 ;
Wcur2(:,k) = R2q/(beta^2)*(abs(I_q_negativo(:,k)))^2;
Wcur(:,k) = .5*(Wcur1(:,k) +  Wcur2(:,k));

% Sendo que
