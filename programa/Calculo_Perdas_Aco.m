%% Perdas no Aco

%% Entrada
P = 1; % Par de polos
CP = 43.5e-3;
re = 65e-3; % Raio interno estator
ri = 42.44e-3; % Raio externo estator
slot_hgt = 11.6e-3; % Altura das ranhuras
fluxo_air_fund = 1.1; % Fluxo fundamental do entreferro
alpha = 10*pi/180; % Angulo de abertura do ima
pho_mat = 8.89; % Densidade do material
f = 50;
Area_ranhura = 1;
Lmre = 6e-3; % Largura media da ranhura do estator
Ns = 24; % Numero de ranhuras do estator
Cfch = .0205744;
Cfa = -1.0848123;
Cfb = 1.1444818;
Cfce = 0.000018;

%% Perdas na coroa
Vol_s = pi*CP*(re*re - ri*ri);  % Volume do estator
Vol_coroa_s = Vol_s - pi*CP*power(ri + slot_hgt,2) + 2*pi*CP*ri*ri; % Volume da coroa do estator
Bpk_coroa = fluxo_air_fund/(2*P*CP*(ri - re - slot_hgt)); % Inducao de pico na coroa
W_parasita_coroa = Cfce*2*pi*Bpk_coroa*f*f;    % Perdas na coroa do estator
W_hist_coroa = Cfch*f*power(Bpk_coroa, Cfa + Cfb*Bpk_coroa);  % 
W_coroa = pho_mat*Vol_coroa_s*alpha*(W_parasita_coroa + W_hist_coroa);

%% Perdas no dente
Vol_dente = pi*CP*(power(ri + max(slot_hgt),2) - ri*ri) - 4*CP*Area_ranhura;
Bpk_dente = fluxo_air_fund*P*2/(CP*Lmre*Ns); % Inducao de pico na coroa
W_parasita_dente = Cfce*2*pi*Bpk_dente*f*f;    % Perdas na coroa do estator
W_hist_dente = Cfch*f*power(Bpk_dente, Cfa + Cfb*Bpk_dente);  % 
W_dente = pho_mat*Vol_dente*(W_parasita_dente + W_hist_dente);

