%% Simulacao - LSPM

% Entradas
P = 1; % Par de polos
Nep = [1;1;1;1;1;1;1]; % Número de Espiras da i-ésima Ranhura
phi_rotor = [1;1;1;1;1;1;1]; % Angulo entre as barras do rotor
Ns = 10; % Numero de ranhuras no estator
Nr = 10; % Numero de ranhuras no rotor
Np = 10; % Numero enrolamento principal
Na = 10; % Numero enrolamento auxiliar
RIE = 1; % Raio interno do estator 
HRp = 1; % Altura media do enrolamento na ranhura
NC = 1; % numero de condutores
MIP = 1; % Fator multiplicativo
HC = 1; % Altura do colarinho
phi_p = 1; % diametro do fio
pho_cond = 1; % Peso especifico do condutor
pho_rotor = 1; % Resistividade do material condutor do rotor
CP = 1; % Comprimento de pacote
Sb = 1; % Area da secao transversal da barra
ra = 1; % raio medio dos aneis do rotor
Sa = 1; % Area da secao transversal do anel
ref = 1; % Coeficiente medido na temperatura de referencia
T = 1; % Temperatura real
CI = 1; % Circuitos internos
pho_enro = 1; % Resistividade do material do enrolamento
bs = 1; % Largura da ranhura
hp1 = 1; % Altura de pescoco de ranhura do estator
g = 1; % Entreferro geometrico
Xxm = 1; % Fator de correcao devido aos slots dos imas
bp1 = 1; % Largura do corpo de ranhura do estator
hs1 = 1; % Altura do corpo de ranhura do estator
RCMp = 1; % Distancia do eixo ate um ponto medio da bobina
RTC = 1; % Raio da area transversal da cabeca
btte = 1; % Largura da sapata do estator
hs2 = 1; % Altura do corpo da ranhura
bs2 = 1; % Largura do corpo de ranhura do rotor
Vd = 110; % Tensao induzida
Bsat = 1; % Inducao de saturacao
hpeq = 1; % Espessura equivalente do pescoco do rotor
alpha = 1; % angulo de fluxo gerado pelo ima
Ml =  1; %Largura do ima
Mh = 1; % altura do ima
Md = 1; % ??????
Lstk = 1; % comprimento do pacote (ima)
Perm = 1; % ??????
Sint = 1; % ??????
ur = 1; % permeabilidade relativa do ima
Ms_hgt = 1; % Altura dos slots dos imas
M_hgt = 1; % Altura dos imas
r_perc = 1; % percentual de perda de remanencia do ima com temperatura
Br = 1; % Inducao remanente do ima
phi_m1 = 1; % Fluxo fundamental devido ao ima
C = 3e-6; % Capacitor
E = 198; % Tensao a vazio
uo = pi*4e-7;
f = 50; % Frequencia
Vlin = 200; % Tensao de rede

%% Entradas Variaveis
delta_carga = 1;
s = 1;

%% Programa
Calculo_Geometrico
Calculo_Resistencias
Calculo_Indutancias
Calculo_Tensao_Induzida_Ima
Parte_Assincrona
Parte_Sincrona