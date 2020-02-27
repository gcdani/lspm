%% Calculo das Resistencias
global CP Sb Sa Nr ra NC pho_rotor T Neep CI Rp Ra;

Rb = pho_rotor*CP/Sb; % Resistencia das barras
Rr = pho_rotor*pi*ra/Sa; % Resistencia dos aneis
Rrp = Rr*.5;

soma = 0;

for i = 1:Nr*.5
    soma = soma + power(sin(2*pi*i/Nr),2);
end

Rbp2 = Rb/soma; % Resistencia em paralelo de cada barra considerando o angulo entre elas
R2 = (2*Rbp2 + 2*Rrp)*power(Neep,2); % Resistencia total da gaiola referida pelo estator

R2t = R2*((ref + T)/(ref + 25)); % Corrigida pela temperatura

% Diferente
% Rp = Ltotp*pho_enrol/(A_fiop*power(NC*CI,2)); % Resistencia do enrolamento principal
Rp = Ltotp*2/(57*A_fiop*power(NC*CI,2));
Ra = Ltota*2/(57*A_fioa*power(NC*CI,2));