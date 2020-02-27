function [V_sp, V_sn, Vo] = componentes_simetricas(Va,Vb,Vc)

% Calcula as componentes simetricas da entrada da funcao

a = exp(1i*2*pi/3);
invA = (1/3).*[1 1 1;1 a a*a; 1 a*a a];
V012 = invA*[Va; Vb; Vc];
V_sp = V012(2);
V_sn = V012(3);
Vo = V012(1);

end