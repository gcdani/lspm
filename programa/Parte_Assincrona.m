%% Circuito Assincrono
global CP raiz Rot_sinc Tind Tmag Tf Tt Tqdi Tqqi Vid Viq Tmotor R2ini R2F X1di X2di
warning off

Vd = Vlin;
V1 = [Vd; 0; 0; -Vq];

int = 0;

for Rot_assinc = 0:Rot_sinc
    
    int = int + 1;
    MIMP = [Z1d + Zmd, -Zmd, 0, 0;
            -Zmd, Zmd + Z2d(int) + RR(int), -RR(int), 0;
            0, -RR(int), RR(int) + Z2q(int) + Zmq, -Zmq;
            0, 0, -Zmq, Zmq + Z1q];
    
    Ia = MIMP\V1;
    I1d(int) = Ia(1);
    I2d(int) = Ia(2);
    I2q(int) = Ia(3);
    I1q(int) = Ia(4)*Neep/(Neea*1i);
    IT(int) = I1d(int) + I1q(int);
    IRR(int) = I2d(int) - I2q(int);
    
    AngI1d(int) = atan(imag(I1d(int))/real(I1d(int)))*180/pi;
    AngI1q(int) = atan(imag(I1q(int))/real(I1q(int)))*180/pi;
    
    Tf(int) = IRR(int)*IRR(int)*.25*R2F(int)/(pi*f);
    Tt(int) = (-power(IRR(int),2)*.5*R2T(int) + (power(I2d(int),2) + power(I2q(int),2))*R2T(int))*.5/(pi*f);
    Tind(int) = abs(Tf(int)) - abs(Tt(int));
    
    freq(int) = Rot_assinc/60;
    
    XCapi(int) = (1/(2*pi*freq(int)*Cap*1i))*power(Neep/Neea,2);
    fluxo_magneto = (CP*6.5e-4)/(43.5e-3);
    R2Ti = R2ini*.5;
    X1di(int) = X1d*freq(int)/f;
    X1qi(int) = X1q*freq(int)/f;
    X2di(int) = X2d*freq(int)/f;
    X2qi(int) = X2q*freq(int)/f;
    Xmdi(int) = Xmd*freq(int)/f;
    Xmqi(int) = Xmdi(int);
    
    MIMPi = [R1d + X1di(int) + Xmdi(int), -Xmdi(int), 0;
            -Xmdi(int), Xmdi(int) + X2di(int) + 2*R2Ti + X2qi(int) + Xmqi(int), -Xmqi(int);
            0, -Xmqi(int), Xmqi(int) + XCapi(int) + R1q + X1qi(int)];
    
    Vi = [2*pi*freq(int)*Neep*fluxo_magneto; 0; -2*pi*freq(int)*Neep*fluxo_magneto];
    Vid(int) = abs(Vi(1));
    Viq(int) = abs(Vi(3))*Neea/Neep;
    Iai = MIMPi\Vi;
    I1di(int) = 2*Iai(1);
    I1qi(int) = 2*Iai(3)*Neep/(Neea*1i);
    Vcapi(int) = I1qi(int)*XCapi(int)*power(Neea/Neep,2);
    
    Tqdi(int) = power(I1di(int),2)*R1d*.5/(pi*freq(int));
    Tqqi(int) = power(I1qi(int),2)*R1q*.5*power(Neea/Neep,2)/(pi*freq(int));
    Tmag(int) = abs(Tqdi(int)) + abs(Tqqi(int));
    
    Tmotor(int) = Tind(int) - Tmag(int);
    
end

warning on

cd(raiz)
cd interface
LSPM_Plot_Manager_Assincrono;