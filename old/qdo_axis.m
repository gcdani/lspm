syms theta real

Ks = (2/3)*[cos(theta) cos(theta - 2*pi/3) cos(theta + 2*pi/3);
            sin(theta) sin(theta - 2*pi/3) sin(theta + 2*pi/3);
            0.5         0.5                         0.5];
        
invKs = [cos(theta)             sin(theta)              1;
         cos(tehta - 2*pi/3)    sin(theta - 2*pi/3)     1;
         cos(theta + 2*pi/3)    sin(theta + 2*pi/3)     1];
     