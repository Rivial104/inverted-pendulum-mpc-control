function [A, B] = pendulumDynamicsLinear(x0, u0, Tc)

    M = 0.5;
    m =0.2;
    l =2;
    g =9.81;
    I = 10e-3;
    b = 0.4;   % współczynnik tarcia wózka (np. tarcie rolkowe)
    b_pend = 0.1; % współczynnik tarcia w przegubie wahadła
    
    dx = zeros(4,1);

    D = (I*(M+m)+M*m*l^2);

    dx(1) = x0(2);
    dx(2) = -(((I+m*l^2)*b)/D)*x0(2) + (((m^2)*g*l^2)/D)*x0(3) + ((I+m*l^2)/D)*u0;
    dx(3) = b_pend*x0(4);
    dx(4) = -((m*l*b)/D)*x0(2) + ((m*g*l)*(M+m)/D)*x0(3) + ((m*l)/D)*u0;


    A1 = -((I+m*l^2)*b)/D;
    A2 = ((m^2)*g*l^2)/D;
    A3 = -(m*l*b)/D;
    A4 = (m*g*l)*(M+m)/D;
    B1 = ((I+m*l^2)/D);
    B2 = ((m*l)/D);

    A = [0, 1, 0, 0;
        0, A1, A2, 0;
        0, 0, 0, 1;
        0, A3, A4, 0];
    
    B = [0; B1; 0; B2];

end
