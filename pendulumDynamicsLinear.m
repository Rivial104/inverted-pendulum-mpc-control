function [A_d, B_d] = pendulumDynamicsLinear(x0, u0, Tc)

    % Parametry systemu
    M = 0.5;
    m = 0.2;
    l = 2;
    g = 9.81;
    I = 10e-3;
    b = 0.4;      % tarcie w wózku
    b_pend = 0.8; % tarcie w przegubie wahadła

    % Stała D - wykorzystana przy wyznaczaniu ciągłych równań
    D = (I*(M+m) + M*m*l^2);

    % Obliczenie ciągłych macierzy stanu (A) i wejścia (B)
    % (wyznaczone w punkcie (x0, u0))
    A1 = -((I + m*l^2)*b)/D;
    A2 = ((m^2)*g*l^2)/D;
    A3 = -(m*l*b)/D;
    A4 = (m*g*l)*(M+m)/D;
    B1 = ((I + m*l^2)/D);
    B2 = ((m*l)/D);

    A = [0, 1, 0, 0;
         0, A1, A2, 0;
         0, 0, 0, 1;
         0, A3, A4, 0];
    
    B = [0; B1; 0; B2];
    
    % Dyskretyzacja metodą Eulera:
    A_d = eye(4) + Tc * A;
    B_d = Tc * B;
end
