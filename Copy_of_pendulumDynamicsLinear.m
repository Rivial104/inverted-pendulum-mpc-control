function dx = pendulumDynamicsLinear(x, u, M, m, l, g, I)
    % x: [x; v; theta; omega]
    % u: siła przyłożona do wózka
    % Parametry tarcia:
    b = 0.4;   % współczynnik tarcia wózka (np. tarcie rolkowe)
    b_pend = 0.1; % współczynnik tarcia w przegubie wahadła
    
    dx = zeros(4,1);
    % Wyznaczamy D - efekt masy układu
    D = M + m - m*cos(x(3))^2;
    
    % Równania dynamiki z uwzględnieniem tarcia:
    dx(1) = x(2);
    % Dodatkowy składnik - tarcie wózka: -b_cart * v
    dx(2) = -((I+m*l^2)*b)/(I*(M+m) + M*m*l^2)*x(2) + ((m^2)*g*l^2)/(I*(M+m) + M*m*l^2)*u;
    dx(3) = x(4);
    % Dodatkowy składnik - tarcie przegubu: -b_pend * omega
    dx(4) = -(m*l*b)/(I*(M+m) + M*m*l^2)*x(2) + (m*g*l)*(M+m)/((I)*(M+m)+ M*m*l^2)*x(3);
end
