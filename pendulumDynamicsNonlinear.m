function x_next = pendulumDynamicsDiscrete(x, u, Tc)
% pendulumDynamicsDiscrete - Dyskretny model układu wózek-wahadło
%
% Wejścia:
%   x  - wektor stanu: [x_cart; xdot_cart; theta; thetadot]
%   u  - sterowanie (siła przyłożona do wózka)
%   Tc - okres próbkowania
%
% Wyjście:
%   x_next - stan w następnym kroku: x(k+1) = x(k) + Tc * f(x(k), u(k))
%
% Parametry (ustawione wewnątrz funkcji, można je przekazywać jako argumenty):
    M = 0.5;      % masa wózka [kg]
    m = 0.2;      % masa wahadła [kg]
    l = 2;        % długość wahadła [m]
    g = 9.81;     % przyspieszenie ziemskie [m/s^2]
    I = 10e-3;    % moment bezwładności wahadła [kg*m^2]
    b = 0.4;      % współczynnik tarcia w wózku
    b_pend = 0.8; % współczynnik tarcia w przegubie wahadła

% Dla uproszczenia korzystamy z poniższego równania,
% gdzie przyjmujemy, że mianownik D jest dany przez:
    D = (I*(M+m) + M*m*l^2);

    dx = zeros(4,1);

    % Ciągłe równania stanu:
    % 1) Pozycja wózka
    dx(1) = x(2);
    % 2) Przyspieszenie wózka
    dx(2) = -(((I + m*l^2)*b)/D)*x(2) + (((m^2)*g*l^2)/D)*x(3) + ((I + m*l^2)/D)*u;
    % 3) Kąt wahadła
    dx(3) = b_pend * x(4);
    % 4) Przyspieszenie kątowe wahadła
    dx(4) = -((m*l*b)/D)*x(2) + ((m*g*l)*(M + m)/D)*x(3) + ((m*l)/D)*u;
    
    % Dyskretyzacja metodą Eulera:
    x_next = x + Tc * dx;
end

