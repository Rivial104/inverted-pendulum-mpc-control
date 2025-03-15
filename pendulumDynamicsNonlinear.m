function dx = pendulumDynamicsNonlinear(x, u)
% PENDULUMDYNAMICSNONLINEAR - nieliniowy model układu wózek-wahadło
%
% Wejścia:
%   x = [ x; xdot; theta; thetadot ] - stan
%   u - sterowanie (siła przyłożona do wózka)
%   M - masa wózka
%   m - masa wahadła
%   l - długość wahadła (od osi obrotu do środka masy)
%   g - przyspieszenie ziemskie
%
% Wyjście:
%   dx = [ xdot; xddot; thetadot; thetaddot ] - pochodna stanu

    M = 0.5;
    m =0.2;
    l =2;
    g =9.81;
    I = 10e-3;
    b = 0.4;
    b_pend = 0.1;

    dx = zeros(4,1);
    u = u(1);
    D = (I*(M+m)+M*m*l^2);

    % Stan
    x_cart     = x(1);  % położenie wózka (nie jest bezpośrednio potrzebne w równaniach)
    xdot_cart  = x(2);  % prędkość wózka
    theta      = x(3);  % kąt wahadła
    thetadot   = x(4);  % prędkość kątowa wahadła

    % Równania stanu
    dx(1) = xdot_cart;
    dx(2) = -(((I+m*l^2)*b)/D)*xdot_cart + (((m^2)*g*l^2)/D)*theta + ((I+m*l^2)/D)*u;
    dx(3) = b_pend*thetadot;
    dx(4) = -((m*l*b)/D)*xdot_cart + ((m*g*l)*(M+m)/D)*theta + ((m*l)/D)*u;
end
