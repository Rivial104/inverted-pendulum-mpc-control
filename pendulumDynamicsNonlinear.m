function x_dot = pendulumDynamicsNonlinear(x, u, Tc)
% Parametry układu
    M      = 0.5;    % masa wózka [kg]
    m      = 0.2;    % masa wahadła [kg]
    l      = 2;    % długość wahadła [m]
    g      = 9.81;   
    I      = 10e-3;  % moment bezwładności wahadła 
    b      = 5;      % współczynnik tarcia w wózku
    b_pend = 0.8;    % współczynnik tarcia w przegubie 

% Ekstrakcja zmiennych stanu:
    x_cart   = x(1);
    x_dot    = x(2);
    theta    = x(3);
    theta_dot= x(4);

% Definicja mianownika D 
    D = M + m - m * cos(theta)^2;

% Równanie na przyspieszenie wózka (x_ddot)
    x_ddot = ( u - b * x_dot + m * l * sin(theta) * theta_dot^2 - m * g * sin(theta) * cos(theta) ) / D;

% Równanie na przyspieszenie kątowe wahadła (theta_ddot)
    theta_ddot = ( -u * cos(theta) - m * l * theta_dot^2 * sin(theta) * cos(theta) + (M + m) * g * sin(theta) - b_pend * theta_dot ) / (l * D);

% Złożenie pochodnych stanu:
    dx = zeros(4,1);
    dx(1) = x_dot;
    dx(2) = x_ddot;
    dx(3) = theta_dot;
    dx(4) = theta_ddot;

    x_dot = dx;
end
