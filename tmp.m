%% MPC dla wahadła odwróconego na wózku z animacją trajektorii
clear; clc; close all;

%% Parametry modelu
M = 0.5;      % Masa wózka [kg]
m = 0.2;      % Masa wahadła [kg]
l = 0.3;      % Długość wahadła [m]
g = 9.81;     % Przyspieszenie ziemskie [m/s^2]

Ts = 0.05;    % Okres próbkowania [s]
N = 10;       % Horyzont predykcji (liczba kroków)
Tsim = 70;     % Czas symulacji [s]
time = 0:Ts:Tsim;

% Wagi w funkcji celu
Q = diag([10, 1, 100, 1]); % [x, v, theta, omega]
R = 0.01;                  % Waga sterowania

% Ograniczenia na sterowanie
F_max = 10;   % [N]
F_min = -10;  % [N]

%% Stan początkowy i punkt odniesienia
% Zakładamy, że chcemy utrzymać wózek w spoczynku i wahadło pionowo (theta=0)
x_current = [2; 0; 0.2; 0];   % niewielkie zaburzenie kąta (0.2 rad)
x_ref = [0; 0; 0; 0];

% Historia symulacji
x_history = x_current;
u_history = [];

% Parametr do linearizacji (finite differences)
delta = 1e-5;

%% Przygotowanie animacji
figure('Name','Trajektoria wahadła');
axis equal; grid on; hold on;
xlim([-1 1]); ylim([-0.5 1.5]);
cart_width = 0.3; cart_height = 0.2;
pendulum_line = line([0,0],[0,0],'LineWidth',2,'Color','b');
cart_patch = patch('XData',[],'YData',[],'FaceColor','r');

%% MPC loop
for k = 1:length(time)-1
    % --- Linearizacja przy użyciu różnic skończonych ---
    nx = length(x_current);
    nu = 1;
    A = zeros(nx,nx);
    B = zeros(nx,nu);
    % Obliczamy f(x,u) dla aktualnego stanu i sterowania (u0 = 0 przy linearizacji)
    u0 = 0;
    f0 = pendulumDynamics(x_current, u0, M, m, l, g);
    % Linearizacja względem stanu
    for j = 1:nx
        dx = zeros(nx,1);
        dx(j) = delta;
        f_plus = pendulumDynamics(x_current+dx, u0, M, m, l, g);
        A(:,j) = (f_plus - f0)/delta;
    end
    % Linearizacja względem sterowania
    f_u = pendulumDynamics(x_current, u0+delta, M, m, l, g);
    B(:,1) = (f_u - f0)/delta;
    
    % Dyskretyzacja (Euler)
    A_d = eye(nx) + Ts * A;
    B_d = Ts * B;
    
    % --- Budowa macierzy predykcyjnych (F i G) ---
    F_mat = zeros(N*nx, nx);
    G_mat = zeros(N*nx, N*nu);
    A_power = eye(nx);
    for i = 1:N
        A_power = A_power * A_d;
        F_mat((i-1)*nx+1:i*nx, :) = A_power;
        for j = 1:i
            A_power_j = eye(nx);
            for p = 1:(i-j)
                A_power_j = A_power_j * A_d;
            end
            G_mat((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = A_power_j * B_d;
        end
    end
    
    % --- Formułowanie problemu QP ---
    Qbar = kron(eye(N), Q);
    Rbar = kron(eye(N), R);
    H = G_mat' * Qbar * G_mat + Rbar;
    H = (H + H')/2;  % Wymuszamy symetrię
    f_qp = G_mat' * Qbar * (F_mat*x_current - repmat(x_ref, N, 1));
    
    % Ograniczenia na sterowanie przez horyzont predykcji
    lb = repmat(F_min, N, 1);
    ub = repmat(F_max, N, 1);
    
    % Rozwiązanie QP
    opts = optimoptions('quadprog','Display','off');
    U_opt = quadprog(H, f_qp, [], [], [], [], lb, ub, [], opts);
    u_opt = U_opt(1);  % tylko pierwsze sterowanie
    
    % --- Aktualizacja stanu za pomocą pełnego, nieliniowego modelu ---
    dx = pendulumDynamics(x_current, u_opt, M, m, l, g);
    x_next = x_current + Ts * dx;
    
    x_current = x_next;
    x_history = [x_history, x_current];
    u_history = [u_history, u_opt];
    
    %% Animacja – rysowanie wózka i wahadła
    % Pozycja wózka (prostokąt)
    cart_x = x_current(1) - cart_width/2;
    cart_y = 0;  % zakładamy, że wózek porusza się po poziomej linii
    cart_patch.XData = [cart_x, cart_x+cart_width, cart_x+cart_width, cart_x];
    cart_patch.YData = [cart_y, cart_y, cart_y+cart_height, cart_y+cart_height];
    
    % Pozycja końca wahadła (punkt na wahadle)
    pendulum_origin = [x_current(1), cart_y+cart_height];
    pendulum_end = pendulum_origin + l * [sin(x_current(3)), cos(x_current(3))];
    set(pendulum_line, 'XData', [pendulum_origin(1), pendulum_end(1)],...
        'YData', [pendulum_origin(2), pendulum_end(2)]);
    
    drawnow;
    pause(0.01);
end

%% --- Wykresy końcowe ---
figure;
subplot(2,1,1);
plot(time, x_history(1,:), 'LineWidth',1.5); hold on;
plot(time, x_ref(1)*ones(size(time)), '--r','LineWidth',1.5);
xlabel('Czas [s]'); ylabel('x (pozycja wózka) [m]');
legend('x','Referencja'); grid on;

subplot(2,1,2);
plot(time, x_history(3,:), 'LineWidth',1.5); hold on;
plot(time, x_ref(3)*ones(size(time)), '--r','LineWidth',1.5);
xlabel('Czas [s]'); ylabel('\theta (kąt wahadła) [rad]');
legend('\theta','Referencja'); grid on;

figure;
plot(time(1:end-1), u_history, 'LineWidth',1.5);
xlabel('Czas [s]'); ylabel('Siła F [N]');
title('Sterowanie'); grid on;

function dx = pendulumDynamics(x, u, M, m, l, g)
    % x: [x; v; theta; omega]
    % u: siła przyłożona do wózka
    dx = zeros(4,1);
    % przyjmujemy, że theta=0 jest pozycją pionową (upright)
    D = M + m - m*cos(x(3))^2;
    dx(1) = x(2);
    dx(2) = ( u + m*l*x(4)^2*sin(x(3)) - m*g*sin(x(3))*cos(x(3)) ) / D;
    dx(3) = x(4);
    dx(4) = ( - u*cos(x(3)) - m*l*x(4)^2*sin(x(3))*cos(x(3)) + (M+m)*g*sin(x(3)) ) / (l*D);
end
