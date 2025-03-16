%% MPC dla wahadła odwróconego na wózku z horyzontem sterowania i animacją trajektorii
clear; clc; close all;

%% Parametry modelu
M = 0.5;      % Masa wózka [kg]
m = 0.2;      % Masa wahadła [kg]
l = 2;        % Długość wahadła [m]
g = 9.81;     % Przyspieszenie ziemskie [m/s^2]
I = 10e-3;
b = 0.4;

Ts = 5*10e-3;      % Okres próbkowania [s]
N = 10;         % Horyzont predykcji (liczba kroków)
M = 4;          % Horyzont sterowania (N_u <= N)
Tsim = 10;      % Czas symulacji [s]
numSteps = round(Tsim/Ts);
time = 0:Ts:Tsim;


% Wagi w funkcji celu
Q = diag([20, 0, 10, 0]); % [x, v, theta, omega]
R_val = 0.01;              % Waga sterowania (skalarny współczynnik)

% MISO 
nx = 4;
nu = 1;  

C = [1 0 0 0;  
     0 0 1 0]; 
ny = size(C,1);

% Aby Rbar miało właściwe wymiary, definiujemy R jako macierz (nu x nu)
R = R_val * eye(nu);

% Ograniczenia na sterowanie
F_max = 20;   % [N]
F_min = -20;  % [N]

default_u = 0;

%% Stan początkowy i punkt odniesienia
x_current = [0; 0; 0.001; 0];   % np. wahadło początkowo odwrócone
% x_ref = [0; 0; 0; 0];       % chcemy dojść do [4; 0; pi; 0]
y_ref_vec = [2; 0];
y_ref = repmat(y_ref_vec, 1, N);

u_current = 0;

% u_min = -30;
% u_max = 30;

% Historia symulacji
x_history = zeros(nx, numSteps+1);
u_history = zeros(nu, numSteps);

x_history(:,1) = x_current;

delta = 1e-5;

%% Przygotowanie animacji
figure('Name','Animacja - trajektoria wahadła');
axis equal; grid on; hold on;
margin = 5; % margines widoku
% Początkowy widok ustalamy, ale będą dynamicznie zmieniane
xlim([-margin, margin]); ylim([-2, 8]);
cart_width = 2; cart_height = 1;
cart_y = 0;
pendulum_line = line([0,0],[0,0],'LineWidth',2,'Color','b');
cart_patch = patch('XData',[],'YData',[],'FaceColor','r');
pred_line = line(nan, nan, 'LineStyle','--','Color','g','LineWidth',1.5);

tip_traj = [];  % Trajektoria końca wahadła



%% MPC loop
for k = 1:numSteps
    nx = length(x_current);
    
   X_nom = zeros(nx,N+1);
   U_nom = zeros(1,N); 
   Y_nom = zeros(ny,N);

   X_nom(:,1) = x_current;
    for p = 1:N
        U_nom(:,p) = u_current;  
        X_nom(:,p+1) = pendulumDynamicsNonlinear(X_nom(:,p), U_nom(:,p), Ts);
    end

    for p = 1:N
        Y_nom(:,p) = C * X_nom(:,p); 
    end

    % Linearyzacja w punkcie równowagi
    A_list = cell(N,1);
    B_list = cell(N,1);
    
    for p = 1:N
     [A_list{p},B_list{p}] = pendulumDynamicsLinear(X_nom(:,p), U_nom(:,p), Ts);
    end  


    % Macierze predykcji
    phi = zeros(nx*N, nx);    % Extended A matrix
    gamma = zeros(nx*N,nu*M); % Extended B matrix

    x_prev = eye(nx);

     for p = 1:N
        % x(k+p) ~ A_{p-1} * x(k+p-1) + B_{p-1} * u(k+p-1)
        % ale A_{p-1}, B_{p-1} zależne od punktu linearyzacji p-1
        x_prev = A_list{p} * x_prev;  % kaskada
        rowStart = (p-1)*nx + 1;
        rowEnd   = p*nx;
        phi(rowStart:rowEnd, :) = x_prev;
        
        A_prod = eye(nx);
        for j = p:-1:1
            if j > M
                A_prod = A_prod * A_list{j-1};
                continue;
            end
            colStart = (j-1)*nu + 1;
            colEnd   = j*nu;
            gamma(rowStart:rowEnd, colStart:colEnd) = A_prod * B_list{j};
            if j > 1
                A_prod = A_prod * A_list{j-1};
            end
        end
    end
    
    r = (ny * N) / 4;

    C_ext = repmat(C', r, 1);

    M_k = kron(eye(N), C) * gamma;
    % M_k = C_ext .* gamma(:,1:2);
    % M_k = C_ext * gamma;

    % Macierze wag
    Psi = eye(ny);
    Lambda = eye(nu);

    Psi_bar = kron(eye(N), Psi);   
    Lambda_bar = kron(eye(M), Lambda);

    E_vector = reshape(y_ref - Y_nom, [], 1);

    H = M_k' * Psi_bar * M_k + Lambda_bar;
    f = - M_k' * Psi_bar * E_vector;

    % --- Formułowanie funkcji celu QP --- 
    % Qbar = kron(eye(N), Q);         % (N*nx x N*nx)
    % Rbar = kron(eye(N_u), R);         % (N_u*nu x N_u*nu)
    % H = G_bar' * Qbar * G_bar + Rbar;
    % H = (H + H')/2; % wymuszenie symetrii
    % f_qp = G_bar' * Qbar * (F_mat*x_current - repmat(x_ref, N, 1));

    % [H,f] = buildCostFunction()
    
    lb = repmat(F_min, M, 1);
    ub = repmat(F_max, M, 1);

    A_ineq = zeros(N, size(gamma,2));  % N nierówności, gamma ma (N*nx) wierszy i (N_u*nu) kolumn
    b_ineq = zeros(N,1);
    theta_max = pi/8;
    
    for p = 1:N
        rowTheta = (p-1)*nx + 3;  % zakładamy, że stan (x, dx, theta, dtheta) => theta jest 3
        A_ineq(p,:) = 1 * gamma(rowTheta, :); % C_theta * gamma()
        b_ineq(p) = theta_max - Y_nom(2,p);  % pi/6 - theta_nom(k+p)
    end
    
    opts = optimoptions('quadprog','Display','off');
    [U_opt, ~, exitflag] = quadprog(H, f, [], [], [], [], lb, ub, [], opts);
    
    if exitflag <= 0 || isempty(U_opt)
        warning('QP infeasible or not converged at step %d, applying default control.', k);
        u_opt = default_u; % np. default_u = zeros(nu,1)
    else
        u_opt = U_opt(1:nu);
    end
    
    disp(u_opt);
    
    % --- Aktualizacja stanu przy użyciu pełnego, nieliniowego modelu ---
    x_next = x_current + Ts * pendulumDynamicsNonlinear(x_current, u_opt, Ts);
    x_current = x_next;
    x_history(:,k+1) = x_current;
    u_history(:,k) = u_opt;
    
    % --- Animacja: rysowanie wózka i wahadła ---
    cart_x = x_current(1) - cart_width/2;
    cart_patch.XData = [cart_x, cart_x+cart_width, cart_x+cart_width, cart_x];
    cart_patch.YData = [cart_y, cart_y, cart_y+cart_height, cart_y+cart_height];
    
    pendulum_origin = [x_current(1), cart_y+cart_height];
    pendulum_end = pendulum_origin + l * [sin(x_current(3)), cos(x_current(3))];
    set(pendulum_line, 'XData', [pendulum_origin(1), pendulum_end(1)], ...
        'YData', [pendulum_origin(2), pendulum_end(2)]);
    
    % % --- Predykcja horyzontu --- (opcjonalnie)
    % x_pred = F_mat * x_current + G_bar * U_opt;
    % pred_tip_x = zeros(N,1);
    % pred_tip_y = zeros(N,1);
    % for i = 1:N
    %     x_i = x_pred((i-1)*nx+1:i*nx);
    %     pendulum_origin_pred = [x_i(1), cart_y+cart_height];
    %     pendulum_end_pred = pendulum_origin_pred + l * [sin(x_i(3)), cos(x_i(3))];
    %     pred_tip_x(i) = pendulum_end_pred(1);
    %     pred_tip_y(i) = pendulum_end_pred(2);
    % end
    % set(pred_line, 'XData', pred_tip_x, 'YData', pred_tip_y);
    
    % Dynamiczne przesuwanie widoku: środek widoku = pozycja wózka
    x_center = x_current(1);
    xlim([x_center - margin, x_center + margin]);
    % Ustawiamy stałe granice osi Y
    ylim([-2, 8]);
    
    drawnow;
    pause(0.01);
    
    tip_traj = [tip_traj, pendulum_end'];
end

%% Wykresy końcowe
figure;
subplot(2,1,1);
plot(time, x_history(1,:), 'LineWidth',1.5); hold on;
plot(time, zeros(size(time)), '--r','LineWidth',1.5);
xlabel('Czas [s]'); ylabel('x (pozycja wózka) [m]');
legend('x','Referencja'); grid on;

subplot(2,1,2);
plot(time, x_history(3,:), 'LineWidth',1.5); hold on;
plot(time, zeros(size(time)), '--r','LineWidth',1.5);
xlabel('Czas [s]'); ylabel('\theta (kąt wahadła) [rad]');
legend('\theta','Referencja'); grid on;

figure;
plot(time(1:end-1), u_history, 'LineWidth',1.5);
xlabel('Czas [s]'); ylabel('Siła F [N]');
title('Sterowanie'); grid on;

figure;
plot(tip_traj(1,:), tip_traj(2,:), 'b.-','LineWidth',1.5, 'MarkerSize',10);
xlabel('x [m]'); ylabel('y [m]');
title('Trajektoria końca wahadła w płaszczyźnie XY');
grid on; axis equal;

function y = pendulumOutput(x)
% pendulumOutput - Funkcja wyjścia układu wózek-wahadło
%
% Wyjście definiujemy jako:
%   y = [ x_cart; theta ]
%
% Wejście:
%   x - wektor stanu: [x_cart; xdot_cart; theta; thetadot]
%
% Wyjście:
%   y - [x_cart; theta]

    y = [x(1); x(3)];
end
