%% MPC dla wahadła odwróconego na wózku z horyzontem sterowania i animacją trajektorii
clear; clc; close all;

%% Parametry modelu
M = 0.5;      % Masa wózka [kg]
m = 0.2;      % Masa wahadła [kg]
l = 1.5;        % Długość wahadła [m]
g = 9.81;     % Przyspieszenie ziemskie [m/s^2]
I = 10e-3;
b = 5;        % współczynnik tarcia w wózku
b_pend = 0.8; % współczynnik tarcia w przegubie

Ts = 10e-3;      % Okres próbkowania [s]
N = 20;         % Horyzont predykcji (liczba kroków)
M = 3;          % Horyzont sterowania (N_u <= N)
Tsim = 10;      % Czas symulacji [s]
numSteps = round(Tsim/Ts);
time = 0:Ts:Tsim;

% MISO 
nx = 4;
nu = 1;  
ny = 2;

X_nom = zeros(nx,N+1);
U_nom = zeros(1,N); 
Y_nom = zeros(ny,N);

x_err_k = zeros(nx,1);
u_err_k = zeros(nu,1);

% Wagi w funkcji celu
Q = diag([10, 0, 100, 0]); % [x, v, theta, omega]
R = 0.01;              % Waga sterowania (skalarny współczynnik)

% Aby Rbar miało właściwe wymiary, definiujemy R jako macierz (nu x nu)
% R = R_val * eye(nu);

% Ograniczenia na sterowanie
F_max = 20;   % [N]
F_min = -20;  % [N]

default_u = 0;

%% Stan początkowy i punkt odniesienia
x_current = [0; 2; 0.3; 0];   % wahadło początkowo odwrócone
x_ref = [0.1; 0; 0.1; 0];   
u_ref = 0.2; 
y_ref_vec = [0.1; 0];
y_ref = repmat(y_ref_vec, 1, N);

u_current = 10e-3;

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
xlim([-margin, margin]); ylim([-2, 8]);
cart_width = 1; cart_height = 0.5;
cart_y = 0;
pendulum_line = line([0,0],[0,0],'LineWidth',2,'Color','b');
cart_patch = patch('XData',[],'YData',[],'FaceColor','r');
pred_line = line(nan, nan, 'LineStyle','--','Color','g','LineWidth',1.5);

tip_traj = [];  % Trajektoria końca wahadła


%% MPC loop single-shooting

for k = 1:numSteps
    nx = length(x_current);

   X_nom(:,1) = x_current;

   disp("Xnom");
   disp(X_nom(:,1));

   C = [1 0 0 0;  
         0 0 1 0]; 
    ny = size(C,1);

    for p = 1:N
        U_nom(:,p) = u_current;  

        x_dot = pendulumDynamicsNonlinear(X_nom(:,p), U_nom(:,p));

        X_nom(:,p+1) = X_nom(:,p) + Ts * x_dot;
        Y_nom(:,p)   = C * X_nom(:,p);

        % Linearyzacja w punkcie równowagi
        [A_analitycal,B_analitycal,A,B] = pendulumDynamicsLinear(Ts);
    end

    % x_k = X_nom(:,k);
    % u_k = U_nom(:,k);

    % Linearyzacja w punkcie równowagi
    % [A,B] = pendulumDynamicsLinear(x_k, u_k, Ts);

    % A_x = zeros((N+1)*nx, nx);
    % B_u = zeros((N+1)*nx, N*nu);

    A_x = A;
    B_u = B;


    % A_X
    for i = 0:N
        A_x(i*nx+1:(i+1)*nx, :) = A^i;
    end

    % B_U
    for i = 2:(N+1)
        for j = 1:(i-1)
            % Blok wiersza i, kolumny j:
            % Powinien być równy A^(i-j-1)*B
            rowStart = (i-1)*nx + 1;
            rowEnd   = i*nx;
            colStart = (j-1)*nu + 1;
            colEnd   = j*nu;
            B_u(rowStart:rowEnd, colStart:colEnd) = A^(i-j-1) * B;
        end
    end

   % Constraints definition
   F_x = zeros(8,4);
   F_x(1:4,:) = eye(4);
   F_x(5:8,:) = -eye(4);
   F_x(1,1) = 0;

   g_x = [0;5;pi/6;pi/12;0;5;pi/6;pi/12];

   F_u = [1;-1];
   g_u = [20;-20];

   F_state = kron(eye(N), F_x);  
   F_input = kron(eye(N), F_u); 

    F_state_ext = [F_state, zeros(8*N, N)];  
    F_input_ext = [zeros(2*N, 4*N), F_input];
    F = [F_state_ext; F_input_ext];

    g = repmat([g_x; g_u], N, 1);

    % Define state and control reference in horizon
    x_ref_k = repmat(x_ref,N);
    u_ref_k = repmat(u_ref, 1, N);

    x_err = X_nom(:,1:N) - x_ref_k;
    u_err = U_nom - u_ref_k; 

    % Macierze wag
    Psi = kron(eye(N), Q);
    Lambda = kron(eye(N), R);

    % U_exp = [U_nom,0];
    % X_short = X_nom(:,1:N);

    % z = [X_short;U_nom];
    z = [x_err;u_err];
    z_vec = reshape(z,[],1);

    H = blkdiag(Psi, Lambda);

    % TODO - poprawić ub i lb
    ub = repmat([g_x(1:4); g_u(1)], N, 1);
    lb = repmat([g_x(5:8); g_u(2)], N, 1);


    opts = optimoptions('quadprog','Display','off');
    [U_opt, ~, exitflag] = quadprog(H, z, [], [], [], [], lb, ub, [], opts);

    if exitflag <= 0 || isempty(U_opt)
        warning('QP infeasible or not converged at step %d, applying default control.', k);
        u_opt = default_u; % np. default_u = zeros(nu,1)
    else
        u_opt = U_opt(2)+u_ref;
    end

    disp(u_opt);


    % --- Aktualizacja stanu przy użyciu pełnego, nieliniowego modelu ---
    x_next = x_current + Ts * pendulumDynamicsNonlinear(x_current, u_opt);
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

