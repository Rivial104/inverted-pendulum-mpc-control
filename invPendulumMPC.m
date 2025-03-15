clear all; close all; clc;

%% Parameters (Same as LQR Code)
r = 0.006; M = 0.135; m = 0.1; I = 0.0007176;
l = 0.2; g = 9.81; b = 0.00007892; c = 0.63;
L = 0.046; Rm = 12.5; kb = 0.031; kt = 0.031;
Er = 2*m*g*l; Ts = 0.01; % Sample time

%% State-Space Matrices
AA = I*(M+m) + M*m*(l^2);
aa = (((m*l)^2)*g)/AA; bb = ((I +m*(l^2))/AA)*(c + (kb*kt)/(Rm*(r^2)));
cc = (b*m*l)/AA; dd = (m*g*l*(M+m))/AA;
ee = ((m*l)/AA)*(c + (kb*kt)/(Rm*(r^2))); ff = ((M+m)*b)/AA;
mm = ((I +m*(l^2))*kt)/(AA*Rm*r); nn = (m*l*kt)/(AA*Rm*r);

A = [0 0 1 0; 0 0 0 1; 0 aa -bb -cc; 0 dd -ee -ff];
B = [0; 0; mm; nn];
C = eye(4); D = [0; 0; 0; 0];

%% Discretization (Zero-Order Hold)
sysd = c2d(ss(A, B, C, D), Ts, 'zoh');
Ad = sysd.A; Bd = sysd.B;

%% MPC Setup
N = 2; % Prediction Horizon
Nu = 1; % Control Horizon
Q = diag([1200 1500 0 0]); R = 0.035; % Cost Weights

% Build Prediction Matrices
[H, F] = build_mpc_matrices(Ad, Bd, Q, R, N, Nu);

%% Simulation
Tf = 10; X0 = [0.2; 179.5*pi/180; 0; 0]; % Initial state
X_des = [0; pi; 0; 0]; % Desired state
u_max = 12*10^20; u_min = -12*10^20; % Input constraints

X = X0; u_seq = [];
for k = 1:round(Tf/Ts)
    % Compute optimal control input
    x_ref = repmat(X_des, N, 1);
    f = F * (X - X_des);
    u_opt = quadprog(H, f, [], [], [], [], u_min, u_max);
    u = u_opt(1); % Apply only first control action
    
    % System Update
    X = Ad * X + Bd * u;
    Xp(k, :) = X'; u_seq(k) = u;
end

%% Plot Results
figure;
subplot(3,1,1); plot(Xp(:,1)); title('Cart Position'); grid on;
subplot(3,1,2); plot(Xp(:,2)*180/pi); title('Pendulum Angle (deg)'); grid on;
subplot(3,1,3); plot(u_seq); title('Control Input (Voltage)'); grid on;

function [H, F] = build_mpc_matrices(A, B, Q, R, N, Nu)
    nx = size(A,1); nu = size(B,2);
    Q_bar = kron(eye(N), Q);
    R_bar = kron(eye(Nu), R);
    
    Phi = []; Gamma = [];
    for i = 1:N
        Phi = [Phi; A^i];
        row = [];
        for j = 1:min(i,Nu)
            row = [A^(i-j)*B, row];
        end
        row = [row, zeros(nx, nu*(Nu-j))];
        Gamma = [Gamma; row];
    end
    
    H = Gamma' * Q_bar * Gamma + R_bar;
    F = Gamma' * Q_bar * Phi;
end