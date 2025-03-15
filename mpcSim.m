clear; clc; close all;

%%  Parametry wahadła
m = 0.2;   % Masa wahadła [kg]
M = 0.5;   % Masa wózka [kg]
l = 0.3;   % Długość wahadła [m]
g = 9.81;  % Przyspieszenie ziemskie [m/s^2]
b = 0.1;   % Tarcie wózka

% Macierze stanu
A = [0 1 0 0;
     0 0 -b/M g/m;
     0 0 0 1;
     0 0 b/(l*M) -g/l];

B = [0; 1/M; 0; -1/(l*M)];
C = [1 0 0 0];
D = 0;

% Funkcja nieliniowa
f = @(x, u) A*x + B*u; % Bez nieliniowości dla porównania

%% Sprawdzenie trajektorii
N = 6;     % Horyzont predykcji
dt = 0.1;   % Krok czasowy
x0 = [0; 0; pi/6; 0];  % Stan początkowy (małe odchylenie)
U0 = zeros(N, 1); % Zerowe sterowanie początkowe

[X0, Y0, A_seq, B_seq, M_seq] = linearizeMPCdynamics(f, C, x0, U0, N, dt);

%% Wykres trajektorii
figure;
subplot(2,1,1);
plot(0:N, X0(1,:), 'b-o', 'DisplayName', 'Pozycja wózka');
hold on;
plot(0:N, X0(3,:), 'r-o', 'DisplayName', 'Kąt wahadła');
xlabel('Krok predykcji k'); ylabel('Stan');
legend; grid on;
title('Trajektoria stanu X_0');

subplot(2,1,2);
plot(0:N-1, Y0, 'g-o', 'DisplayName', 'Wyjście (pozycja wózka)');
xlabel('Krok predykcji k'); ylabel('Wyjście');
legend; grid on;
title('Trajektoria wyjść Y_0');

%% 🛠 Sprawdzenie macierzy A_seq i B_seq (powinny być bliskie A, B)
disp('Macierz A_seq(:,:,1) (powinna być bliska A):');
disp(A_seq(:,:,1));

disp('Macierz B_seq(:,:,1) (powinna być bliska B):');
disp(B_seq(:,:,1));
