clc;clear;

Tc = 0.1;       % okres próbkowania [s]
Tfinal = 1;      % czas symulacji [s]
time = 0:Tc:Tfinal;
N = length(time);

% Warunki początkowe
x0 = [1; 0; 10e-6; 0];  % mały kąt 0.1 rad
x_nonlinear = zeros(4, N);
x_linear    = zeros(4, N);
x_nonlinear(:,1) = x0;
x_linear(:,1)    = x0;
x_history = zeros(4, N);
x_compare = zeros(8, N);
u = 0;  

for k = 1:N-1
    
    % u = 3*sin(1/10*k); 

    % Model nieliniowy: obliczenie pochodnej i dyskretyzacja Eulera
    x_dot = pendulumDynamicsNonlinear(x_nonlinear(:,k), u);
    % disp(x_dot);
    x_history(:,k) = x_dot;
    x_nonlinear(:,k+1) = x_nonlinear(:,k) + Tc * x_dot;
    % x_nonlinear(2,k+1) = saturateVelocity(x_nonlinear(2,k+1));
    % x_nonlinear(3,k+1) = saturateAngle(x_nonlinear(3,k+1));
    
    % Model liniowy: dyskretyzacja przy danym punkcie pracy
    [A, B, A_d, B_d] = pendulumDynamicsLinear(Tc);
    x_linear(:,k+1) = A_d* x_linear(:,k) + B_d * u;
    % x_linear(2,k+1) = saturateVelocity(x_linear(2,k+1));
    % x_linear(3,k+1) = saturateAngle(x_linear(3,k+1));
end

x_compare(1:4,:) = x_nonlinear;
x_compare(5:8,:) = x_history;

figure;
subplot(2,1,1);
plot(time, x_nonlinear(1,:), 'LineWidth',1.5); hold on;
plot(time, x_linear(1,:), '--k', 'LineWidth',1.5);
xlabel('Czas [s]'); ylabel('x (pozycja wózka) [m]');
legend('Model nieliniowy', 'Model liniowy'); grid on;

subplot(2,1,2);
plot(time, x_nonlinear(3,:), 'LineWidth',1.5); hold on;
plot(time, x_linear(3,:), '--k', 'LineWidth',1.5);
xlabel('Czas [s]'); ylabel('\theta (kąt wahadła) [rad]');
legend('Model nieliniowy', 'Model liniowy'); grid on;

function theta_sat = saturateAngle(theta)
    % Definicja ograniczeń:
    maxTheta = pi/4;
    minTheta = -pi/4;
    
    % Ograniczenie na kąt:
    if theta > maxTheta
        theta_sat = maxTheta;
    elseif theta < minTheta
        theta_sat = minTheta;
    else
        theta_sat = theta;
    end
end

function vel_sat = saturateVelocity(V)
    % Definicja ograniczeń:
    maxV = 4;
    minV= -4;
    
    % Ograniczenie na kąt:
    if V > maxV
        vel_sat = maxV;
    elseif V < minV
        vel_sat = minV;
    else
        vel_sat = V;
    end
end
