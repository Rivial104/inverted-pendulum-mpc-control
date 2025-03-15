function [X0, Y0, A_seq, B_seq, M_seq] = linearizeMPCdynamics(f, C, x0, U0, N, dt)
    % Wejścia:
    % f   - Funkcja nieliniowej dynamiki (anonimowa @(x,u))
    % C   - Macierz wyjść
    % x0  - Początkowy stan
    % U0  - Początkowe sterowania na horyzoncie [Nu x 1]
    % N   - Horyzont predykcji
    % dt  - Krok czasowy
    
    nx = length(x0);  % Liczba stanów
    nu = length(U0);  % Liczba wejść
    
    % Inicjalizacja trajektorii
    X0 = zeros(nx, N+1);
    Y0 = zeros(size(C,1), N);
    A_seq = zeros(nx, nx, N);
    B_seq = zeros(nx, nu, N);
    M_seq = zeros(size(C,1), nu, N);
    
    % Symulacja układu nieliniowego w otwartej pętli
    X0(:,1) = x0;
    for k = 1:N
        X0(:,k+1) = f(X0(:,k), U0(k));  % Symulacja kroku
        Y0(:,k) = C * X0(:,k);  % Wyjście
    end
    
    % Linearyzacja wokół trajektorii X0, U0
    epsilon = 1e-5;
    for k = 1:N
        A_k = zeros(nx, nx);
        B_k = zeros(nx, nu);
        
        % Aproksymacja różnicowa dla A_k
        for i = 1:nx
            dx = zeros(nx,1);
            dx(i) = epsilon;
            A_k(:,i) = (f(X0(:,k) + dx, U0(k)) - f(X0(:,k), U0(k))) / epsilon;
        end
        
        % Aproksymacja różnicowa dla B_k
        for j = 1:nu
            du = zeros(nu,1);
            du(j) = epsilon;
            B_k(:,j) = (f(X0(:,k), U0(k) + du) - f(X0(:,k), U0(k))) / epsilon;
        end
        
        % Zapisanie macierzy
        A_seq(:,:,k) = A_k;
        B_seq(:,:,k) = B_k;
        M_seq(:,:,k) = C * A_k;  % Macierz M(k)
    end
end
