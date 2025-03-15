function [H, f] = buildCostFunction(Y_nom, Y_ref, M_k, Psi_bar, Lambda_bar)
    % Y_nom: wektor (ny*N,1) - zebrane nominalne wyjścia
    % Y_ref: wektor (ny*N,1) - zebrane zadane wyjścia
    % M_k: macierz (ny*N x nu*N_u) - Toeplitzowa zależność \DeltaY od \DeltaU
    % Psi_bar: (ny*N x ny*N) macierz wag
    % Lambda_bar: (nu*N_u x nu*N_u) macierz wag
    %
    % Zwracamy H, f dla quadprog: 0.5 * dU^T * H * dU + f^T * dU
    
    E = Y_ref - Y_nom;       % (ny*N,1) wektor błędu
    H = M_k' * Psi_bar * M_k + Lambda_bar;
    H = 0.5 * (H + H');      % symetryzacja
    f = - M_k' * Psi_bar * E; % UWAGA: w quadprog MATLAB jest postać 0.5*dU^T*H*dU + f^T*dU
end