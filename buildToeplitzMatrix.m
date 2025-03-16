function [Phi_y, Gamma_y] = buildToeplitzMatrix(A, B, C, N)

    nx = size(A,1);
    nu = size(B,2);
    ny = size(C,1);

    Phi_y   = zeros(ny*N, nx);
    Gamma_y = zeros(ny*N, nu*N);

    for p = 1:N
        rowStart = (p-1)*ny + 1;
        rowEnd   = p*ny;

        % A^p
        A_p = A^p;   
        % C*A^p -> wpływ stanu początkowego na wyjście y(k+p)
        Phi_y(rowStart:rowEnd, :) = C * A_p;

        % Toeplitzowa część sterowań
        for j = 0:(p-1)
            colStart = j*nu + 1;
            colEnd   = (j+1)*nu;

            A_j = A^(p-1-j); 
            Gamma_y(rowStart:rowEnd, colStart:colEnd) = C * A_j * B;
        end
    end
end
