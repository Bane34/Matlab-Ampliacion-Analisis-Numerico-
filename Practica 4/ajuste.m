%% (k, 1/k + 1/k^2 + 1/k^3)

function [A, b, aEN, aQR] = ajuste(k0)
    for k = 1 : 50
        A(k, 1) = 1 / (k0 + k);
        A(k, 2) = 1 / (k0 + k)^2;
        A(k, 3) = 1 / (k0 + k)^3;

        b(k, 1) = ((k0 + k)^2 + (k0 + k) + 1) / ((k0 + k)^3);
    end

    % Calculo del ajuste por minimos cuadrados
    B = A' * A;
    b_prim = A' * b;

    aEN = B \ b_prim;

    % Calculo mediante la factorización QR
    [Q, R] = qr(A);

    c = Q \ b;
    aQR = R \ c;   
end