function [d, l, u] = facLDLtcasitri(a, b, c)
    % Programa que calcula la factorización LDLt para una matriz
    % casitridiagonal simétrica donde:
    %   El vector a contiene la diagonal principal
    %   El vector b contiene la subdiagonal principal (que también es la supradiagonal principal de la matriz)
    %   El vector c contiene la última fila (que también es la última columna de la matriz)
    n = numel(a); %dimension de la matriz

    % Primera fila
    d(1,1) = a(1);

    % Filas 2 a n-1
    for j = 2 : n - 1
        l(j - 1, 1) = b(j - 1) / d(j - 1, 1);
        d(j, 1) = a(j) - d(j - 1, 1) * l(j - 1, 1)^2 ;
    end

    % Última fila, primera columna
    u(1, 1)= c(1, 1) / d(1, 1);

    % Última fila, columnas 2 a n-2
    for j = 2 : n - 2
        u(j, 1) = (c(j, 1) - u(j - 1) * l(j - 1, 1) * d(j - 1, 1)) / d(j, 1);
    end

    % Última fila, últimos dos elementos
    l(n - 1, 1) = (b(n - 1, 1) - u(n - 2, 1) * l(n - 2, 1) * d(n - 2, 1)) / d(n - 1, 1);
    d(n, 1) = a(n, 1) - sum((u(1 : n - 2, 1).^2) .* d(1 : n - 2, 1)) - d(n - 1, 1) * l(n - 1, 1)^2;
end