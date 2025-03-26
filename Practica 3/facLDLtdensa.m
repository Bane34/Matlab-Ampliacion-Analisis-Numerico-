function [L, D, tiempos, memoria]=facLDLtdensa(mat)
% No tocar las sentencias escritas ni usar las variables que empleo yo.
% No se pueden borrar variables dentro del programa (comando clear o similares).
% No se puede alterar la matriz mat.

% Se asume que mat es una matriz que admite la descomposición LDLt, es decir,
% se sabe que es simétrica.
% Además, se asume que mat es densa. Es decir, todos sus elementos pueden ser no nulos.
    tic
    n = size(mat, 1);
    d = zeros(n, 1);
    l = eye(n, n);

    d(1, 1) = mat(1, 1);
    
    for i = 2 : n
        l(i, 1) = mat(i, 1) / d(1);

        % Calculamos los l_ij
        for j = 2 : i - 1
            suma = 0;
    
            for k = 1 : j - 1
                suma = suma + l(j, k) * d(k) * l(i, k);
            end

            l(i, j) = (mat(i, j) - suma) / (d(j));
        end

        % Calculamos los d_ii
        suma = 0;
        for k = 1 : i - 1
            suma = suma + (l(i, k)^2) * d(k);
        end

        d(i) = mat(i, i) - suma;
        
    end

    L = l;
    D = diag(d);

    tiempos = toc;
    memauxil = whos;
    memoria = sum([memauxil.bytes]) - 8 * (size(mat, 1)^2);
end
