function [] = compardensavscasitri

    % Algoritmo compacto de factorizacion LDL^T 
    % para matrices casitridiagonales simetricas definidas positivas

    % Compara resolverlo con la fórmula tal cual o usando una implementación
    % adaptada a que la matriz sea casitridiagonal.

    disp('Programa que analiza distintas implementaciones') 
    disp('de LDL^t para una matriz casitridiagonal')
    disp(' ')
    disp('Para seguir pulse cualquier tecla')
    pause

    nexp = 6; % numero de veces que se duplica la dimension
    %n=8 % Caso de prueba
    n = 64;

    for j = 1 : nexp
        dimen(j) = n;
        n
    
        % a,b,c vectores COLUMNA
        a =  [2 : n + 1]; %diagonal principal (completar)
        b =  -ones(n - 1, 1); %subdiagonal y superdiagonal (completar)
        c =  -ones(n - 2, 1); % fila inferior (completar)
    
        A = diag(a) + diag(b, -1) + diag(b, 1); % Incluir la diagonal principal a y la
         % subdiagonal y supradiagonal principal b
         % con la función diag

        % Se añade la última fila y última columna
        A(n, 1 : n - 2) = c;
        A(1 : n - 2, n) = c';
    
        % Descomentar/comentar para ver que funciona
        % Dejarlo comentado una vez acabado 
        A
        pause

        % La factorización LDLt con la fórmula tal cual
        tic
        [L D tiempos memoria] = facLDLtdensa(A);
        tiempo1(j) = toc;
       
        % La factorización LDLt afinada a una casitridiagonal
        tic
        [d, l, u] = facLDLtcasitri(a,b,c);
        tiempo2(j) = toc;
    
        % Se duplica la dimensión y se vuelve a analizar
        n = n * 2;
    end

% Representación gráfica (no tocar)
    loglog(dimen,tiempo1,'ro-',dimen,tiempo2,'b*:')
    legend('mi LDL^T (densa)','mi LDL^T (afinada a tridiagonal)','Location','NorthWest')
    title(['Comparativa de tiempos'])
    xlabel('Dimension de la matriz')
    ylabel('Tiempo de CPU')
    
end

function [L D tiempos memoria]=facLDLtdensa(mat)

% Incluir el programa facLDLtdensa construido
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

function [d, l, u] = facLDLtcasitri(a,b,c)
    % Programa que calcula la factorización LDLt para una matriz
    % casitridiagonal simétrica donde:
    %   El vector a contiene la diagonal principal
    %   El vector b contiene la subdiagonal principal (que también es la supradiagonal principal de la matriz)
    %   El vector c contiene la última fila (que también es la última columna de la matriz)
    n = numel(a); %dimension de la matriz

    % Primera fila
    d(1) = a(1);

    % Filas 2 a n-1
    for j = 2 : n - 1
        l(j - 1) = b(j - 1) / d(j - 1);
        d(j) = a(j) - d(j - 1) * l(j - 1)^2 ;
    end

    % Última fila, primera columna
    u(1)= c(1) / d(1);

    % Última fila, columnas 2 a n-2
    for j = 2 : n - 2
        u(j) = (c(j) - u(j) * l(j - 1) * d(j - 1)) / d(j);
    end

    % Última fila, últimos dos elementos
    l(n - 1) = (b(n - 1) - u(n - 2) * l(n - 2) * d(n - 2)) / d(n - 1);
    d(n) = a(n) - sum((u(1 : n - 2).^2) .* d(1 : n - 2)) - d(n - 1) * l(n - 1)^2;
end