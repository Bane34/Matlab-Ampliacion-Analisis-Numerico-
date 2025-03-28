function [] = practica3(ind_disp)

% Algoritmo compacto de factorizacion LDL^T 
% para matrices casitridiagonales simetricas definidas positivas

% Si ind_disp=1 la matriz, aunque tridiagonal, se considera llena 
% Si ind_disp=2 la matriz se considera dispersa 

disp('Problema para ilustrar la adaptacion del algoritmo') 
disp('compacto de factorizacion LDL^T al caso de matrices')
disp('tridiagonales simetricas definidas positivas') 
disp(' ')
disp('Para seguir pulse cualquier tecla')
pause

switch ind_disp
  case 1
    nexp = 9; % numero de veces que se duplica la dimension
    texto = ' MATLAB no sabe que A es dispersa';
  case 2
    nexp = 12; % numero de veces que se duplica la dimension
    texto = ' MATLAB sabe que A es dispersa';
end

%n=8;    % Caso de prueba
n = 64;

for j = 1 : nexp
    dimen(j) = n;
    n
    % a,b,c vectores COLUMNA
    a = [2 : n + 1]'; %diagonal principal (completar)
    b = -ones(n - 1, 1); %subdiagonal y superdiagonal (completar)
    c = -ones(n - 2, 1); % fila inferior (completar)
    
    switch ind_disp
        case 1
            A = diag(a) + diag(b, -1) + diag(b, 1); % Incluir la diagonal principal a y la
                 % subdiagonal y supradiagonal principal b
                 % con la función diag

            % Se añade la última fila y última columna
            A(n, 1 : n - 2) = c';
            A(1 : n - 2, n) = c';
            
            % Descomentar/comentar para ver que funciona
            % Dejarlo comentado una vez acabado

             %A
             pause
        
        case 2
            A = spdiags(a, 0, n, n);
            A = spdiags(b, -1, A);
            A = spdiags([1 b']', 1, A);
            A(n, 1 : n - 2) = c';
            A(1 : n - 2, n) = c';
        
            % Descomentar/comentar para ver que funciona
            % Dejarlo comentado una vez acabado

            % A
            % full(A)
            % pause
    end

    % La factorización LDLt afinada a una casitridiagonal
    tic
    [d, l, u] = facLDLtcasitri(a,b,c);
    tiempo1(j) = toc;
    
    % La factorización LU de Matlab
    tic
    [L, U, P] = lu(A);
    tiempo2(j) = toc;

    % La factorización de Cholesky de Matlab
    tic
    % G = chol(A, 'lower');
    G = chol(A, 'lower');
    tiempo3(j) = toc;

    % Se duplica la dimensión y se vuelve a analizar
    n = n * 2;
end

% Representación gráfica (no tocar)
    figure(ind_disp)
    loglog(dimen,tiempo1,'ro-',dimen,tiempo2,'b*:',dimen,tiempo3,'kv--')
    legend('mi LDL^T (afinada casitridiagonal)','LU','Cholesky','Location','NorthWest')
    title(['Resultados cuando ', texto])
    xlabel('Dimension de la matriz')
    ylabel('Tiempo de CPU')

end

function [d, l, u] = facLDLtcasitri(a,b,c)

% Incluir el programa facLDLtcasitri construido
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
