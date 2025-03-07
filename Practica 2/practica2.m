function [] = practica2(Nexp, n, ind, norma)

% Nexp: es el numero de sistemas lineales que se van a resolver con la misma matriz (numero de experimentos)
% n: dimension de la matriz

% ind: un indice para elegir el tipo de matriz con el que se trabaja. Valor
% de 1 a 5 para los 5 tipos de matrices consideradas en la práctica

% norma: norma en la que se quiere trabajar.

% disp('Problema para ilustrar el papel del numero de condicion en las')
% disp('cotas para la diferencia entre las soluciones de dos sistemas')
% disp('lineales con la misma matriz y terminos independientes perturbados.')

% disp('Se incluyen, entre otras, las matrices de Hilbert, tridiagonales,')
% disp('de Vandermonde, de Householder, ...')
% disp(' ')
% disp('Para seguir pulse cualquier tecla')
% pause

format short e

% Buscar si hay funciones preprogramadas en Matlab para alguno de los casos

% Se valorará el menor número de sentencias empleadas en toda la
% implementación.

switch ind 
    case 1
        A = hilb(n);
        texto = ' de Hilbert';
    case 2
        A = ((n + 1)^2) * (diag(ones(n - 1, 1), -1) + diag(-2 * ones(n, 1)) + diag(ones(n - 1, 1), 1));
        texto = ' tridiagonal';
    case 3
        A = vander((1 / n) * [1 : n]);
        texto=' de Vandermonde';
    case 4
        v = rand(n, 1);
        A = eye(n) - 2 * (v * v') / (v' * v);
        texto = ' de Householder';
    case 5
        A = triu(-ones(n), 1) + eye(n);
        texto = ' triangular';
end

mu = cond(A, norma); %numero de condicion de la matriz en la norma indicada como dato de entrada

for ii = 1 : Nexp
    % Completar lo que se pide en la práctica
        % 1. Generación aleatoria de la solución 
        % 2. Cálculo del término independiente exacto para esa solución
        % 3. Perturbación del término independiente
        % 4. Cáclulo de la solución del sistema perturbado y cociente buscado
        
        x = rand(n, 1);
        b = A * x;
        delta_b = 1.0e-3 * rand(n, 1) ;
        delta_x = (A \ (delta_b + b));
        delta_x = delta_x - x;

        C(ii,1)= (norm(delta_x, norma) / norm(x, norma)) / (norm(delta_b, norma) / norm(b, norma));
end

L = log(mu);
N = histcounts(C, [exp(L * [-1 : 1 / 10 : 1])])

% pause

figure(ind)
clf
semilogx(exp(L * [-1 : 1 / 10 : 1]), zeros(21, 1), '+', 'LineWidth', 2)
hold on

for ii = 1 : length(N)
    semilogx([exp(L * [-1 + 1 / 20 + (ii - 1) / 10]) exp(L * [-1 + 1 / 20 + (ii - 1) / 10])], [0 N(ii)], 'r-', 'LineWidth', 4)
    if N(ii) > 0
        text(exp(L * [-1 + 1 / 40 + (ii - 1) / 10]), N(ii) + floor(Nexp / 10), num2str(N(ii)))
    end
end

axis([1 / 2 / mu 2 * mu 0 Nexp])
title(['matriz', texto, ', 1/\mu(A) = ', num2str(1 / mu), ', \mu(A) = ', num2str(mu)])
xlabel('Intervalos logaritmicamente equiespaciados en [1/\mu(A), \mu(A)]')
ylabel('Numero de experimentos con cociente en cada intervalo')

