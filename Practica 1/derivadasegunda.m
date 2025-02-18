% Programa para calcular la derivada segunda usando
% diferencias divididas de segundo orden

function [x, fdotdot, DD] = derivadasegunda(fun, a, b, N)
    h = (b - a) / N;
    x = a + h .* [1 : N - 1]';

    DD = (1 / (h^2)) * (diag(ones(N - 2, 1), -1) + diag(-2 * ones(N - 1, 1)) + diag(ones(N - 2, 1), 1));
    fdotdot = DD * fun(x);
end