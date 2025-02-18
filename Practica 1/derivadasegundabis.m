function [xcompar, fdotdotbis] = derivadasegundabis(fun, a, b, N)
    h = (b - a) / N;
    x = a + h * [0 : N]';
    
    fdotdotbis = zeros(N - 1, 1);

    for j = 2 : N
        fdotdotbis(j - 1, 1) = fun(x(j - 1)) - 2 * fun(x(j)) + fun(x(j + 1));
    end
    
    xcompar = x(2 : N - 1, 1);
    fdotdotbis = (1 / h^2) * fdotdotbis;
end