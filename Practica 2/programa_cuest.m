n = [4 : 4 : 40];
n_exp = 1000;

CC = zeros(n_exp, 6);

fin_vals = zeros(6, 10);
    
for jj = 1 : 10
    A = triu(-ones(n(jj)), 1) + eye(n(jj));

     % Realizamos los experimentos
    for ii = 1 : n(jj)
        x = rand(n(jj), 1);
        b = A * x;
        delta_b = 1.0e-3 * rand(n(jj), 1) ;
        delta_x = (A \ (delta_b + b));
        delta_x = delta_x - x;

        CC(ii, 1) = norm(delta_x, Inf);
        CC(ii, 2) = norm(delta_b, Inf);
        CC(ii, 3) = CC(ii, 2) / norm(b, Inf);
        CC(ii, 4) = norm(A, Inf);
        CC(ii, 5) = norm(inv(A), Inf);
        CC(ii, 6) = cond(A, Inf);
    end
    
    fin_vals(1, jj) = sum(CC(:, 1)) / 1000; % norma x
    fin_vals(2, jj) = sum(CC(:, 2)) / 1000; % norma db
    fin_vals(3, jj) = sum(CC(:, 3)) / 1000; % db / b
    fin_vals(4, jj) = sum(CC(:, 4)) / 1000; % Norma mat
    fin_vals(5, jj) = sum(CC(:, 5)) / 1000; % Norma inv
    fin_vals(6, jj) = sum(CC(:, 6)) / 1000; % Num de condicion
end

loglog(n, fin_vals(1, :), n, fin_vals(2, :), n, fin_vals(3, :), n, fin_vals(4, :), n, fin_vals(5, :), n, fin_vals(6, :));
l = legend({'$\vert\vert\mathbf{x}\vert\vert_\infty$', '$\vert\vert\delta\mathbf{b}\vert\vert_\infty$', ...
    '$\vert\vert\delta\mathbf{b}\vert\vert_\infty / \vert\vert\mathbf{b}\vert\vert_\infty$', ...
    '$\vert\vert A_n\vert\vert_\infty$', ...
    '$\vert\vert A_n^{-1}\vert\vert_\infty$', ...
    '$\mu(A_n)$'}, 'Location', 'northwest', 'NumColumns', 2);
set(l, 'Interpreter', 'latex')
