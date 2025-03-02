function [] = practica1
    % Practica 1 a completar
    format short e

    % Definimos el intervalo y el numero inicial de nodos N
    a = -2;
    b = 2;
    N = 20;

    % Fijamos el número de veces que se van a duplicar el número de nodos [20,40,80,...,1280]
    numrep = 7;

    % Inicializamos donde almacenar las normas de los errores y de las
    % matrices.
    errmatvserrfor = zeros(numrep, 1);

    error1 = zeros(numrep, 1);
    error2 = zeros(numrep, 1);
    errorinf = zeros(numrep, 1);

    norma1 = zeros(numrep, 1);
    norma2 = zeros(numrep, 1);
    normainf = zeros(numrep, 1);
    normafro = zeros(numrep, 1);

    % Para cada cantidad de nodos, calcular y almacenar en dos tablas:
    % Las normas vectoriales de los errores cometidos
    % Las normas matriciales de las matrices de derivación
    for i = 1 : numrep
         % Se almacena el número de nodos (dimensión del vector)
        dimen(i, 1) = N * ( 2^(i - 1) );

        % Se aproxima la derivada segunda mediante la fórmula y un método matricial
        [x, fdotdot, DD] = derivadasegunda(@mifuncion, a, b, dimen(i, 1));

        % Se aproxima la derivada segunda mediante la fórmula punto a punto
        [xcompar, fdotdotbis] = derivadasegundabis(@mifuncion, a, b, dimen(i, 1));

        % Diferencia entre el método matricial y el método punto a punto
        errmatvserrfor(i, 1) = norm(fdotdotbis - fdotdot);

        % Valor exacto de la derivada segunda
        fsegexacta = dersegex(x);

        % Calcular las normas vectoriales del error
        error1(i)   = norm(fdotdot - fsegexacta, 1);       % Norma 1
        error2(i)   = norm(fdotdot - fsegexacta, 2);      % Norma 2
        errorinf(i) = norm(fdotdot - fsegexacta, "inf");  % Norma Infinito

        % Calcular las normas matriciales de la matriz de derivación
        norma1(i)   = norm(DD, 1);           % Norma 1
        norma2(i)   = norm(DD, 2);          % Norma 2
        normainf(i) = norm(DD, "inf");      % Norma Infinito
        normafro(i) = norm(DD, "fro");      % Norma Frobenius

    end

        %%%%%%%%%
        %%% Sacar por pantalla Error Método Matricial vs Error Método fórmula elemento a elemento
        disp('Error Método Matricial vs Error Método fórmula elemento a elemento');
        tabla_dim_err = [dimen, errmatvserrfor];
        disp(tabla_dim_err);

        % Descomentar el pause para parar y analizar el error 
        pause
        
        % Pulsar una tecla para continuar

        %%% Figura 1 %%%
        % Tabla con la primera columna dimen y los errores al aproximar f''
        disp('Tabla con el valor de N y los errores en las normas 1, 2 e infinito respectivamente');
        tabla_dim_err12inf = [dimen error1, error2, errorinf];
        disp(tabla_dim_err12inf);

        % Tabla que muestre el decrecimiento de los errores en cada norma al duplicar N (orden)
        disp('Tabla con el decrecimiento del error al duplicar N');
        tabla_crec_err = [error1(1 : numrep - 1) ./ error1(2 : numrep), ...
                          error2(1 : numrep - 1) ./ error2(2 : numrep), ...
                          errorinf(1 : numrep - 1) ./ errorinf(2 : numrep)];
        disp(tabla_crec_err);
        % Descomentar el pause para parar y analizar el error 
        pause

        figure(1)
        clf
        loglog(dimen, error1 ./ error2, 'ro-', dimen, error2 ./ errorinf, 'b*-', dimen, error1 ./ errorinf, 'kv-')

        % Completar las etiquetas de los ejes y la leyenda adecuadamente
        xlabel('N');
        ylabel('Errores');
        legend('Error norma 1 / norma 2', 'Error norma 2 / norma inf', 'Error norma 1 / norma inf');

        disp('>----------------------------------------<o>----------------------------------------<')

    %%% Figura 2 %%%
        % Tabla con la primea columna dimen y las normas de las matrices de derivacion
        disp('Norma 1, 2, infinito y de frobenius de la matriz DD respectivamente')
        tabla_norm_mat = [dimen norma1 norma2 normainf normafro];
        disp(tabla_norm_mat);

        % Tabla que muestre el crecimiento de dichas normas con la dimension
        disp('Tabla con el crecimiento de las normas con la dimensión');
        tabla_crec_norm = [norma1(2 : numrep) ./ norma1(1 : numrep - 1), ...
                           norma2(2 : numrep) ./ norma2(1 : numrep - 1), ...
                           normainf(2 : numrep) ./ normainf(1 : numrep - 1), ...
                           normafro(2 : numrep) ./ normafro(1 : numrep - 1)];
        disp(tabla_crec_norm);

        % Descomentar el pause para parar y analizar el comportamiento
        pause

        figure(2)
        clf
        loglog(dimen, norma1, 'ro-', dimen, norma2, 'b*-', dimen, normainf, 'kv-', dimen, normafro, 'gd-');

        % Completar las etiquetas de los ejes y leyenda adecuadamente
        xlabel('N');
        ylabel('Valores de las normas matriciales');
        legend('Norma 1', 'Norma 2', 'Norma infinito', 'Norma de Frobenius');  
end

% Programas auxiliares 

function [x, fdotdot, DD] = derivadasegunda(fun, a, b, N)

    % Programa que calcula la aproximación central a la derivada segunda
    % matricialmente

    % Se trabaja con una malla equiespaciada del intervalo [a, b].
    % x(j) = a + jh, h = (b - a) / N, j = 1, ..., N - 1.
    
    % Como datos de salida da:
    % x(j), j=1, ..., N - 1
    % DD la matriz de derivación de forma que al multiplicar 
    %       fdotdot= DD*fun(x)
    % fdotdot = la aproximación a la derivada segunda para x(j)=a+jh, j=(b-a)/N, j=1,...,N-1.
   
    h = (b - a) / N;
    x = a + h .* [1 : N - 1]';

    DD = (1 / (h^2)) * (diag(ones(N - 2, 1), -1) + diag(-2 * ones(N - 1, 1)) + diag(ones(N - 2, 1), 1));
    fdotdot = DD * fun(x);
end

function [xcompar, fdotdotbis] = derivadasegundabis(fun, a, b, N)

    % Programa que calcula la aproximación central a la derivada segunda
    % con la fórmula elemento a elemento

    % Se asume que es en una malla equiespaciada del intervalo [a,b].

    % La derivada se calcula con la fórmula de orden 2 usando los nodos
    % x(j)=a + j * h, h = (b - a) / N, j=0, ..., N.

    % Como datos de salida da:
    % xcompar=x(j), j=1,...,N-1
    % fdotdotbis = la aproximación a la derivada segunda con la fórmula

    h = (b - a) / N;
    x = a + h * [0 : N]';
    
    fdotdotbis = zeros(N - 1, 1);

    for j = 2 : N
        fdotdotbis(j - 1, 1) = fun(x(j - 1)) - 2 * fun(x(j)) + fun(x(j + 1));
    end
    
    xcompar = x(2 : N - 1, 1);
    fdotdotbis = (1 / h^2) * fdotdotbis;
end

% Programa auxiliar con la función a analizar
function f = mifuncion(x)
    % Poner la función que se desee analizar (Ej.: (x^2-4)(x^2+1))
    f = (x.^2 - 4) .* (x.^2 + 1);
    % Comentar la anterior y poner otra según se solicite 
    %f = exp(- x.^2);
end

% Programa auxiliar con la derivada segunda exacta de la función a analizar
function f = dersegex(x)
    % Poner la derivada segunda exacta de la función que se desee analizar (Ej.: ((x^2-4)(x^2+1))'' )
    f = 12 * (x.^2) - 6;
    % Comentar la anterior y poner otra función según se solicite
    %f = (4 * x.^2 - 2) .* exp(-x.^2);
end 

