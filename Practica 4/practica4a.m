function [] = practica4a()

    disp('Practica 4a')
    disp('Problema para ilustrar la no conveniencia de formar y resolver')
    disp('las ecuaciones normales cuando se quiere resolver un sistema')
    disp('lineal en el sentido de minimos cuadrados.')
    disp(' ')
    disp('Para seguir pulse cualquier tecla')

    pause

    clear variables

    k0 = 0;

    for i = 1 : 10,
        [A, b, aEN, aQR] = ajuste(k0);
        condA(i,1) =  cond(A, 2)% Número de condicion de A (norma 2).
        condAtA(i,1) = cond(A' * A, 2) % Número de condicion de A^TA.

        % Cálculo de los errores en la solución (norma(2) del vector con los errores en las soluciones).
        normENerr(i, 1) = norm([1; 1; 1] - aEN, 2); % Error con las ecuaciones normales.
        normQRerr(i, 1) = norm([1; 1; 1]- aQR, 2);  % Error con la factorización QR.

        % Cálculo del error en el residuo (norma(2) del vector residuo = b-A*(solución obtenida)).
        normENres(i, 1) =  norm(A * aEN - b, 2); % Residuo con las ecuaciones normales.
        normQRres(i, 1) =  norm(A * aQR -b, 2); % Residuo con la factorización QR.

        kk(i, 1) = k0;
        k0 = k0 +100;
    end 

% Tabla con los resultados
disp(' ')
disp('  k0        mu_2(A)       mu_2(A^TA)       Error          Error         Residuo       Residuo')
disp('                                        Ec. Normales     Fact. QR     Ec. Normales    Fact. QR')
disp(' ')
fprintf(' %3i     %2.4e      %2.4e      %2.4e     %2.4e     %2.4e    %2.4e \n',  [kk';condA';condAtA';normENerr';normQRerr';normENres';normQRres'])

disp(' ')
disp('Para seguir pulse cualquier tecla')
pause

figure(1)
clf
loglog(condA,normENerr,'-*',condA,normQRerr,':o')
title('Error en la solución (- ecuaciones normales, : factorizacion QR)')
xlabel('Numero de condicion de la matriz')
ylabel('Norma 2 del error en la solucion')

end

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