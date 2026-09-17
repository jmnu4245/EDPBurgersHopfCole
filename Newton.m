function [x, incr, I, iter] = Newton_trabajo(U,k,h,E,tol, maxiter)

% U es tan solo dos columnas, contiene todas las filas
iter=0;
incr=tol+1;
I=[];
% Se toma de estimación inicial la columna de un instante anterior
x0=U(2:end-1,1);

while incr>tol && iter<maxiter
    % Se añade a la estimación inicial los valores conocidos de la cc1 y cc2
    Ux0=[U(1,2); x0; U(end,2)];
    [fx,dfx] = funcionpr([U(:,1) Ux0],k,h,E);
    % Se resuelve el sistema
    z=dfx\fx;
    x=x0-z;
    % Vector de diferencias entre iteraciones
    I=abs(x-x0);
    % Valor máximo de ese vector de diferencias
    incr=max(I);
    % Criterio de parada: el del paper?
%     I=[I, incr];
    iter=iter+1;
    % Actualizar x0
    x0=x;
end

%ACOC=3;
% if incr>tol
%     disp('Se necesitan más iteraciones')
% else
%     ACOC=log(I(3:end)./I(2:end-1))./log(I(2:end-1)./I(1:end-2));
%     
%     x=vpa(x,6);
%     incr=vpa(incr, 6);
%     ACOC=vpa(ACOC, 6);
% end

end