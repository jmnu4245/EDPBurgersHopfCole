function [x, incr, I, iter] = KLAM5S(U,k,h,E,tol, maxiter)

lambda=-5;
phi=0;
% U es tan solo dos columnas, contiene todas las filas
iter=0;
incr=tol+1;
I=[];
% Se toma de estimación inicial la columna de un instante anterior
x0=U(2:end-1,1);
% Se añade a la estimación inicial los valores conocidos de la cc1 y cc2
Ux0=[U(1,2); x0; U(end,2)];
[fx,dfx] = funcionpr([U(:,1) Ux0],k,h,E);

while incr>tol && iter<maxiter
    % Se resuelve el sistema: [F'(x^k)]^-1 * F(x^k)
    z=dfx\fx;
    % y^k = x^k - [F'(x^k)]^-1 * F(x^k)
    y=x0-z;
    % Se añade a la estimación y los valores conocidos de la cc1 y cc2
    Uy=[U(1,2); y; U(end,2)];
    [fy,~] = funcionpr([U(:,1) Uy],k,h,E);
    % Cálculo de valores
    vk=(fy'*fy)/(fx'*fx);
    Kk=1/(1+lambda*vk);
    pk=Kk*(1+phi*vk);
    qk=2*Kk*vk;
    % A = pk F(y^k) + qk F(x^k)
    A=pk*fy+qk*fx;
    % Se resuelve el sistema: [F'(x^k)]^-1 * (pk F(y^k) + qk F(x^k))
    z=dfx\A;
    % x^(k+1) = y^k - [F'(x^k)]^-1 * (pk F(y^k) + qk F(x^k))
    x=y-z;
    % Se añaden los valores conocidos de la cc1 y cc2
    Ux=[U(1,2); x; U(end,2)];
    [fx,dfx] = funcionpr([U(:,1) Ux],k,h,E);
    % Criterio de parada
    incr=norm(x-x0, inf) + norm(fx, inf);
    I=[I, incr];
    iter=iter+1;
    % Actualizar x0
    x0=x;
    
end
end