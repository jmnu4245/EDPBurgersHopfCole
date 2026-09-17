function [x, t, U, viter, mI,maxCond] = CrankNicolson_trabajo(cc1,cc2,ci,a,b,nx,T,nt,E, tol, maxiter)
h=(b-a)/nx; x=a:h:b;
k=T/nt; t=0:k:T;
%matriz de nx+1 filas y nt+1 columnas
U=zeros(nx+1,nt+1);
%fijada la t para ci
cix=feval(ci,x);
U(:,1)=cix;
cc1t=feval(cc1,t);
U(1,:)=cc1t;
cc2t=feval(cc2,t);
U(nx+1,:)=cc2t;
C=2:nx;
viter=zeros(1,nt);
mI = zeros(maxiter,nt);

for j=1:nt
    % Newton devuelve nx-1 filas 
    %[vx, incr, I, iter] = Newton_trabajo(U(:,j:j+1),k,h,E,tol,maxiter);
    %[vx, incr, I, iter] = Kurchatov(U(:,j:j+1),k,h,E,tol, maxiter);
    %[vx, incr, I, iter] = TJF4S(U(:,j:j+1),k,h,E,tol, maxiter);
    %[vx, incr, I, iter] = KLAM5(U(:,j:j+1),k,h,E,tol, maxiter);
    %v2 implementa diferencias divididas analíticas
    %[vx, incr, I, iter] = M4_D_v2(U(:,j:j+1),k,h,E,tol, maxiter);
    [vx, incr, I, iter] = M4_D(U(:,j:j+1),k,h,E,tol, maxiter);
    % Actualización de valores
    U(C,j+1)=vx;
    viter(j)=iter;
    mI(1:length(I), j) = I(:);
end
end