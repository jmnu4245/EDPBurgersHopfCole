function [fx,dfx] = funcionpr(U,k,h,E)

%u es para el trabajo, x es para la prueba
%fx=[x(1)^2+x(2)^3-5;x(1)*sin(x(1))-x(2)*cos(x(1))-1];
%dfx=[2*x(1) 3*x(2)^2;sin(x(1))+x(1)*cos(x(1))+x(2)*sin(x(1)) -cos(x(1))];

%j indice del tiempo
j=1;

%f tiene nx-1 funciones

nx=length(U(:,1))-1;

%un vector que solo almacena las incógnitas
fx=zeros(nx-1,1); %vector columna
A=k*E/(2*h^2);
B=k/(4*h);

for i=2:nx
    fx(i-1,1)=(1+2*A)*U(i,j+1)+(-1+2*A)*U(i,j)-A*(U(i+1,j)+U(i-1,j)+U(i+1,j+1)+U(i-1,j+1))+B*(U(i,j)*(U(i+1,j)-U(i-1,j))+U(i,j+1)*(U(i+1,j+1)-U(i-1,j+1)));
end

dfx=zeros(nx-1,nx-1); %matriz
ds=-A+B*U(2:nx-1,j+1);
dp=1+2*A+B*(U(3:nx+1,j+1)-U(1:nx-1,j+1));
di=-A-B*U(3:nx,j+1);
dfx=diag(dp)+diag(ds,1)+diag(di,-1);


end