function U = EsquemaB_Hopf_Cole(x,t,v,cc1,cc2,ci)

m=length(x);
n=length(t);
U=zeros(m,n);

cix=feval(ci,x);
U(:,1)=cix;
cc1t=feval(cc1,t);
U(1,:)=cc1t;
cc2t=feval(cc2,t);
U(m,:)=cc2t;

h=abs(x(2)-x(1));
k=t(2)-t(1);

if h^2/(4*E*k)<1
    disp('No se cumple la condición')
end


for j=1:n-1
    for i=2:m-1
        Ui_menos=U(i-1,j);
        Ui=U(i,j);
        Ui_mas=U(i+1,j);

        a=Ui;
        b=(Ui_mas-Ui_menos)/(2*h);
        c=(Ui_mas+Ui_menos-2*Ui)/h^2;

        k0=6*v;
        k1=24*b*v;
        k2=(36*b^2+6*a*c)*v;
        k3=a^3*c+(24*b^3+12*a*b*c)*v;
        k4=a^3*b*c+6*(b^4+a*b^2*c)*v;

        s1=-6*a*b*v+12*c*v^2;
        s2=(-18*a*b^2+6*a^2*c)*v+24*b*c*v^2;
        s3=-18*a*b^3*v+12*b^2*c*v^2;
        s4=-a^4*b*c-6*(a*b^4+a^2*b^2*c)*v;

        U(i,j+1)=U(i,j)+(s1*k+s2*k^2+s3*k^3+s4*k^4)/(k0+k1*k+k2*k^2+k3*k^3+k4*k^4);

    end
end

surf(t,x,U)
xlabel('tiempo (t)')
ylabel('espacio (x)')
zlabel('solución (u)')

end