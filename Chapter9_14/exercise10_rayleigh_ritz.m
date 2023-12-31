%Rayleigh Ritz method
clear;
clc;

%Declaration of variables
syms x;
syms y;
syms N;
N=3;                %Order of Approximation
rho=-1;

%Initialization of derivative matrices
u=sym('u',[N,N]);
ux=sym('ux',[N,N]);         %First derivative matrix
uxx=sym('uxx',[N,N]);       %Second derivative matrix
uy=sym('uy',[N,N]);
uyy=sym('uyy',[N,N]);

%Initialization of intermediate matrices
a=sym('a',[N,N]);
b=sym('b',[N,N]);
c=sym('c',[N,N]);
d=sym('d',[N,N]);
e=sym('e',[N,N]);
f=sym('f',[N,N]);
psi=sym('psi',[N,N]);       %Coefficients matrix
phi=0;                      %Solution matrix

%Formation of matrices
for m=0:N-1
    for i=0:N-1
        u(i+1,m+1)=(1-x^2)*(1-y^2)*(x^(2*m)*y^(2*i)+...
                    x^(2*i)*y^(2*m));
    end
end

%Calculation of derivative of matrices
for m=0:N-1
    for i=0:N-1
        ux(i+1,m+1)=diff(u(i+1,m+1),'x');
        uy(i+1,m+1)=diff(u(i+1,m+1),'y');
        uxx(i+1,m+1)=diff(ux(i+1,m+1),'x');
        uyy(i+1,m+1)=diff(uy(i+1,m+1),'y');
    end
end

%Calculation of inner products on LHS
for i=1:N
    for j=1:N
        a(i,j)=(uxx(i,j)+uyy(i,j))*u(i,j);
        c(i,j)=int(a(i,j),x,-1,1);
        d(i,j)=int(c(i,j),y,-1,1);
    end
end

%Calculation of inner products on RHS
for i=1:N
    for j=1:N
        b(i,j)=rho*u(i,j);
        e(i,j)=int(b(i,j),x,-1,1);
        f(i,j)=int(e(i,j),y,-1,1);
    end
end

%Matrix inversion
for i=1:N
    for j=1:N
       psi(i,j)=d(i,j)\f(i,j);
    end
end

%Calculation for solution
for i=1:N
    for j=1:N
       phi=phi+psi(i,j)*u(i,j);
    end
end