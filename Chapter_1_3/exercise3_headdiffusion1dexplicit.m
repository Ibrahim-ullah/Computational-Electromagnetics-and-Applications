% Define parameters for diffusion equation and the range in space and time

L= 1.;      %length of wire
T= 1.;      %final time


%parameters needed to solve the equation within explicit method

maxk = 2500;    %number of time steps
dt = T/maxk;
n=50;       %number of space steps
nint=50;    %the wavefront:intermediate point from which u=0(nint<n)!!

dx = L/n;
cond = 1/4;  %conductivity
b=2.*cond*dt/(dx*dx)    %stability parameter (b=<1)
%the stability condition depends more the on solution propagating from the point
%(n,j), hence we use b = 2*r for the stability parameter


%initial temperature of the wire a sinus

for i = 1:(n+1) 
        x(i) = (i-1)*dx;
        u(i,1)  = sin(pi*x(i));
end

%Temperature at the boundary (T=0)

for k=1:maxk+1
    u(1,k) = 0.;
    u(n+1,k) = 0.;
    time(k) = (k-1)*dt;
end

%Implementation of the explicit method


for i=2:n   
    for k=1:maxk
        u(i,k+1)=u(i,k) + 0.5*b*(u(i-1,k) +u(i+1,k)-2.*u(i,k));
    end
end



%Graph of temperature at different selected times

u=u';
figure(2)
mesh(x,time,u)
title('Temperature using the explicit method')
xlabel('X')
ylabel('Temperature')


% figure(1)
% plot(x,u(:,1),'-',x,u(:,100),'-',x,u(:,300),'-',x,u(:,600),'-')
% title('temerature within the explicit method')
% xlabel('X')
% ylabel('T')

    