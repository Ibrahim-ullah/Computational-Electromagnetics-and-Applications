% Define parameters for advection eqn & the range in space and time

Lmax= 1.;      %maximum length
Tmax= 1.;      %final time
c =1;           %Advection velocity

%parameters needed to solve the equation within explicit method

maxt = 2500;    %number of time steps
dt = Tmax/maxt;
n=50;       %number of space steps
nint=50;    %the wavefront:intermediate point from which u=0(nint<n)!!

dx = Lmax/n;
%r=c*dt/(2.*dx);  
r=0.3;

% the lax method is stable for r<1/2 but it gets diffused unless r=1/2
% initial value of the function u(amplitude of wave)

for i =1:(n+1)
    if i < nint
        u(i,1) = 1.;
    else u(i,1)=0.;
    end
    x(i) = (i-1)*dx;
end

%value of amplitude at boundary

for k=1:maxt+1
    u(1,k)=1.;
    u(n+1,k)=0;
    time(k) = (k-1)*dt;
end


%Implementation of the lax method
for k=1:maxt %time loop
    for i=2:n  %space loop
        u(i,k+1) = 0.5*((u(i+1,k) + u(i-1,k)) - r*(u(i+1,k) - u(i-1,k)));
    end
end


%Graph of temperature at different selected times

u=u';
figure(2)
mesh(x,time,u)
title('Temperature using the explicit method')
xlabel('X')
ylabel('Temperature')




