close all;
clear all;
clc;
%Define Parameters for simulation

%Computational Domain in x-direction
c=3e8;       %Speed of light 
freq=3e7;    %Frequency in Hz
eps_r=1;     %Relative dielectric constant of medium
lambda=(c/freq)/sqrt(eps_r); %Wavelength
xdim=100;
dx=lambda/10;  %x-position step
x=0:dx:xdim;
xsteps=length(x);

%Total simulation time
time_tot=400;
Rx=0.5;
dt=Rx*dx/c; %Courant Stability condition
%t=0:dt:time_tot;
tsteps=time_tot

%Grid dimension in y-direction
ydim=100;
Ry=0.5;
dy=c*dt/Ry;
y=0:dy:ydim;
ysteps=length(y);

%position of source
xsource=floor(0.5*xsteps);
ysource=floor(0.5*ysteps);

%Initialize field vectors
Ez=zeros(ysteps,xsteps);
Hx=zeros(ysteps,xsteps);
Hy=zeros(ysteps,xsteps);

Ex2=zeros(tsteps,xsteps);
Exlast_1=zeros(tsteps,xsteps);
Ey2=zeros(tsteps,ysteps);
Eylast_1=zeros(tsteps,xsteps);

for n=1+ceil(1/min(Rx,Ry)):tsteps
    for l=1:xsteps
        for m=1:ysteps-1
           Hx(m,l)=Hx(m,l)-Ry*(Ez(m+1,l)-Ez(m,l));
        end
    end
    
     for m1=1:ysteps
         for l1=1:xsteps-1
             Hy(m1,l1)=Hy(m1,l1)-Rx*(Ez(m1,l1+1)-Ez(m1,l1));
         end
     end
     
     for m2=2:ysteps
         for l2=2:ysteps
             Ez(m2,l2)=Ez(m2,l2)-Rx*(Hy(m2,l2)-Hy(m2,l2-1))-Ry*(Hx(m2,l2)-Hx(m2-1,l2));
         end
     end
     
     %Absorbing boundary conditions using one-way wave equation
     
     %In x-direction
     Ex2(n,:)=Ez(:,2);
     Ez(:,1)=Ex2(n-1/Rx,:);
     Exlast_1(n,:)=Ez(:,xsteps-1);
     Ez(:,xsteps)=Exlast_1(n-1/Rx,:);
   
     %In y-direction
     Ey2(n,:)=Ez(2,:);
     Ez(1,:)=Ey2(n-1/Ry,:);
     Eylast_1(n,:)=Ez(ysteps-1,:);
     Ez(ysteps,:)=Eylast_1(n-1/Ry,:);

     %Defining source
     source=sin((2*pi*freq*dt));
     
     %Assigning source
     Ez(ysource,xsource)=source;
     
     %Plotting Ez wave
     surf(x,y,Ez,'Linewidth',2);
     xlabel('X \rightarrow');
     ylabel('\leftarrow Y');
     zlabel('E_z \rightarrow');
     titlestring=['\fontsize{20} Plot of Ez vs X&Y for 2D FDTD'];
     title(titlestring,'color','k')
     axis([0 100 0 100 -1 1]);
     colormap('winter')
     pause(0.2)
end