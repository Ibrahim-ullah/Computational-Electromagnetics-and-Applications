clear all;
clc;

%Define parameters for the simulation

c = 3e8; %speed of light
freq_in = 3e7; %frequency in Hz
eps_r =1; %relative dielectric constant of medium

lambda = (c/freq_in)/sqrt(eps_r);
xdim=100;
dx=lambda/10; %x-position step
x=0:dx:xdim;
xsteps=length(x);

%total simulation time 
time_tot=400;
Rx=0.5; %courant constant of stability R < 1
dt=Rx*dx/c;
%t = 0:dt:time_tot;
tsteps = time_tot;

%Grid dimension in  y(ydim) direction
ydim = 100;
Ry = 0.5;
dy = c*dt/Ry;
y = 0:dy:ydim;
ysteps = length(y);

%position of source
xsource = floor(xsteps/2);
ysource = floor(ysteps/2);

%Initialize of field vectors
Ez = zeros(ysteps,xsteps);
Hx =  zeros(ysteps,xsteps);
Hy =zeros(ysteps,xsteps);

Ex2 = zeros(tsteps,xsteps);
Exlast_1 = zeros(tsteps,xsteps);
Ey2 = zeros(tsteps,ysteps);
Eylast_1 = zeros(tsteps,ysteps);


for n = 1+ceil(1/min(Rx,Ry)):tsteps
    
for l=1:xsteps
    for m = 1:ysteps-1
        Hx(m,l) = Hx(m,l)-Ry*(Ez(m+1,l)-Ez(m,l));
    end
end



for m1=1:ysteps
    for l1 = 1:xsteps-1
        Hy(m1,l1) = Hy(m1,l1)-Rx*(Ez(m1,l1+1)-Ez(m1,l1));
    end
end


%%absorbing boundary condition using one way wave equation
%%performance deteriorte for off-normal incidence

%in x-direction
    Ex2(n,:) = Ez(:,2);
    Ez(:,1) = Ex2(n-1/Rx,:);
    Exlast_1(n,:) = Ez(:,xsteps-1);
    Ez(:,xsteps) = Exlast_1(n-1/Rx,:);
    
% in y-direction
    Ey2(n,:) = Ez(2,:);
    Ez(1,:) = Ey2(n-1/Rx,:);
    Exlast_1(n,:) = Ez(ysteps-1,:);
    Ez(ysteps,:) = Eylast_1(n-1/Ry,:);
    
%% Define Source

source = sin(((2*pi*(freq_in)*n*dt)));

%for a single pulse gaussian force
% if n<=42
%       pulse = (10-15*cos(n*pi/20)+6*cos(2*n*pi/20)-cos(3*n*20);   
% end
%
%for a periodic gaussian source
%pulse = (10-15*cos(n*pi/20)+6*cos(2*n*pi/20) - cos(3*n*20);
%
%%Assigining source
    Ez(ysource,xsource) = source;   %hard sources acts as Metal
%  Ez(ysource,xsource) = Ez(ysource,xsource)+source; %soft source
%  Ez(ysource+ysource/2,xsource+xsource/2) = souce;
%  Ez(ysource-ysource/2,xsource-xsource/2) = -1*souce;
%  Hy(xsource) = source;
%
%% Plotting Ez-wave
    surf(x,y,Ez,'linewidth',2);
    xlabel('X \rightarrow');
    ylabel('\leftarrow Y');
    zlabel('E_z \rightarrow');
    titlestring=['\fontsize{10}Plot of E_z vs X & Y for 2D FDTD',num2str(n)];
    title(titlestring,'color','k');
    %axis([0 100 0 100 -1 1]);
    colormap(summer)
    shading interp; % Smooth shading for a better visual appearance
    view(30, 45); % Adjust the view angle for better perspective
    pause(0.2)
    
    
end
    