% matlab code for one slit diffraction experiment

close all;
clear all;
clc;
set(0,'DefaultFigureRenderer','OpenGL');
fig1 = figure('color','w');
set(fig1,'Name','FDTD analysis with PML');
set(fig1,'NumberTitle','off');

%%constants
e0=8.852e-12;   %permittivity of free space

%%source parameters
c = 3e8;    %speed of EM wave
freq = 3e9; %frequency of EM wve = 3GHz
lambda = c/freq;    %wavelength of EM wave

%%Grid parameters
a=2;    % x length of the box = 2m
b=4;    % y length of the box = 4m
dx = lambda/5;  %mesh size along x-direction
dy = lambda/5;  %mesh size along y-direction

x=0:dx:a;
Nx=length(x);

y=0:dy:b;
Ny=length(y);

%%time parameters
time_tot = 9000;

R = 0.5;
%CFL stability condition <= 0.7 for 2D wave equation

dt=R*dx/c;
tsteps=time_tot;

% compute PML parameters
Nx2=2*Nx;
Ny2=2*Ny;
NPML = [0 20 20 20];

sigx = zeros(Nx2, Ny2);
for i=1:2*NPML(1)
    i1 = 2*NPML(1) - i + 1;
    sigx(i1,:) = (0.5*e0/dt)*(i/2/NPML(1))^3;
end


for i=1:2*NPML(2)
    i1=Nx2 - 2*NPML(2) + i;
    sigx(i1,:) = (0.5*e0/dt)*(i/2/NPML(2))^3;
end


sigy = zeros(Nx2, Ny2);
for j=1:2*NPML(3)
    j1 = 2*NPML(3) - j + 1;
    sigy(:,j1) = (0.5*e0/dt)*(j/2/NPML(3))^3;
end


for i=1:2*NPML(4)
    j1=Ny2 - 2*NPML(4) + j;
    sigy(:,j1) = (0.5*e0/dt)*(j/2/NPML(4))^3;
end

URxx=1;
URyy=1;
ERzz=1;


%% Compute update coefficients

sigHx = sigx(1:2:Nx2,2:2:Ny2);
sigHy = sigy(1:2:Nx2,2:2:Ny2);
mHx0 = (1/dt) + sigHy/(2*e0);
mHx1 = ((1/dt) - sigHy/(2*e0))./mHx0;
mHx2 = -c./URxx./mHx0;
mHx3 = -(c*dt/e0) * sigHx./URxx./mHx0;
sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Nx2,1:2:Ny2);
mHy0 = (1/dt) + sigHx/(2*e0);
mHy1 = ((1/dt) - sigHx/(2*e0))./mHy0;
mHy2 = -c./URyy./mHy0;
mHy3 = -(c*dt/e0) * sigHy./URyy ./mHy0;
sigDx = sigx(1:2:Nx2, 1:2:Ny2);
sigDy = sigy(1:2:Nx2, 1:2:Ny2);
mDz0 = (1/dt) + (sigDx + sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2);
mDz1 = (1/dt) - (sigDx + sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2);
mDz1 = mDz1 ./ mDz0;
mDz2 = c./mDz0;
mDz4 = -(dt/e0^2)*sigDx.*sigDy./mDz0;
mEz1 = 1/ERzz;

%% Initiate field matrices
Dz = zeros(Nx,Ny);
Ez = zeros(Nx,Ny);
Hx = zeros(Nx,Ny);
Hy = zeros(Nx,Ny);

%% Initialize curl matirces

CEx = zeros(Nx,Ny);
CEy = zeros(Nx,Ny);
CHz  = zeros(Nx,Ny);

%% Initialize integration matrices

ICEx = zeros(Nx,Ny);
ICEy = zeros(Nx,Ny);
IDz  = zeros(Nx,Ny)

%%starting from the main FDTD loop
for t=1:9000
    
    % Defining slit
    Ez(5,1:(floor(Ny/2)-2)) = 0;
    Ez(5,(floor(Ny/2)+2):Ny) = 0;
    
    %Calculating CEx
    for i=1:Nx
        for j=1:Ny-1
            CEx(i,j)=(Ez(i,j+1)-Ez(i,j))/dy;
        end
        CEx(i,Ny)=(0-Ez(i,j))/dy;  %for tackling Y-high side. in the boundary we can only get highest value of j
    end
    
    %Calculating CEy
    for j=1:Ny
        for i=1:Nx-1
            CEy(i,j) =-(Ez(i+1,j)-Ez(i,j))/dx;
        end
        CEy(Nx,j) = -(0-Ez(Nx,j))/dx; %for tackling X-High side
    end
    
    %update H integrations
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;
    
    %update H fields
    Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
    Hy = mHy1.*Hy + mHy2.*CEy + mHy3.*ICEy;
    
    %compute CHz
    %curl equations automatically include PEC BC
    for j=2:Ny
        CHz(1,j) = ((Hy(1,j)-0)/dx)-((Hx(1,j)-Hx(1,j-1))/dy);
        for i=2:Nx
            CHz(i,j) = ((Hy(i,j)-Hy(i-1,j))/dx)-((Hx(i,j)-Hx(i,j-1))/dy);
        end
    end
    
    %update d integration
    IDz = IDz + Dz;
    
    %update D field
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;
    
    %for sine wave mode
    source=sin(((2*pi*(freq)*t/2*dt)));
    
    %Assigning source on -low edge
    for i=1:Ny
        Dz(1,j) = sin((pi*(j-1)*dy)/b)*source; %soft source.mode function
    end
    
    %update Ez field
    Ez = mEz1.*Dz;
    Ez(:,1) = 0;
    
    
    %Plotting Ez wave
    [yy , zz] = meshgrid(y,x);
    pcolor(x,y,Ez');
    shading interp;
    xlabel('X \rightarrow');
    ylabel('\leftarrow Y');
    zlabel('E_z \rightarrow');
    titlestring=['\fontsize{20}Slit Experiment at timestep=',num2str(t)];
    title(titlestring,'color','k');
    axis([0 a 0 b -1 1]);
    view(0,90)
    caxis([-1, 1])
    colormap(jet);
    colorbar;
    getframe;
    
end
      
        
    
    
    





