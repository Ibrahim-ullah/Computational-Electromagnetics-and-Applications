clf;
clear;
clc;

%Physical domain
a=50;
b=50;

%Constants
eps0=8.85e-12;       %Permittivity of free space
mu0=4*pi*1e-7;       %Permeability of free space

%Source parameters
c=3e8;               %Speed of light
freq=3e7;            %Frequency of the source
lambda=c/freq;       %Wavelength of the source
k0=(2*pi)/lambda;    %Wave number
w=2*pi*freq;

%Grid parameters
dx=lambda/20;        %Step size along x-direction
dy=lambda/20;        %Step size along y-direction

x=0:dx:a;
nx=length(x);

y=0:dy:b;
ny=length(y);

%Time parameters
sim_time=120;
R=0.5;               %CFL stability criterion
dt=dx*R/c; 
tsteps=sim_time;

%PML parameters
pmlp=2;              %profile of conductivity inside the PML
refl_th=-70 ;        %Theoretical reflection in dB

sig0x=0.01;
sig0y=0.01;

d_xpml=8;            %pml of 8 cells thick in the x-direction
d_ypml=8;            %pml of 8 cells thick in the y-direction

leftrightpml=[1 1];
lowuppml=[1 1];

%Initializing parameters of pml
sigx=zeros(nx,ny);
sigy=zeros(nx,ny);
mur=ones(nx,ny);
epsr=ones(nx,ny);
mu=zeros(nx,ny);
eps=zeros(nx,ny);
z=zeros(nx,ny);
oneM=ones(nx,ny);
dtM=dt*ones(nx,ny);

%Filling up conductivity values
for i=1:nx
    for j=1:ny
        if leftrightpml(1)>0 && (i<=d_xpml)
            sigx(i,j)=sig0x*(((d_xpml-i)*dx)/(d_xpml*dx))^pmlp;
        end
        
        if leftrightpml(2)>0 && (i>=(nx-d_xpml))
            sigx(i,j)=sig0x*((i-(nx-d_xpml))/(d_xpml*dx))^pmlp;
        end
        
        if lowuppml(1)>0 && (j<=d_ypml)
           sigy(i,j)=sig0y*(((d_ypml-j)*dy)/(d_ypml*dy))^pmlp;
        end
        
        if lowuppml(2)>0 && (j>=(ny-d_ypml))
            sigx(i,j)=sig0y*((j-(ny-d_ypml))/(d_ypml*dy))^pmlp;
        end  
    mu(i,j)=mur(i,j)*mu0;
    eps(i,j)=epsr(i,j)*eps0;
    z(i,j)=sqrt(mu(i,j)/eps(i,j));    
    end
end

%Initialize field matrices
Hx=zeros(nx,ny);
Hy=zeros(nx,ny);
Ez=zeros(nx,ny);
Kx=zeros(nx,ny);
Ky=zeros(nx,ny);

%Initialize curl matrices
CHx=zeros(nx,ny);
CHy=zeros(nx,ny);
CEz=zeros(nx,ny);

%Initialize integration matrices
ICHx=zeros(nx,ny);
ICHy=zeros(nx,ny);
ICEz=zeros(nx,ny);
ICKx=zeros(nx,ny);
ICKy=zeros(nx,ny);

%FDTD loop for the E-field

for t=1:tsteps
    
    %Boundary condition PMC
    Hx(:,1)=0;
    Hx(:,ny)=0;
    Hy(1,:)=0;
    Hy(nx,:)=0;
    
    %Boundary conditions PEC
    Ez(:,1)=0;
    Ez(:,ny)=0;
    Ez(1,:)=0;
    Ez(nx,:)=0;
     
    %Defining a scatterer
    
%     for i=73:74
%         for j=48:52
%             Ez(i,j)=0;
%         end
%     end
    
    %Update Equations
 
    for i=1:nx
        for j=1:ny-1
            Hx(i,j)=Hx(i,j)-((dt/mu(i,j))*((Ez(i,j+1)-Ez(i,j))/dy))+...
                    ((dt*sigx(i,j)/mu(i,j))*Kx(i,j))+((dt*sigy(i,j)/mu(i,j))*Hx(i,j));
            Kx(i,j)=Kx(i,j)-(dt*(Ez(i,j+1)-Ez(i,j)))/dy;
        end
    end
    
    for i=1:nx-1
        for j=1:ny
            Hy(i,j)=Hy(i,j)+((dt/mu(i,j))*((Ez(i+1,j)-Ez(i,j))/dx))+...
                    ((dt*sigy(i,j)/mu(i,j))*Ky(i,j))-((dt*sigx(i,j)/mu(i,j))*Hy(i,j));
            Ky(i,j)=Ky(i,j)+(dt*(Ez(i+1,j)-Ez(i,j)))/dx;
        end
    end
    
    for i=2:nx
        for j=2:ny
            Ez(i,j)=Ez(i,j)+(dt/eps(i,j))*((Hy(i,j)-Hy(i-1,j))/dx)-...
                    (dt/eps(i,j))*((Hx(i,j)-Hx(i,j-1))/dy)-((dt*sigx(i,j)/eps(i,j))*Ez(i,j))-...
                    ((dt*sigy(i,j)/eps(i,j))*Ez(i,j));
        end
    end
    
    Ez(1,1)=Ez(1,1)+(dt/eps(1,1))*((Hy(1,1)+Hy(1,1))/dx)-...
                    (dt/eps(i,j))*((Hx(1,1)+Hx(1,1))/dy)-sigx(1,1)*Ez(1,1)/eps(1,1)-sigy(1,1)*Ez(1,1)/eps(1,1);
    
    %Defining the Source 
    source=sin(2*pi*freq*t*dt);
    Ez(floor(nx/2),floor(ny/2))=2*source;
    
    %Plotting the Ez wave
    [xx,yy]=meshgrid(x,y);
    pcolor(x,y,Ez');
    rectangle('position',[36 23.5 2 2],'facecolor','m');
    shading interp
    xlabel('X \rightarrow');
    ylabel('\leftarrow Y');
    zlabel('E_z \rightarrow');
    titlestring=['\fontsize{14} FDTD Setup with UPML at timestep(',num2str(t),')'];
    title(titlestring,'color','k');
    
    axis([0 a 0 b -1 1]);
    view(0,90);
    caxis([-1,1]);
    colormap('jet');
    colorbar;
    getframe;
    
end