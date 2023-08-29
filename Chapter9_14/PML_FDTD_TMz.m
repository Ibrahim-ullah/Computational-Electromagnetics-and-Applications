%Grid dimension in x (xdim) direction

c=3e8; %speed of light
freq_in=3e7; %freq in Hz
eps_r = 1; %relative dielectric constant of medium

lamda = (c/freq_in)/sqrt(eps_r);
xdim=100;
dx=lambda/10; %x-position step
x=0:dx:xdim;
xsteps = length(x);


%Grid dimension in  y(ydim) direction
ydim = 100;
Ry = 0.5;
dy = c*dt/Ry;
y = 0:dy:ydim;
ysteps = length(y);

%position of source
xsource = floor(xsteps/2);
ysource = floor(ysteps/2);

%Defining PML parameters

length_pml_x = xsteps/5; %width of pml along x-direction
length_pml_y = ysteps/5; %width of pml along y-direction

xn = zeros(1,xsteps);
yn = zeros(1,ysteps);
fi1 = zeros(1,xsteps);
fi2 = ones(1,xsteps);
fi3 = ones(1,xsteps);
gi2 = ones(1,xsteps);
gi3 = ones(1,xsteps);
fj1 = zeros(1,ysteps);
fj2 = ones(1,ysteps);
fj3 =  ones(1,ysteps);
gj2 = ones(1,ysteps);
gj3 = ones(1,ysteps);


gaz = ones(ysteps,xsteps);
gbz = zeros(ysteps,xsteps);

for ii=1:length_pml_x
    xnum = length_pml_x-ii; %making the transition at problem
    xn(ii) = 0.333*((xnum/length_pml_x)^3);
    fi1(ii)=xn(ii); %for 1-side of PML in x-direction
    fi1(xsteps-ii)=xn(ii); %for 2nd-side of PML in x-direction
    fi2(ii) = 1/(1+xn(ii));
    fi2Ixsteps-ii)=1/(1+xn(ii));
    fi3(ii)=(1-xn(ii))/(1+xn(ii)); %for 1-side of PML in x-direction   
    gi2(ii)=1/(1+xn(ii));
    gi2(xsteps-ii)=(1-xn(ii))/(1+xn(ii));
    gi3(ii)=(1-xn(ii))/(1+xn(ii));
    gi3(xsteps-ii)=(1-xn(ii))/(1+xn(ii));
    gaz(:,ii)=1/(1+2*xn(ii));
    gaz(:,xsteps-ii)=1/(1+2*xn(ii));
    gbz(:,ii)=2*xn(ii);
    gbz(:,xsteps-ii)=2*xn(ii);
end


for jj=1:length_pml_y
    ynum = length_pml_y-jj; %making the transition at problem
    yn(jj) = 0.333*((ynum/length_pml_y)^3);
    fj1(jj)=yn(jj); %for 1-side of PML in x-direction
    fj1(ysteps-jj)=yn(jj); %for 2nd-side of PML in x-direction
    fj2(jj) = 1/(1+yn(jj));
    fj2(ysteps-jj)=1/(1+yn(jj));
    fi3(jj)=(1-yn(jj))/(1+yn(jj)); %for 1-side of PML in x-direction   
    gi2(jj)=1/(1+yn(jj));
    gi2(ysteps-jj)=(1-yn(jj))/(1+yn(jj));
    gi3(jj)=(1-yn(jj))/(1+yn(jj));
    gi3(ysteps-jj)=(1-yn(jj))/(1+yn(jj));
    gaz(:,jj)=1/(1+2*yn(jj));
    gaz(:,ysteps-jj)=1/(1+2*yn(jj));
    gbz(:,jj)=2*yn(jj);
    gbz(:,ysteps-jj)=2*yn(jj);
end


%Initialize of field vectors
Dz=zeros(ysteps,xsteps);
Ez=zeros(ysteps,xsteps);
Hx=zeros(ysteps,xsteps);
Hy=zeros(ysteps,xsteps);

IHx=zeros(ysteps,xsteps);
IHy=zeros(ysteps,xsteps);
Iz=zeros(ysteps,xsteps);

for n=1:steps
%%Implementing 2D-FDTD
%Calculating Hx
for l =1:steps
    for m = 1:ysteps-1
        curl_ex=Ez(m,l)-Ez(m+1,l);
        IHx(m,l)=IHx(m,l)+fi1(l)*curl_ex;
        Hx(m,l)=fj3(m)*Hx(m,l)+fj2(m)*Ry*(curl_ex+IHx(m,l));
    end
end


%calculating Hy
for m1 =1:ysteps
    for l1 = 1:xsteps-1
        curl_ey=Ez(m,l1+1)-Ez(m1,l1);
        IHy(m1,l1)=IHy(m1,l1)+fj1(m1)*curl_ey;
        Hy(m1,l1)=fj3(l1)*Hy(m1,l1)+fj2(l1)*Rx*(curl_ey+IHy(m1,l1));
    end
end

%Calculating Dz

for m2=2:ysteps
    for l2=2:xsteps
        Dz(m2,l2)=gi3(l2)*gj3(m2)*Dz(m2,l2)+gi2(l2)*gj2(m2)*Rx
    end
end

%Calculating Ez
Ez = gaz.*(Dz-Iz);
Iz = Iz + gbz.*Ez;


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
    










