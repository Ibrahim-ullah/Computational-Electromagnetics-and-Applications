%specifying parameters

lx =10;      %length in x-direction
ly=10;       %length in y-direction
nx=100;       %number of steps in space(x)
ny=100;       %number of steps in space(y)
niter=1000;        %number of iterations
dx=lx/(nx-1);   %width of space step(x)
dy=ly/(ny-1);   %width of space step(y)
x = 0:dx:lx;    %range of x(0,10) specifying the given range
y = 0:dy:ly;     %range of x(0,10) specifying the given range


%%
%source terms

f= zeros(ny,nx)     %preallocation p
f(round(nx/4),round(ny/4))  = 3000;
f(round(3*nx/4),round(3*ny/4))  = -3000;


%intial conditions
p= zeros(ny,nx)     %preallocation p
pn=zeros(ny,nx);    %preallocating pn


% Boundary conditions (Dirichlet conditions)

 p(:,1) = 0;
 p(:,nx) =0;
 p(1,:)  =0;
 p(ny,:) =0;

%%
%Explicit iterative CD scheme in space (5-point difference). by 5 point
%difference we meant we are taking four variable into consideration to
%calculate variable.

i = 2:nx-1;   %we start from 2 because of the boundary conditions
j = 2:ny-1;

for it = 1:niter
    pn = p;
    p(i,j) = ((dy^2*(pn(i+1,j)+pn(i-1,j)))+ (dx^2*(pn(i,j+1)+pn(i,j-1)))-(f(i,j)*dx^2*dy^2))/(2*(dx^2+dy^2));
       
    
end
 
%Plotting the solution

surf(x,y,p,'Edgecolor','none');
shading interp
title({'2-D Poissons equation';
        ['{\itNumber of iterations} = ',num2str(it)]})
    
xlabel('Spatial co-ordinate (x) \rightarrow')
ylabel('{\leftarrow} Spatial co-ordinate (y)')
zlabel('Solution profile (p) \rightarrow')

view(0,90);




 
 
 
    
    
    
    