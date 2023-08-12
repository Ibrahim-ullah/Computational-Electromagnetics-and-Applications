function cap = capacitor(a,b,c,d,n,tol)

% Arguments:
% a width of inner conductor
% b height of inner conductor
% c width of outer conductor
% d height of outer conductor
% n number of points in the x-direction (horizontal)
% tol relative tolerance for capacitance (1e-9)
% rel relaxation parameter = 2-c/n
% value calculated based on input arguments
% Returns cap = capacitance per unit length [pF/m]

%make grids

rel = 2-c/n;    %relaxation parameter
h = 0.5*c/n;    %grid size
na = round(0.5*a/h);    %number of segments on a
x = linspace(0,0.5*c,n+1);  %Grid points along x-axis
m = round(0.5*d/h); %number of segments on d
mb = round(0.5*b/h);    %number of segments on b
y = linspace(0,0.5*d, m+1); %grid points along y-axis


%initialize potential and mask array
f= zeros(n+1,m+1);  %2D array with solution
mask = ones(n+1,m+1)*rel;   %2D-array with relaxation

for i = 1:na+1
    for j = 1:mb+1
        mask(i,j) = 0;
        f(i,j) =1;
    end
end


%Gauss Seidel iteration

oldcap = 0;
for iter = 1:1000 %maximum number of iterations
    f = seidel(f,mask,n,m); %perform gauss-seidel iteration
    cap = gauss(n,m,h,f);   %compute the capacitance
    if (abs(cap-oldcap)/cap<tol)
        break %stop if change in capacitance is sufficiently small
    else
    oldcap = cap; %continue untill converged
    end
end


str = sprintf('Number of iterations =%4i',iter);
disp(str)


%make one seidel iteration

function f = seidel(f,mask,n,m)
% arguments
% f = 2D array with solution
% mask = 2D array with relaxation
% n = number of points in the x-direction
% m = number of points in the y-direction
% Returns:
% f = 2D array with solutions after gauss-seidel iteration
% Gauss seidel iteration

for i =2:n
    for j =2:m
        f(i,j) = f(i,j) + mask(i,j)*(0.25*(f(i-1,j) + f(i+1,j) + f(i,j-1) + f(i,j+1)) - f(i,j));
    end
end

%symmetry on left boundary i-1 -> i+1

i=1;
for j=2:m
    f(i,j) = f(i,j) +mask(i,j)*(0.25*(f(i+1,j) + f(i+1,j) + f(i,j-1) + f(i,j+1)) - f(i,j));
end

%symmetry on lower boundary j-1 -> j+1
j=1
for i=2:n
    f(i,j) = f(i,j) + mask(i,j)*(0.25*(f(i-1,j) + f(i+1,j) + f(i,j+1) + f(i,j+1)) -f(i,j));
end

    
%compute capcitance from the potential

function cap = gauss(n,m,h,f)
% arguments
% n = number of points in the x-direction (horizontal)
% m = number of points in the y-direction (vertical)
% h = cell size
% f = 2D array with solution
% Returns:
% cap = capacitance per unit length (pF/m)

q = 0;
for i = 1:n
    q = q + (f(i,m)+f(i+1,m))*0.5 %integrate alogn upper boundary
end

for j=1:m
    q = q + (f(n,j)+f(n,j+1))*0.5 %integrate alogn upper boundary
end

cap = q*4; % 4 quadrants
cap = cap*8.854187; %epsilon0*1e12 gives answer pF/m










