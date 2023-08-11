delx = 0.5; % resolution size
r = 1;      %aspect ration
u = 1;      % constant of wave equation
delt = r^2*delx/2;  % time step size
Tsteps = round(1/delt);   % number of time steps



% X1 is the potential grid of the simulation, due to symmetry only half of
% the field is calculated


X1 = zeros(Tsteps, 1/(2*delx)+2);  %initialize X1

%Initial conditions and reflection line defined

x = 0:delx:0.5+delx
X1(1,:) = sin(pi*x);
X1(2,2:end-1) = .5*(X1(1,1:end-2) + X1(1,3:end));
X1(2,end) = X1(2,end-2); %Reflection line


for row = 3:size(X1,1)
    for col = 2:size(X1,2)-1
        X1(row,col) = X1(row-1, col-1) + X1(row-1, col+1) - X1(row-2, col);
    end
    X1(row,end) = X1(row, end-2) %reflect line
end


% Use symmetry condition to create entire field


X2 = [X1, fliplr(X1(:,1:end-3))];
AX2 = zeros(Tsteps,101)


for n = 1:1:200
        for x = 1:1:101
            AX2(n,x) = sin(pi*x/100)*cos(pi*n/100);
        end
end


figure(1), imagesc(0:delx:1,(0:delt:Tsteps*delt),X2), colorbar
ylabel('\leftarrow time (sec)')
xlabel('x')
title('Numerical Hyperbolic PDE')

figure(2), imagesc(0:delx:1,(0:delt:Tsteps*delt),AX2), colorbar
ylabel('\leftarrow time (sec)')
xlabel('x')
title('Numerical Hyperbolic PDE')



if(delx ==.1)
    dispmat = [X1(1:8, 1:7)];
    disp(sprintf('\nCompare to ta, Solution of the wave equation'))
    disp(num2str(dispmat))
end




