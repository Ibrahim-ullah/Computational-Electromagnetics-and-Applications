function k = HFD1D(a, N)

%arguements:
% a = lenght of interval
% N = number of subintervals (equal length)
% Returns
% k = eigenvalue
h = a/N; %grid size
A = spalloc(N-1, N-1,3*(N-1)) % Allocate sparse matrix with 3*(N-1) nonzeros
d = -2/h^2;  %value of diagonal entries
s = 1/h^2; %value of upper and lower diagonal entries

%initialize the diagonal entries
for i = 1:N-1
    A(i,i) = d; %diagonal entries
end

%initialize the upper and lower diagonal entries

for i = 1:N-2
    A(i,i+1)=s; %upper diagonal entries
    A(i+1,i)=s; %lower diagonal entries
end

%computing the eigenvalues
lambda = eig(A);  %the eig function is useful for 1D and 2D problem but not for 3d problems for 3d problems we have to use eigs
k = sqrt(sort(-lambda));

