function M = ComputeM(nbpoints_discr)

% not very confident with this. Dive a bit more into formulas

N = nbpoints_discr;
M = zeros(2*N,2*N);
M(1,1) = 1;
M(N+1,N+1)=1;

for j = 2:N
    M(j,j)=1/sqrt(2);
    M(j,N+j)=1i/sqrt(2);
end
for j = 1:N-1
    M(2*N-j+1,j+1)=1/sqrt(2);
    M(2*N-j+1,N+j+1)=-1i/sqrt(2);
end
