function M = ComputeM(N)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;

N2 = N-1;
M = zeros(2*N2,2*N2);
M(1,1) = 1;
M(N2+1,N2+1)=1;

for j = 2:N2
    M(j,j)=1/sqrt(2);
    M(j,N2+j)=1i/sqrt(2);
end
for j = 1:N2-1
    M(2*N2-j+1,j+1)=1/sqrt(2);
    M(2*N2-j+1,N2+j+1)=-1i/sqrt(2);
end
