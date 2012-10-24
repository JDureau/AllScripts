function W = MultByM(nbpoints_discr,Z)

% multiplyign by M in O(N) operations

N = nbpoints_discr;

W = zeros(2*N,1);
W(1,1) = Z(1);
W(N+1,1) = Z(N+1);
for i = 2:N
    W(i,1)       = 1/sqrt(2) * (Z(i) + 1i*Z(i+N));
    W(2*N-i+2,1) = 1/sqrt(2) * (Z(i) - 1i*Z(i+N));
end
