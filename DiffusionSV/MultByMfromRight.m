function W = MultByMfromRight(N,Z)

% multiplying by M from the rightin O(N) operations

N2 = N;

W = zeros(2*N2,1);
W(1,1) = Z(1);
W(N2+1,1) = Z(N2+1);
for i = 2:N2
    W(N2+i,1) = 1i/sqrt(2) * (Z(i)    - Z(2*N2-i+2));
    W(i,1)    = 1/sqrt(2)  * (Z(i)    + Z(2*N2-i+2));
end
