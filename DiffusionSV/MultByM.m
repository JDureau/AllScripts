function W = MultByM(N,Z)

% multiplyign by M in O(N) operations

N2 = N;

W = zeros(2*N2,1);
W(1,1) = Z(1);
W(N2+1,1) = Z(N2+1);
% for i = 2:N2
%     
%     W(i,1)       = 1/sqrt(2) * (Z(i) + 1i*Z(i+N2));
%     W(2*N2-i+2,1) = 1/sqrt(2) * (Z(i) - 1i*Z(i+N2));
% end
W(2:N2,1) = 1/sqrt(2) *  (Z(2:N2) + 1i*Z((N2+2):2*N2));
W(2*N2:-1:(N2+2),1) = 1/sqrt(2) *  (Z(2:N2) - 1i*Z((N2+2):2*N2));