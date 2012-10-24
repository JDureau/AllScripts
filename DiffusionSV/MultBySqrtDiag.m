function W = MultVectBySqrtDiag(A,B)

% W : sqrt(A) * B where A is diagonal and B is a vector

N = size(A,1);
W = zeros(N,1);

% tmp = chol(A);
for i = 1:N
    W(i) = sqrt(A(i,i))*B(i);
end