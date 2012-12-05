function X = Bh_to_X(B,sigma_X)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;

X = sigma_X*cumsum(B);% sigma_X is missing here