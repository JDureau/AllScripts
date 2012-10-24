function Bh = X_to_Bh(X)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


Bh = diff([0 X]);% accounting for X(0) = 0;