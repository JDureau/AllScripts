function Z = Sample_Z(N)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;
    
    
    Z = randn(2*(N-1),1);