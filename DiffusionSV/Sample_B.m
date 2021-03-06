function B = Sample_B(N,step,H)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;
    
    Z = Sample_Z(N);
    Lambda = ComputeLambda(N,step,H);
    Z2 = MultByM(N,Z);
    Z2 = MultBySqrtDiag(Lambda,Z2);
    Bh = ifft(Z2);
    Bh = Bh(1:(N-1));
    
    
    