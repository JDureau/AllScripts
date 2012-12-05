function Bh = Z_to_Bh(Z,N,step,H)

    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;

   
%     Lambda = ComputeLambda(N,step,H)';
    DiagOfLambda = ComputeDiagOfLambda(N,step,H)';
    Z2 = MultByM(N,Z);
    Z2 = sqrt(real(DiagOfLambda)).*Z2;
%     Z2 = MultBySqrtDiag(Lambda,Z2);
    Bh = ifft(Z2);
    Bh = Bh(1:(N-1));

