function Bh = Z_to_Bh(Z,N,step,Par)

    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;

   H = Par.H.Value;

    if 0%H == 0.5
    %     Lambda = ComputeLambda(N,step,H)';
        Bh = sqrt(step)*Z(1:N);
    else
        DiagOfLambda = ComputeDiagOfLambda(N,step,H)';
        Z2 = MultByM(N,Z);
        Z2 = sqrt(2*N)*sqrt(real(DiagOfLambda)).*Z2; % sqrt(2N) is because Matlab's inverse fft divides by N...
    %     Z2 = MultBySqrtDiag(Lambda,Z2);
        Bh = ifft(Z2);
        Bh = real(Bh(1:(N)));
    end