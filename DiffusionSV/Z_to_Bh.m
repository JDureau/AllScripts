function Bh = Z_to_Bh(Z,N,step,H)

    Lambda = ComputeLambda(N,step,H);
    Z2 = MultByM(N,Z);
    Z2 = MultBySqrtDiag(Lambda,Z2);
    Bh = ifft(Z2);
    Bh = Bh(1:N);

