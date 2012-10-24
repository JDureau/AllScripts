function B = Sample_B(nbpoints_discr,step,H)

    Z = Sample_Z(nbpoints_discr);
    Lambda = ComputeLambda(nbpoints_discr,step,H);
    Z2 = MultByM(nbpoints_discr,Z);
    Z2 = MultBySqrtDiag(Lambda,Z2);
    B = ifft(Z2);
    Bh = Bh(1:N);
    
    
    