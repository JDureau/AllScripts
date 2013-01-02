function Y = SimDatafBM(N,step,sigma,H)

    % note: sigma is a function, that computes the volatility term of the
    % price
    
    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,H);
    X = Bh_to_X(Bh);
    Y = SampleObs(X,step,sigma);