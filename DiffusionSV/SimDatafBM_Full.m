function SimDatafBM_Full(data_file,N,step,Vol,H,sigma_X,mu,rho,kappa)

    % note: sigma is a function, that computes the volatility term of the
    % price
    
    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,H);
    X = Bh_to_X_Full(Bh,H,step,sigma_X,kappa);
    Y = SampleObs_Full(X,Bh,H,step,Vol,rho,mu);
    
    nobs = length(Y);
    obsstep = N/(nobs-1);
    
    Data = struct();
    Data.X = X;
    Data.Y = Y;
    Data.N = N;
    Data.step = step;
    Data.Htrue = H;
    Data.Ztrue = Z;
    Data.sigma_Xtrue = sigma_X;
    Data.nobs = nobs;
    Data.obsstep = obsstep;
    Data.step = step;
    
    SavePath = '/Users/dureaujoseph/Documents/PhD_Data/fBM/';
    save([SavePath data_file],'Data');
