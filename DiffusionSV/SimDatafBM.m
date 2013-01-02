function SimDatafBM(data_file,N,step,Vol,H,sigma_X)

    % note: sigma is a function, that computes the volatility term of the
    % price
    
    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,H);
    X = Bh_to_X(Bh,sigma_X);
    Y = SampleObs(X,step,Vol);
    
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
