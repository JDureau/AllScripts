function Data = SimDatafBM_Full(N,step,Vol,Par)

    % note: sigma is a function, that computes the volatility term of the
    % price
    
    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,Par);
    X = Bh_to_X_Full(Bh,step,Par);
    Y = SampleObs_Full(X,Bh,step,Vol,Par);
    
    nobs = length(Y);
    obsstep = N/(nobs-1);
    
    Data = struct();
    Data.Z = Z;
    Data.X = X;
    Data.Y = Y;
    Data.N = N;
    Data.step = step;
    Data.ParTrue = Par;
    Data.nobs = nobs;
    Data.obsstep = obsstep;
    Data.step = step;
    
    
    
    SavePath = '/Users/dureaujoseph/Documents/PhD_Data/fBM/';
%     save([SavePath data_file],'Data');
