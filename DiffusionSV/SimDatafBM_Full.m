function Data = SimDatafBM_Full(N,step,Vol,Par)

    % note: sigma is a function, that computes the volatility term of the
    % price
    
    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,Par);
    X = Bh_to_X_Full(Bh,step,Par);
    Obss = SampleObs_Full(X,Bh,step,Vol,Par);
    
    nobs = length(Obss.Y);
    obsstep = N/(nobs-1);
    
    
    Data = struct();
    Data.Z = Z;
    Data.X = X;
    Data.Y = Obss.Y;
    try
        Data.Yx = Obss.Yx;
    end
    Data.Obss = Obss; 
    Data.N = N;
    Data.step = step;
    Data.ParTrue = Par;
    Data.nobs = nobs;
    Data.obsstep = obsstep;
    Data.step = step;
    
    [LogLik LogTerm1] = ComputeLogLikZ_Full(Z,Obss,Vol,Par);
    Data.LogTerm1 = LogTerm1;
    Data.LogLik = LogLik;
    
    SavePath = '/Users/dureaujoseph/Documents/PhD_Data/fBM/';
%     save([SavePath data_file],'Data');
