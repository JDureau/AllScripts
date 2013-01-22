function [] = LaunchingfBMRealData(ind)

ind = ind+1;

cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/DiffusionSV'])

SavePath = '/users/ecologie/dureau/src/AllData/fBM/';


files = {'15sep08_15sep09.csv','15mar07_15mar08.csv'};

step = 0.05;
Data = LoadYahooData([SavePath '/' files{ind}],step);
Par.N = Data.N;
Par.nobs = Data.nobs;
Par.step = Data.step;
Par.npoints = Data.npoints;


DataSeries = 'SP500';

Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative

Par.Qsampler='HybridMC';
Par.theta_sampler='JointHMC'; % JointHMC or GibbsRW
Par.loop = 40000;
Par.nsteps = 10;
Par.h=0.015;

Par.H.Value = 0.5;
Par.sigma_X.Value = 0.08;
Par.rho.Value = -0.1;
Par.mu_Y.Value = -0.0014;
Par.mu_X.Value = 0.1;
Par.X0.Value = 0.8;
Par.kappa.Value = 0.027;
Data.ParTrue = Par;
Data.Z = Sample_Z(Data.N);

Par.Names.All = {'H','sigma_X','mu_Y','rho','kappa','mu_X','X0'};

Par.H.MinLim = 0.4;
Par.H.MaxLim = 1;
Par.H.Transf = @logit;
Par.H.InvTransf = @invlogit;
Par.H.Corr = @logitCorr;
Par.sigma_X.MinLim = 0;
Par.sigma_X.MaxLim = 10;
Par.sigma_X.Transf = @logit;
Par.sigma_X.InvTransf = @invlogit;
Par.sigma_X.Corr = @logitCorr;
Par.mu_Y.MinLim = -10;
Par.mu_Y.MaxLim =  10;
Par.mu_Y.Transf = @logit;
Par.mu_Y.InvTransf = @invlogit;
Par.mu_Y.Corr = @logitCorr;
Par.mu_X.MinLim = -10;
Par.mu_X.MaxLim =  10;
Par.mu_X.Transf = @logit;
Par.mu_X.InvTransf = @invlogit;
Par.mu_X.Corr = @logitCorr;
Par.X0.MinLim = -10;
Par.X0.MaxLim =  10;
Par.X0.Transf = @logit;
Par.X0.InvTransf = @invlogit;
Par.X0.Corr = @logitCorr;
Par.rho.MinLim = -1;
Par.rho.MaxLim = 0;
Par.rho.Transf = @logit;
Par.rho.InvTransf = @invlogit;
Par.rho.Corr = @logitCorr;
Par.kappa.MinLim = 0;
Par.kappa.MaxLim = 1;
Par.kappa.Transf = @logit;
Par.kappa.InvTransf = @invlogit;
Par.kappa.Corr = @logitCorr;
for k2 = 1:length(Par.Names.All)
   Par.(Par.Names.All{k2}).Estimated = 1;
end
%Par.mu_X.Estimated = 0;
Par = DefineIndexes(Par);
Par = NoTransfToTransf(Par);

Res = RunJointMCMC_Full(Data,Par);


save([SavePath '/' DataSeries '_' num2str(ind)],'Res')




