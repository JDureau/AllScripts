function [] = LaunchingfBMExperiments(ind)


cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/DiffusionSV'])

SavePath = '/users/ecologie/dureau/src/AllData/fBM/';

hs = 0.5:0.05:0.95;
Par.H.Value = hs(ind);

SimSeries = 'TestingHidentifiability';

nobs = 250; % Y(0) = 0 is counted as an observation
step = 0.05;
N = (nobs-1)/step;
Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative

Qsampler='HybridMC';
theta_sampler='JointHMC'; % JointHMC or GibbsRW
Par.loop = 100;
Par.Qsampler = Qsampler;
Par.nsteps = 10;
Par.h=0.015;


Par.sigma_X.Value = 0.08;
Par.rho.Value = -0.1;
Par.mu_Y.Value = -0.0014;
Par.mu_X.Value = 0.1;
Par.X0.Value = 0.8;
Par.kappa.Value = 0.027;

Par.H.MinLim = 0;
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

Data = SimDatafBM_Full(N,step,Vol,Par);
Res = RunJointMCMC_Full(Data,Par);
save([SavePath '/' SimSeries '_' num2str(ind)],'Res')




