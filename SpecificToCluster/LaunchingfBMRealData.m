function [] = LaunchingfBMRealData(ind,MorePars)

ind = ind+1;

cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/DiffusionSV'])

SavePath = '/users/ecologie/dureau/src/AllData/fBM/';
2*log(std(diff(Data.Y)))

% clf
% for i = 1:6
%     Data = LoadYahooData([SavePath '/' files{i}],step);
%     hold on
%     plot(Data.Y-Data.Y(1),'k')
% end
% hold off


DataSeries = 'SP500';

Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative

Par.Qsampler='HybridMC';
Par.theta_sampler='JointHMC'; % JointHMC or GibbsRW
Par.loop = 40000;
Par.nsteps = 10;
Par.h=0.01;

Par.thetafixed = 0;
Par.Zfixed = 0;
% PARAMETERS
Par.H.Value = 0.5;
Par.sigma_X.Value = 2;
Par.rho.Value = 0;
Par.mu_Y.Value = -0.002;
Par.mu_X.Value = -3;
Par.X0.Value = -3;
Par.kappa.Value = 1;
Par.Names.All = {'H','sigma_X','mu_Y','rho','kappa','mu_X','X0'};

Par.H.MinLim = 0.4;
Par.H.MaxLim = 1;
Par.H.Transf = @logit;
Par.H.InvTransf = @invlogit;
Par.H.Corr = @logitCorr;
Par.H.CorrDer = @logitCorrDer;
Par.H.TransfType = 'Logit';
Par.sigma_X.Transf = @mylog;
Par.sigma_X.InvTransf = @invlog;
Par.sigma_X.Corr = @logCorr;
Par.sigma_X.CorrDer = @logCorrDer;
Par.sigma_X.TransfType = 'Log';
Par.mu_Y.Transf = @myid;
Par.mu_Y.InvTransf = @invid;
Par.mu_Y.Corr = @idCorr;
Par.mu_Y.CorrDer = @idCorrDer;
Par.mu_Y.TransfType = 'Id';
Par.mu_X.Transf = @myid;
Par.mu_X.InvTransf = @invid;
Par.mu_X.Corr = @idCorr;
Par.mu_X.CorrDer = @idCorrDer;
Par.mu_X.TransfType = 'Id';
Par.X0.Transf = @myid;
Par.X0.InvTransf = @invid;
Par.X0.Corr = @idCorr;
Par.X0.CorrDer = @idCorrDer;
Par.X0.TransfType = 'Id';
Par.rho.MinLim = -1;
Par.rho.MaxLim = 1;
Par.rho.Transf = @logit;
Par.rho.InvTransf = @invlogit;
Par.rho.Corr = @logitCorr;
Par.rho.CorrDer = @logitCorrDer;
Par.rho.TransfType = 'Logit';
Par.kappa.MinLim = 0;
Par.kappa.MaxLim = 253;
Par.kappa.Transf = @logit;
Par.kappa.InvTransf = @invlogit;
Par.kappa.Corr = @logitCorr;
Par.kappa.CorrDer = @logitCorrDer;
Par.kappa.TransfType = 'Logit';

for k2 = 1:length(Par.Names.All)
   Par.(Par.Names.All{k2}).Estimated = 1;
end
%Par.mu_X.Estimated = 0;
% if not(MorePars)
%     Par.X0.Estimated = 0;
%     Par.mu_X.Estimated = 0;
% end    
Par.rho.Estimated = 1;
Par.H.Estimated = 1;
Par.sigma_X.Estimated = 1;
Par = DefineIndexes(Par);
Par = NoTransfToTransf(Par);

Par.nsteps = 1;
Par.loop = 2000;
Par.h = 0.01;
Par.NbZpar = 0;
Par.RManif = 0;
Res = RunJointMCMC_Full(Data,Par);
Data.ParTrue = Res.Par;
Data.Z = Res.Z;
Par = Data.ParTrue;
Par.nsteps = 10;
Data.Cov = cov(Res.TransfThetas');
Par.loop = 30000;
Par.h = 0.01;
% Data = BostCov([3 4],Data,2000);
tic
Res = RunJointMCMC_Full(Data,Par);
toc


Data.ParTrue = Res.Par;
Data.Z = Res.Z;
save([SavePath '/Data_SP500_' num2str(ind) '_' num2str(MorePars) '.mat'],'Data');

load([SavePath '/Data_SP500_' num2str(ind) '_' num2str(MorePars) '.mat'])

Par = Data.ParTrue;
Par.nsteps = 10;
Par.loop = 20000;
Par.h = 0.002;
Res = RunJointMCMC_Full(Data,Par);

save([SavePath '/' DataSeries '_' num2str(ind) '_' num2str(MorePars) '_HMC.mat'],'Res')



% Par = Res.Par;
% Par.loop = 2000;
% Data.Z = Res.Z;
% Res = RunJointMCMC_Full(Data,Par);
% 
% 
% save([SavePath '/' DataSeries '_' num2str(ind)],'Res')




