% validating the code for fractional BM with the full fractional volatility
% model





cd('/Users/dureaujoseph/AllScripts/')
addpath([pwd '/DiffusionSV/'])

SavePath = '/Users/dureaujoseph/Documents/PhD_Data/fBM';

%% sampling from BM - comparing with Wikipedia plots$

N = 4000;
step = 100/N/253;
Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative

Par.H.Value = 0.6;
Par.sigma_X.Value = 0.08;
Par.rho.Value = 0;
Par.mu_Y.Value = -0.0014;
Par.mu_X.Value = -0.0014;
Par.X0.Value = 0;
Par.kappa.Value = 0.023;
Par.obstep = 1/253;
Par.loop=20000;

Z = Sample_Z(N);
Bh = Z_to_Bh(Z,N,step,Par);
X = Bh_to_X_Full(Bh,step,Par);
Obss = SampleObs_Full(X,Bh,step,Vol,Par);
subplot(3,1,1)
plot(X)
ylabel('Vol')
subplot(3,1,2)
plot(Obss.Y)
ylabel('Asset Price')
subplot(3,1,3)
plot(Obss.Yx)
ylabel('Volatility')



%% validating score computation with numerical computation

Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative
nobs = 250;
step = 0.05/253;
N = nobs/step/253;
Par.H.Value = 0.6;
Par.sigma_X.Value = 0.08;
Par.rho.Value = -0.1;
Par.mu_Y.Value = -0.0014;
Par.mu_X.Value = 0.1;
Par.X0.Value = 0.8;
Par.kappa.Value = 0.027;
Par.tau.Value = 0.1;
Par.theta_sampler='JointHMC'; % JointHMC or GibbsRW


Par.Names.All = {'tau','H','sigma_X','mu_Y','rho','kappa','mu_X','X0'};
Par.tau.MinLim = 0.4;
Par.tau.MaxLim = 1;
Par.tau.Transf = @mylog;
Par.tau.InvTransf = @invlog;
Par.tau.Corr = @logCorr;
Par.tau.CorrDer = @logCorrDer;
Par.tau.TransfType = 'Log';
Par.H.MinLim = 0.4;
Par.H.MaxLim = 1;
Par.H.Transf = @logit;
Par.H.InvTransf = @invlogit;
Par.H.Corr = @logitCorr;
Par.H.CorrDer = @logitCorrDer;
Par.H.TransfType = 'Logit';
Par.sigma_X.MinLim = 0;
Par.sigma_X.MaxLim = 0.5;
Par.sigma_X.Transf = @logit;
Par.sigma_X.InvTransf = @invlogit;
Par.sigma_X.Corr = @logitCorr;
Par.sigma_X.CorrDer = @logitCorrDer;
Par.sigma_X.TransfType = 'Logit';
Par.mu_Y.MinLim = -2;
Par.mu_Y.MaxLim =  2;
Par.mu_Y.Transf = @logit;
Par.mu_Y.InvTransf = @invlogit;
Par.mu_Y.Corr = @logitCorr;
Par.mu_Y.CorrDer = @logitCorrDer;
Par.mu_Y.TransfType = 'Logit';
Par.mu_X.MinLim = -10;
Par.mu_X.MaxLim =  10;
Par.mu_X.Transf = @logit;
Par.mu_X.InvTransf = @invlogit;
Par.mu_X.Corr = @logitCorr;
Par.mu_X.CorrDer = @logitCorrDer;
Par.mu_X.TransfType = 'Logit';
Par.X0.MinLim = -10;
Par.X0.MaxLim =  10;
Par.X0.Transf = @logit;
Par.X0.InvTransf = @invlogit;
Par.X0.Corr = @logitCorr;
Par.X0.CorrDer = @logitCorrDer;
Par.X0.TransfType = 'Logit';
Par.rho.MinLim = -1;
Par.rho.MaxLim = 0;
Par.rho.Transf = @logit;
Par.rho.InvTransf = @invlogit;
Par.rho.Corr = @logitCorr;
Par.rho.CorrDer = @logitCorrDer;
Par.rho.TransfType = 'Logit';
Par.kappa.MinLim = 0;
Par.kappa.MaxLim = 10;
Par.kappa.Transf = @logit;
Par.kappa.InvTransf = @invlogit;
Par.kappa.Corr = @logitCorr;
Par.kappa.CorrDer = @logitCorrDer;
Par.kappa.TransfType = 'Logit';

% SimDatafBM_Full('scoreTest.mat',N,step,Vol,Par);
for k = 1:length(Par.Names.All)
   Par.(Par.Names.All{k}).Estimated = 1;
end

Par.obsx = 1;
Par.tau.Estimated = 1;
Par = DefineIndexes(Par);
Par = NoTransfToTransf(Par);


ratios = [];
inds = [];
Par.GradCorr = 1;
Par.Prior = 1;

for i = 1:50
    % sample data

    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,Par);
    X = Bh_to_X_Full(Bh,step,Par);
    Obss = SampleObs_Full(X,Bh,step,Vol,Par);

    % dL/dZ 
    ind = ceil(rand(1,1)*N);
    inds(i) = ind;
    LogLik1 = ComputeLogLikZ_Full(Z,Obss,Vol,Par);
    [LogLik LogTerm1 LogTerm2 LogTerm31] = ComputeLogLikZ_Full(Z,Obss,Vol,Par);
    LogPrior1 = ComputeLogPriorZ_Full(Par);
    epsil = 0.000001;    
    Z2 = Z;
    Z2(ind) = Z2(ind) + epsil;
    LogLik2 = ComputeLogLikZ_Full(Z2,Obss,Vol,Par);
    [LogLik LogTerm1 LogTerm2 LogTerm32] = ComputeLogLikZ_Full(Z2,Obss,Vol,Par);
    LogPrior2 = ComputeLogPriorZ_Full(Par);
    ScoreNumerical = (LogLik2 + LogPrior2 - LogLik1 - LogPrior1)/epsil;
    ScoreAnalytic = ComputeScore_Full(Z,Obss,Vol,VolDer,Par);
    disp(['ratio between the two gradient estimates (comp ' num2str(ind) '): ' num2str(ScoreAnalytic(ind)/ScoreNumerical,10)]);
    ratios(i)=ScoreAnalytic(ind)/ScoreNumerical;
    if isnan(ratios(i))
        die
    end
end

clf
[fi,xi] = ksdensity(ratios);
plot(xi,fi)
title(['Analytic/Numerical Gradient ratios (Z)'],'FontSize',20)
disp([median(ratios) std(ratios)])


% parameters

Par.obsx = 1;

Par.tau.Value = 0.1;
Par.H.Value = 0.6;
Par.sigma_X.Value = 0.08;
Par.rho.Value = -0.1;
Par.mu_Y.Value = -0.0014;
Par.mu_X.Value = 0.1;
Par.X0.Value = 0.8;
Par.kappa.Value = 0.027;
Par.GradCorr = 1;
Par.Prior = 1;
Par.theta_sampler='JointHMC'; 
Par = NoTransfToTransf(Par);

Par.Names.All = {'tau','H','sigma_X','mu_Y','rho','kappa','mu_X','X0'};
% Par.Names.All = {'H','sigma_X','mu_Y','rho','kappa','mu_X','X0'};

for k = 3:length(Par.Names.All)
    for k2 = 1:length(Par.Names.All)
        Par.(Par.Names.All{k2}).Estimated = 0;
    end
    Par.(Par.Names.All{k}).Estimated = 1;
    Par = DefineIndexes(Par);
    Par = NoTransfToTransf(Par);
    ratios = [];
  
    for i = 1:20
        
        Z = Sample_Z(N);
        Bh = Z_to_Bh(Z,N,step,Par);
        X = Bh_to_X_Full(Bh,step,Par);
        Y = SampleObs_Full(X,Bh,step,Vol,Par);
    
        
        ScoreAnalytic = ComputeScore_Full(Z,Y,Vol,VolDer,Par);

        LogLik1 = ComputeLogLikZ_Full(Z,Y,Vol,Par);
        LogPrior1 = ComputeLogPriorZ_Full(Par);
        epsil = 0.000001;    
        Par2 = Par;
        Par2.(Par.Names.All{k}).TransfValue = Par.(Par.Names.All{k}).TransfValue + epsil;
        Par2 = TransfToNoTransf(Par2);
%         Par2.(Par.Names.All{k}).Value = Par.(Par.Names.All{k}).Value + epsil;
%         Par2 = NoTransfToTransf(Par2);
        LogLik2 = ComputeLogLikZ_Full(Z,Y,Vol,Par2);
        LogPrior2 = ComputeLogPriorZ_Full(Par2);
        ScoreNumerical = (LogLik2 + LogPrior2 - LogLik1 - LogPrior1)/epsil;
%         ScoreNumerical = (LogLik2  - LogLik1 )/epsil;
        disp(['ratio between the two gradient estimates: ' num2str(ScoreAnalytic(end)/ScoreNumerical,10)]);
        ratios(i)=ScoreAnalytic(length(Z)+Par.(Par.Names.All{k}).Index)/ScoreNumerical;
    end
    
    clf
    [fi,xi] = ksdensity(ratios);
    plot(xi,fi)
    title(['Density of numerical/analytic' Par.Names.All{k} '- gradient ratios'],'FontSize',20)
    
    pause()
end

    
    
%% playing with HMC


%Create data and likelihood components objects etc
nobs = 250; % Y(0) = 0 is counted as an observation
step = 0.05/253;
N = (nobs-1)/step/253;
Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative


%%% MODES SELECTION

% Qsampler='MALA';
Qsampler='HybridMC';

theta_sampler='JointHMC'; % JointHMC or GibbsRW
% theta_sampler='GibbsRW';
% theta_sampler='GibbsHMC';
% theta_sampler='Blocks';

Par.thetafixed = 0;
Par.Zfixed = 0;
% PARAMETERS

Par.obsx = 0;

Par.tau.Value = 0.1;
Par.H.Value = 0.5;
Par.sigma_X.Value = 0.14%2;
Par.rho.Value = -0.76;
Par.mu_Y.Value = 0.246;
Par.mu_X.Value = -3.28;
Par.X0.Value = -3.28;
Par.kappa.Value = 0.03;%0.6;
Par.Names.All = {'tau','H','sigma_X','mu_Y','rho','kappa','mu_X','X0'};

Par.tau.Transf = @mylog;
Par.tau.InvTransf = @invlog;
Par.tau.Corr = @logCorr;
Par.tau.CorrDer = @logCorrDer;
Par.tau.TransfType = 'Log';
Par.H.MinLim = 0.4;
Par.H.MaxLim = 1;
Par.H.Transf = @logit;
Par.H.InvTransf = @invlogit;
Par.H.Corr = @logitCorr;
Par.H.CorrDer = @logitCorrDer;
Par.H.TransfType = 'Logit';
Par.sigma_X.MinLim = 0;
Par.sigma_X.MaxLim = 20;
Par.sigma_X.Transf = @logit;
Par.sigma_X.InvTransf = @invlogit;
Par.sigma_X.Corr = @logitCorr;
Par.sigma_X.CorrDer = @logitCorrDer;
Par.sigma_X.TransfType = 'Logit';
Par.mu_Y.MinLim = -2;
Par.mu_Y.MaxLim =  2;
Par.mu_Y.Transf = @logit;
Par.mu_Y.InvTransf = @invlogit;
Par.mu_Y.Corr = @logitCorr;
Par.mu_Y.CorrDer = @logitCorrDer;
Par.mu_Y.TransfType = 'Logit';
Par.mu_X.MinLim = -20;
Par.mu_X.MaxLim =  20;
Par.mu_X.Transf = @logit;
Par.mu_X.InvTransf = @invlogit;
Par.mu_X.Corr = @logitCorr;
Par.mu_X.CorrDer = @logitCorrDer;
Par.mu_X.TransfType = 'Logit';
Par.X0.MinLim = -20;
Par.X0.MaxLim =  20;
Par.X0.Transf = @logit;
Par.X0.InvTransf = @invlogit;
Par.X0.Corr = @logitCorr;
Par.X0.CorrDer = @logitCorrDer;
Par.X0.TransfType = 'Logit';
Par.rho.MinLim = -1;
Par.rho.MaxLim = 1;
Par.rho.Transf = @logit;
Par.rho.InvTransf = @invlogit;
Par.rho.Corr = @logitCorr;
Par.rho.CorrDer = @logitCorrDer;
Par.rho.TransfType = 'Logit';
Par.kappa.MinLim = 0;
Par.kappa.MaxLim = 20;
Par.kappa.Transf = @logit;
Par.kappa.InvTransf = @invlogit;
Par.kappa.Corr = @logitCorr;
Par.kappa.CorrDer = @logitCorrDer;
Par.kappa.TransfType = 'Logit';
for k2 = 1:length(Par.Names.All)
    Par.(Par.Names.All{k2}).Estimated = 1;
end
Par.rho.Estimated = 1;
Par.H.Estimated = 1;
Par.tau.Estimated = 0;
Par = DefineIndexes(Par);
Par = NoTransfToTransf(Par);

clf
Data = SimDatafBM_Full(N,step,Vol,Par);
plot(Data.X)
 

if strcmp(Qsampler,'MALA')
    Par.nsteps = 1;
else
    Par.nsteps = 10;
end
Par.theta_sampler='JointHMC'


Par.nsteps = 1;
Par.loop = 1000;
Par.h = 0.01;
Par.NbZpar = 0;
Par.RManif = 0;
Res = RunJointMCMC_Full(Data,Par);
Data.ParStart = Res.Par;
Data.Z = Res.Z;
Par = Data.ParTrue;
Par.nsteps = 20;
Data.Cov = cov(Res.TransfThetas');
Par.loop = 500;
Par.h = 0.07;
Res = RunJointMCMC_Full(Data,Par);


Par.loop = 500;
Par.Qsampler = Qsampler;
if strcmp(Qsampler,'MALA')
    Par.nsteps = 1;
else
    Par.nsteps = 10;
end

theta_sampler='JointHMC'
Par = DefineIndexes(Par);
Par = NoTransfToTransf(Par);

Par.theta_sampler = theta_sampler;
if strcmp(theta_sampler,'JointHMC')
    Par.h=0.004;
    Par.h=0.2;
elseif strcmp(theta_sampler,'GibbsHMC')
    Par.hZ=1;%0.08;
    Par.hP=0.11;
elseif strcmp(theta_sampler,'GibbsRW')
    Par.hH = 0.5;
    Par.hZ = 0.1;
    Par.hmu = 0.1;
    Par.hrho = 0.1;
    Par.hkappa = 0.1;
    Par.hsig = 1;
elseif strcmp(theta_sampler,'Blocks')
    Par.d = 10;
end


Res = RunJointMCMC_Full(Data,Par)
PlotfBMoutput(Res)

Data.Cov = diag(diag(cov(Res.TransfThetas')));


save([SavePath '/testingHMC.mat'],'Res')


load([SavePath '/Test.mat'])
save([SavePath '/Test.mat'],'Data')


load([SavePath '/D' num2str(3) '_Exp' num2str(1) '_' num2str(6) '_' num2str(1) '_' num2str(1)])
PlotfBMoutput(Res)
Cov = cov(Res.TransfThetas');
Data = Res.Data;
Par.H.Estimated=0;
Par = DefineIndexes(Par);
Par = NoTransfToTransf(Par);


% Data = SimDatafBM_Full(N,step,Vol,Par);
% Par = Data.ParTrue;
Par.Epsil = 1;
Par.MCMCType = 'Rand';
Par.G = eye(length(Par.Names.Estimated));
% Par.G = Cov^(-1);
Par.ModelType='SMC';
Par.NbVariables = 3;
Par.NbParticules = 200;
Par.NoPaths = 0;
Par.PathsToKeep = [1];
Par.NbParsEstimated  =length(Par.Names.Estimated);
Par.ComputationTStep = Data.step;
Par.Vol = Vol;
Par.Problem = 'vol';
Data.ObservedVariables = 1;
Par.AdaptC = 0.999;
Par.GMeth =  'cst given';
Data.NbComputingSteps = [0 Data.obsstep*ones(1,Data.nobs)] ;
fullvolModel.InitializeParameters = @fullvolInitialize;
fullvolModel.SMC_projection = @fullvol_SMC_projection;
fullvolModel.LikFunction = 'normpdf(Data.Y(IndTime)-Data.Y(IndTime-1),Variables(:,2),sqrt(Variables(:,3)))';
TempPar = ProposeInitialParameter(Data, fullvolModel, Par);
Res = RunEstimationMethod(Data, fullvolModel, Par, TempPar, 20000);
PlotfBMpmcmc(Res)

ns = 600:100:1000;
Ress = {};
for i = 1:length(ns)
    Par.NbParticules = 9000;%ns(i);
    Par.G = Cov^(-1);
    Par.AdaptC = 0;
    Par.NoPaths = 1;
    Res = RunEstimationMethod(Data, fullvolModel, Par, TempPar, 1000);
    Ress{i} = Res;
end

save([SavePath '/D3_pMCMC_' num2str(Par.NbParticules) '_Cov0' '.mat'],'Res');


Par.G = cov(Res.TransfThetas')^(-1);
TempPar = Res.TempPar;
Res = RunEstimationMethod(Data, fullvolModel, Par, TempPar, 10000);
PlotfBMpmcmc(Res)


Data.Cov = cov(Res.TransfThetas')^(-1);

Par.NbParticules = 100;
Res = EstimationSMCsmoothGen(Data, fullvolModel, Par);
Res.LogLik

clf
xis = step:step:nobs-1-step;
try
    plot(Data.X(obsstep:obsstep:N-1),'g')
catch
    q50 = quantile(squeeze(Res.CompletePaths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.5);
end    
hold on
try
    q2p5 = quantile(squeeze(Res.CompletePaths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.025);
    q25 = quantile(squeeze(Res.CompletePaths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.25);
    q50 = quantile(squeeze(Res.CompletePaths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.5);
    q75 = quantile(squeeze(Res.CompletePaths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.75);
    q97p5 = quantile(squeeze(Res.CompletePaths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.975);
    plot(q2p5,':')
    plot(q25,'--')
    plot(q50,'-')
    plot(q75,'--')
    plot(q97p5,':')
end
legend('True volatility','95% c.i.','50% c.i.','Posterior median')
hold off
ylabel('Volatility','FontSize',12)
xlabel('Time','FontSize',12)





Data.ParTrue = Par;
save([SavePath '/DataSet3.mat'],'Data')


load([SavePath '/DataSet4.mat'])



load([SavePath '/StateBug.mat'])
Data = state.Data;
Par = state.Par;
Data.ParTrue = Par;
Data.Z  = state.Z;
Par.loop = 200;
Par.h = 0.02;
Par.nsteps = 10;
Res = RunJointMCMC_Full(Data,Par)
PlotfBMoutput(Res)



save([SavePath '/Data_Hsims_0.9.mat'],'Data')


load([SavePath 'Res_JointHMC10_Allest_fBMData_H=0.6_sigma=0.5_mu=0.1_rho=0.3_kappa=-0.2.mat'])


N = 1000;
step = 0.01;
Par.H.Value = 0.5;
Z = Sample_Z(N);
Bh1 = Z_to_Bh(Z,N,step,Par);
Par.H.Value = 0.75;
Z = Sample_Z(N);
Bh2 = Z_to_Bh(Z,N,step,Par);
Par.H.Value = 0.25;
Z = Sample_Z(N);
Bh3 = Z_to_Bh(Z,N,step,Par);

clf
subplot(1,3,1)
plot(Par.sigma_X.Value*cumsum(Bh2),'k')
title('H=0.25','FontSize',20)
xlabel('time','FontSize',16)
ylabel('B^H_t','FontSize',16)
subplot(1,3,2)
plot(Par.sigma_X.Value*cumsum(Bh1),'k')
title('H=0.5','FontSize',20)
xlabel('time','FontSize',16)
ylabel('B^H_t','FontSize',16)
subplot(1,3,3)
plot(Par.sigma_X.Value*cumsum(Bh3),'k')
title('H=0.75','FontSize',20)
xlabel('time','FontSize',16)
ylabel('B^H_t','FontSize',16)




data_file=strcat('fBMData_H=',num2str(Htrue),'_sigma=',num2str(sigma_Xtrue),'_mu=',num2str(mutrue),'_rho=',num2str(rhotrue),'_kappa=',num2str(kappatrue),'.mat');

data_files = {'fBMData_H=0.6_sigma=0.1_mu=0.1_rho=0.9_kappa=0.1.mat'}

% DATA
SimDatafBM_Full(data_file,N,step,Vol,Par);

SavePath = '/Users/dureaujoseph/Documents/PhD_Data/fBM/';


% % datafiles = {'fBMData_H=0.6_sigma=0.1.mat','fBMData_H=0.6_sigma=0.01.mat','fBMData_H=0.8_sigma=0.1.mat'};
% % method = {'MALA','HMC'};
% % 
% % indData = 1;
% load([SavePath '/' datafiles{indData}]);

load([SavePath '/' data_file]);


Par = struct();
Par.Names.All = {'H','sigma_X','mu_Y','rho','kappa','mu_X','X0'};
Par.H.MinLim = 0.4;
Par.H.MaxLim = 1;
Par.H.Transf = @logit;
Par.H.InvTransf = @invlogit;
Par.H.Corr = @logitCorr;
Par.H.Estimated = 0;
Par.sigma_X.MinLim = 0;
Par.sigma_X.MaxLim = 100;
Par.sigma_X.Transf = @mylog;
Par.sigma_X.InvTransf = @invlog;
Par.sigma_X.Corr = @logCorr;
Par.sigma_X.Estimated = 0;
Par.mu.MinLim = -1;
Par.mu.MaxLim =  1;
Par.mu.Transf = @logit;
Par.mu.InvTransf = @invlogit;
Par.mu.Corr = @logitCorr;
Par.mu.Estimated = 0;
Par.rho.MinLim = -1;
Par.rho.MaxLim = 1;
Par.rho.Transf = @logit;
Par.rho.InvTransf = @invlogit;
Par.rho.Corr = @logitCorr;
Par.rho.Estimated = 0;
Par.kappa.MinLim = -1;
Par.kappa.MaxLim = 1;
Par.kappa.Transf = @logit;
Par.kappa.InvTransf = @invlogit;
Par.kappa.Corr = @logitCorr;
Par.kappa.Estimated = 0;
Par = DefineIndexes(Par);
Par.H.Value = Data.Htrue;
Par.sigma_X.Value = Data.sigma_Xtrue;
Par.mu.Value = Data.mutrue;
Par.rho.Value = Data.rhotrue;
Par.kappa.Value = Data.kappatrue;

for i = 1:length(Par.Names.All)
    Par.(Par.Names.All{i}).Estimated = 1;
end
Par = DefineIndexes(Par);
Par = NoTransfToTransf(Par);


Par.loop = loop;
Par.Qsampler = Qsampler;
if strcmp(Qsampler,'MALA')
    Par.nsteps = 1;
else
    Par.nsteps = 10;
end

Par.theta_sampler = theta_sampler;
if strcmp(theta_sampler,'JointHMC')
    Par.h=0.0001;
elseif strcmp(theta_sampler,'GibbsHMC')
    Par.hZ=0.19;
    Par.htheta=0.13;
elseif strcmp(theta_sampler,'GibbsRW')
    Par.hH = 0.5;
    Par.hZ = 0.1;
    Par.hmu = 0.1;
    Par.hrho = 0.1;
    Par.hkappa = 0.1;
    Par.hsig = 1;
elseif strcmp(theta_sampler,'Blocks')
    Par.d = 10;
end



Res = RunJointMCMC_Full(Data,Par)


save([SavePath '/Res_JointHMC10_Allest_' data_file],'Res')


est = 'No';
load([SavePath '/Res_RW_' est 'est_' data_file])
ResRW = Res;
load([SavePath '/Res_MALA_' est 'est_' data_file])
ResMALA = Res;
load([SavePath '/Res_HMC5_' est 'est_' data_file])
ResHMC5 = Res;
load([SavePath '/Res_HMC10_' est 'est_' data_file])
ResHMC10 = Res;
load([SavePath '/Res_HMC20_' est 'est_' data_file])
ResHMC20 = Res;

clf
[fi,xi] = ksdensity(ResRW.out_Hs);
plot(xi,fi,'k')
hold on
[fi,xi] = ksdensity(ResMALA.out_Hs);
plot(xi,fi,'g')
[fi,xi] = ksdensity(ResHMC.out_Hs);
plot(xi,fi,'b')
hold off
ylabel('p(H|y)','FontSize',20)
xlabel('H','FontSize',20)
h = legend('Gibbs-RW: ESS=0.4%','Joint-MALA: ESS=0.5%','Joint-HMC: ESS=6.4%');
set(h,'FontSize',20)

load([SavePath 'Res_JointHMC10_Allest_fBMData_H=0.6_sigma=0.5_mu=0.1_rho=0.3_kappa=-0.2.mat'])

Ress = {ResRW,ResMALA,ResHMC5,ResHMC10,ResHMC20};
Ress = {Res,Res};
ESSs = {};

for k = 2:2
    k
%     Ress{k}.h
    %figure(11); for i=10:10:min(500,loop); plot(out_Q(i,:),'b');hold on; end; plot(out_Q(1,:),'r');hold off
    ESS=zeros(size(Ress{k}.out_Zs,2),1);
    for i=1:size(Ress{k}.out_Zs,2)
        r=sum(autocorr(Ress{k}.out_Zs(:,i),100));
        ESS(i)=100/(1+2*r);
    end
    ESSs{k} = ESS;
    
    figure(12);plot(ESS);title('ESS (%) for each v over time')   
    disp([min(ESS),median(ESS),max(ESS)])   
    
    Names = Res.Par.Names.Estimated;
    for i = 1:length(Names)
        ind = Res.Par.(Names{i}).Index;
        r=sum(autocorr(Ress{k}.Thetas(ind,:),1200));
        disp(100/(1+2*r));
    end
%     r=sum(autocorr(Ress{k}.out_Hs,1200));
%     disp(100/(1+2*r));
%     r=sum(autocorr(Ress{k}.out_sigs,1200));
%     disp(100/(1+2*r));
end
    
plot(ESSs{1},'k')
hold on
plot(ESSs{2},'g')
plot(ESSs{3},'b')
hold off
ylabel('ESS(Z_t)','FontSize',20)
xlabel('t','FontSize',20)
h = legend('Gibbs-RW: ESS_{min}=0.1%','Joint-MALA: ESS_{min}=0.1%','Joint-HMC: ESS_{min}=1.1%');
set(h,'FontSize',20)




load([SavePath 'Res_MALA_Hsigest_fBMData1000_H=0.6_sigma=0.1.mat'])
Res1000 = Res;
load([SavePath 'Res_JointHMC_Hsigest_fBMData_H=0.6_sigma=0.1.mat'])

clf
subplot(3,1,3)
[fi,xi]=ksdensity(Res.out_Hs);
plot(xi,fi);
hold on
[fi,xi]=ksdensity(Res1000.out_Hs);
plot(xi,fi,'g');
hold off
title('H','FontSize',20)
subplot(3,1,3)
[fi,xi]=ksdensity(Res.out_sigs);
plot(xi,fi);
hold on
[fi,xi]=ksdensity(Res1000.out_sigs);
plot(xi,fi,'g');
hold off
title('\sigma_X','FontSize',20)
h = legend('250 observations','1000 observations');
set(h,'FontSize',20)
subplot(3,1,3)
plot(Res.out_Hs,Res.out_sigs,'.')
hold on
plot(Res1000.out_Hs,Res1000.out_sigs,'.g')
hold off
xlabel('H')
ylabel('\sigma_X')

% r=sum(autocorr(out_kappa,100));
% ESSk=100/(1+2*r);
% %disp(ESSk)
% disp([min(ESS),median(ESS),max(ESS),ESSk])

