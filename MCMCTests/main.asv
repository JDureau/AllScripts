% Methodological tests


addpath('H:\My Documents\PhD Work\Matlab Scripts\MCMCTests\')
addpath('H:\My Documents\PhD Work\Matlab Scripts\General Tools\') 
addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes\mcmcdiag\') 



% testing on a GMM(2)





Epss = [ 0.2 0.6  1 1.4 1.8 2];
AccRates = [];
RelESSs = [];
for i = 1:length(Epss)
    i
    Parameters.Epsil = Epss(i);
    Res = RunMCMC(Parameters.ArgMax,Parameters,3000);
    AccRates(i) = Res.AccRate;
    RelESSs(i,:) = Res.RelESS;
end

figure()
subplot(2,1,1)
plot(Epss,AccRates)
subplot(2,1,2)
plot(Epss,RelESSs)


%% Asymptotic results: suppose we have learnt the optimal GMM:

% Cross
dim = 2;
mu1 = zeros(1,dim);
tmp = eye(dim);
tmp(1,1) = 1;
sigma1 = tmp;
mu2 = zeros(1,dim);
tmp = eye(dim);
tmp(2,2) = 1;
sigma2 = tmp;
X = [mvnrnd(mu1,sigma1,10000);mvnrnd(mu2,sigma2,10000)];
scatter(X(:,1),X(:,2),10,'ko')
Parameters.OptDens = gmdistribution.fit(X,2);
Parameters.f = @fGMM;


% Banana
dim = 2;
Parameters.B = 0.15;
Parameters.f = @fBanana;
Parameters.LogRatioFun = @LogRatioGMMrand;
Parameters.SampleFun = @SampleGMMrand;
Parameters = FindFisherInfMat(zeros(Parameters.Dim,1),Parameters) ;
Parameters.CholCov = chol((-Parameters.Hess)^-1);
X = mvnrnd(Parameters.ArgMax,Parameters.CholCov*Parameters.CholCov',1000);
Parameters.Dens = gmdistribution.fit(X,1);
Res = RunMCMC(Parameters.ArgMax,Parameters,100000);

SavePath = 'S:\ResultsMCMC\';
save([SavePath 'BananaSamples_2d_0p15B.mat'],'Res');

plot(Res.Vals(1,:),Res.Vals(2,:),'.')

load([SavePath 'BananaSamples_2d_0p15B.mat'])

X = Res.Vals';
AICs = [];
BICs = [];
for i = 1:15
    i
    tmp = gmdistribution.fit(X,i);
    AICs(i) = tmp.AIC;
    BICs(i) = tmp.BIC;
end
plot(AICs)
hold on
plot(BICs,'k')
hold off

[b,i] = min(AICs)
[b,j] = min(BICs)
% both say 14

Parameters.OptDens = gmdistribution.fit(X,14);
plot(Res.Vals(1,:),Res.Vals(2,:),'.k')
hold on
ezcontour(@(x,y)pdf(Parameters.OptDens,[x y]),[-20 20],[-40 20],500);
hold off

Parameters.Dim = dim;

Epss = [ 0.5  1 1.5 2 2.5 3]; % as chosing eps is hard, we try for a series of values


% GMM1 - Fish
Parameters = FindFisherInfMat(zeros(Parameters.Dim,1),Parameters) ;
Parameters.CholCov = chol((-Parameters.Hess)^-1);
X = mvnrnd(Parameters.ArgMax,Parameters.CholCov*Parameters.CholCov',1000);
Parameters.Dens = gmdistribution.fit(X,1);
Parameters.LogRatioFun = @LogRatioGMMrand;
Parameters.SampleFun = @SampleGMMrand;

AccRates = [];
RelESSs = [];
for i = 1:length(Epss)
    i
    Parameters.Epsil = Epss(i);
    
        Res = RunMCMC(Parameters.ArgMax,Parameters,3000);
        AccRates(i) = Res.AccRate;
        RelESSs(i,:) = Res.RelESS;
    
end
TempRes = struct();
TempRes.AccRates = AccRates;
TempRes.RelESSs = RelESSs;
ResGMM1_Fish = TempRes;

subplot(2,1,1)
plot(Epss,AccRates)
title('Acc Rates')
subplot(2,1,2)
plot(Epss,RelESSs)
title('ESS Rel')


% GMM1 - 1GM
X = [mvnrnd(mu1,sigma1,10000);mvnrnd(mu2,sigma2,10000)];
Parameters.Dens = gmdistribution.fit(X,1);
Parameters.LogRatioFun = @LogRatioGMMrand;
Parameters.SampleFun = @SampleGMMrand;

AccRates = [];
RelESSs = [];
for i = 1:length(Epss)
    i
    Parameters.Epsil = Epss(i);
    Res = RunMCMC(Parameters.ArgMax,Parameters,3000);
    AccRates(i) = Res.AccRate;
    RelESSs(i,:) = Res.RelESS;
end
TempRes = struct();
TempRes.AccRates = AccRates;
TempRes.RelESSs = RelESSs;
ResGMM1_1GM = TempRes;


subplot(2,1,1)
plot(Epss,AccRates)
title('Acc Rates')
subplot(2,1,2)
plot(Epss,RelESSs)
title('ESS Rel')

% GMM2 - RW
Parameters.Dens = Parameters.OptDens;
Parameters.LogRatioFun = @LogRatioGMMrand;
Parameters.SampleFun = @SampleGMMrand;

AccRates = [];
RelESSs = [];
for i = 1:length(Epss)
    i
    Parameters.Epsil = Epss(i);
    Res = RunMCMC(Parameters.ArgMax,Parameters,3000);
    AccRates(i) = Res.AccRate;
    RelESSs(i,:) = Res.RelESS;
end
TempRes = struct();
TempRes.AccRates = AccRates;
TempRes.RelESSs = RelESSs;
ResGMM2_RW = TempRes;


subplot(2,1,1)
plot(Epss,AccRates)
title('Acc Rates')
subplot(2,1,2)
plot(Epss,RelESSs)
title('ESS Rel')


% GMM2 - Lang
Parameters.Dens = Parameters.OptDens;
Parameters.LogRatioFun = @LogRatioGMMLang;
Parameters.SampleFun = @SampleGMMLang;

AccRates = [];
RelESSs = [];
for i = 1:length(Epss)
    i
    Parameters.Epsil = Epss(i);
    Res = RunMCMC(Parameters.ArgMax,Parameters,3000);
    AccRates(i) = Res.AccRate;
    RelESSs(i,:) = Res.RelESS;
end
TempRes = struct();
TempRes.AccRates = AccRates;
TempRes.RelESSs = RelESSs;
ResGMM2_Lang = TempRes;


subplot(2,1,1)
plot(Epss,AccRates)
title('Acc Rates')
subplot(2,1,2)
plot(Epss,RelESSs)
title('ESS Rel')



% MALA - Fish
Parameters.ScalingCov = (-Parameters.Hess)^-1;
Parameters.LogRatioFun = @LogRatioMALA;
Parameters.SampleFun = @SampleMALA;

AccRates = [];
RelESSs = [];
for i = 1:length(Epss)
    i
    Parameters.Epsil = Epss(i);
    Res = RunMCMC(Parameters.ArgMax,Parameters,3000);
    AccRates(i) = Res.AccRate;
    RelESSs(i,:) = Res.RelESS;
    Samples =  Res.Vals;
end
TempRes = struct();
TempRes.AccRates = AccRates;
TempRes.RelESSs = RelESSs;
TempRes.Samples = Samples;

ResMALA_Fish = TempRes;


subplot(2,1,1)
plot(Epss,AccRates)
title('Acc Rates')
subplot(2,1,2)
plot(Epss,RelESSs)
title('ESS Rel')


plot(ResMALA_Fish.Samples(1,:),ResMALA_Fish.Samples(2,:),'.')


