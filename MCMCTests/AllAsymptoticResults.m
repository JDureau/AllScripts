function Res = AllAsymptoticResults(Parameters,Epss)



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
    disp(['GMM1-Fish ' num2str(i)])
    Parameters.Epsil = Epss(i);
    Res = RunMCMC(Parameters.ArgMax,Parameters,3000);
    AccRates(i) = Res.AccRate;
    RelESSs(i,:) = Res.RelESS;
end
TempRes = struct();
TempRes.AccRates = AccRates;
TempRes.RelESSs = RelESSs;
ResGMM1_Fish = TempRes;

% subplot(2,1,1)
% plot(Epss,AccRates)
% title('Acc Rates')
% subplot(2,1,2)
% plot(Epss,RelESSs)
% title('ESS Rel')


% GMM1 - 1GM
X = [mvnrnd(mu1,sigma1,10000);mvnrnd(mu2,sigma2,10000)];
Parameters.Dens = gmdistribution.fit(X,1);
Parameters.LogRatioFun = @LogRatioGMMrand;
Parameters.SampleFun = @SampleGMMrand;

AccRates = [];
RelESSs = [];
for i = 1:length(Epss)
    disp(['GMM1-1GM ' num2str(i)])
    Parameters.Epsil = Epss(i);
    Res = RunMCMC(Parameters.ArgMax,Parameters,3000);
    AccRates(i) = Res.AccRate;
    RelESSs(i,:) = Res.RelESS;
end
TempRes = struct();
TempRes.AccRates = AccRates;
TempRes.RelESSs = RelESSs;
ResGMM1_1GM = TempRes;


% subplot(2,1,1)
% plot(Epss,AccRates)
% title('Acc Rates')
% subplot(2,1,2)
% plot(Epss,RelESSs)
% title('ESS Rel')

% GMMn - RW
Parameters.Dens = Parameters.OptDens;
Parameters.LogRatioFun = @LogRatioGMMrand;
Parameters.SampleFun = @SampleGMMrand;

AccRates = [];
RelESSs = [];
Samples = {};
for i = 1:length(Epss)
    disp(['GMMn-RW ' num2str(i)])
    Parameters.Epsil = Epss(i);
    Res = RunMCMC(Parameters.ArgMax,Parameters,1000);
    AccRates(i) = Res.AccRate;
    RelESSs(i,:) = Res.RelESS;
    Samples{i} =  Res.Vals;
end
TempRes = struct();
TempRes.AccRates = AccRates;
TempRes.RelESSs = RelESSs;
ResGMMn_RW = TempRes;


subplot(2,1,1)
plot(Epss,AccRates)
title('Acc Rates')
subplot(2,1,2)
plot(Epss,RelESSs)
title('ESS Rel')


% GMMn - Lang
Parameters.Dens = Parameters.OptDens;
Parameters.LogRatioFun = @LogRatioGMMLang;
Parameters.SampleFun = @SampleGMMLang;

AccRates = [];
RelESSs = [];
Samples = {};
for i = 1:length(Epss)
    disp(['GMMn-Lang ' num2str(i)])
    Parameters.Epsil = Epss(i);
    Res = RunMCMC(Parameters.ArgMax,Parameters,1000);
    AccRates(i) = Res.AccRate;
    RelESSs(i,:) = Res.RelESS;
    Samples{i} =  Res.Vals;
end
TempRes = struct();
TempRes.AccRates = AccRates;
TempRes.RelESSs = RelESSs;
ResGMMn_Lang = TempRes;


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
    disp(['MALA-Fish ' num2str(i)])
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


Res.ResMALA_Fish = ResMALA_Fish;
Res.ResGMMn_Lang = ResGMMn_Lang;
Res.ResGMMn_RW = ResGMMn_RW;
Res.ResGMM1_1GM = ResGMM1_1GM;
Res.ResGMM1_Fish = ResGMM1_Fish;






