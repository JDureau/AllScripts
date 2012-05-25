
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DObs = struct();

DObs.InitializeParameters = @DObsInitialize;
DObs.LikFunction = 'normpdf(Variables(:,1),Data.Observations(IndTime),Parameters.SigmaObs.Value)';
DObs.SMC_projection = @DObs_SMC_projection;
DObs.PathLik = @DObs_PathLik;

clf
Parameters.ComputationTStep = 1;
betas =  (0.8 + cumsum(randn(1,40+1)*0.1*sqrt(Parameters.ComputationTStep)));
Parameters.betainit.Value = 0.8;% betas(1);
Parameters.betainit.TransfValue = log(0.8);

Data.Instants = [0:40]/Parameters.ComputationTStep;
Data.ObservedVariables = 1*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];

Parameters.ObsNoise =  0.15 ;
Parameters.SigmaObs.Value = Parameters.ObsNoise;
DataGen = DObs_CreateData(betas,Parameters,Data,DObs);
Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaObs.Estimated = 1;
Parameters.SigmaObs.Value = 0.15;
Parameters.SigmaObs.Min = -10^14;
Parameters.SigmaObs.Max = 10^14;
Parameters.SigmaObs.MinLim = 0;
Parameters.SigmaObs.MaxLim = 50;
Parameters.SigmaObs.Estimated = 1;
Parameters.SigmaObs.TransfType = 'Log';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


Parameters.Betas = betas;
Parameters.NoPaths = 0;
Parameters.GibbsAdaptC = 0.99;
Parameters.MCMCType = 'Rand';
Parameters.ModelType = 'SMC';
Parameters.NbParticules = 1000;
Parameters.GMeth = 'frfr';
Parameters.PathsToKeep = [1 2]';
Parameters.PMCMC = 'Gibbs';
Parameters.GMeth = 'cst given';
Parameters.G = 0.05^(-1)*eye(2);
Parameters.SigmaRW.Value = 0.1;
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.Epsil = 1;
TempPar = ProposeInitialParameter(DataGen, DObs, Parameters);
TempPar.SigmaRW.TransfValue = log(0.1);
Parameters.SigmaRW.Value = 0.1;
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.GibbsReparType = 'classic';
Parameters.Model = 'DObs';
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.IndBetaVar = 1;
Parameters.PMCMC = 'PMMH';
Parameters.GMeth = 'cst given';
Parameters.AdaptC = 0.99;
ResPMMH = RunEstimationMethod(DataGen, DObs,Parameters,TempPar,20000);
Parameters.PMCMC = 'Gibbs';
Parameters.GMeth = 'frfr';
ResGibbs = RunEstimationMethod(DataGen, DObs,Parameters,TempPar,2000);
Parameters.PMCMC = 'GibbsRepar';
Parameters.GMeth = 'frfr';
Parameters.SigmaRW.GibbsSigma = 0.000001;
ResGibbsRepar = RunEstimationMethod(DataGen, DObs,Parameters,TempPar,2000);
Parameters.Lambda = 0.009;
ResGibbsRepar2  = SimpleGibbs(betas-betas(1),DataGen.Observations-betas(1),Parameters,2000);

clf
subplot(2,1,1)
plot(autocorrelation(ResGibbsRepar2.SigmasRecord,300),'g')
hold on
plot(autocorrelation(ResGibbsRepar.Thetas,300))
hold off
subplot(2,1,2)
plot(ResGibbsRepar2.SigmasRecord,'g')
hold on
plot(ResGibbsRepar.Thetas)
hold off



clf
subplot(2,1,1)
plot(ResGibbsRepar.Thetas)
% hold on
% plot(ResGibbs.Thetas,'g')
% plot(ResPMMH.Thetas,'r')
% hold off
subplot(2,1,2)
plot(ResGibbsRepar.Paths(:,1,ceil(length(betas)/2)))
% hold on
% plot(ResGibbs.Paths(:,1,ceil(length(betas)/2)),'g')
% plot(ResPMMH.Paths(:,1,ceil(length(betas)/2)),'r')
hold off


plot(ResGibbs.Thetas(ind,:),ResGibbs.Paths(:,1,ceil(length(betas)/2)),'.')



clf
subplot(2,2,1)
plot(autocorrelation(ResGibbsRepar.Thetas,300))
hold on
plot(autocorrelation(ResGibbs.Thetas,300),'g')
plot(autocorrelation(ResPMMH.Thetas,300),'r')
hold off
subplot(2,2,2)
[F,XI]=KSDENSITY(ResGibbsRepar.Thetas(ind,:)) ;
plot(XI,F,'b')
hold on
[F,XI]=KSDENSITY(ResGibbs.Thetas(ind,:)) ;
plot(XI,F,'g')
[F,XI]=KSDENSITY(ResPMMH.Thetas(ind,:)) ;
plot(XI,F,'r')
hold off
subplot(2,2,3)
plot(autocorrelation(ResGibbsRepar.Paths(:,1,ceil(length(betas)/2)),30))
hold on
plot(autocorrelation(ResGibbs.Paths(:,1,ceil(length(betas)/2)),30),'g')
plot(autocorrelation(ResPMMH.Paths(:,1,ceil(length(betas)/2)),30),'r')
hold off
subplot(2,2,4)
[F,XI]=KSDENSITY(ResGibbsRepar.Paths(:,1,ceil(length(betas)/2))) ;
plot(XI,F,'b')
hold on
[F,XI]=KSDENSITY(ResGibbs.Paths(:,1,ceil(length(betas)/2))) ;
plot(XI,F,'g')
[F,XI]=KSDENSITY(ResPMMH.Paths(:,1,ceil(length(betas)/2))) ;
plot(XI,F,'r')
hold off




% DataGen = SEIR_CreateData(betas,Parameters,Data,SEIRModel);
% PathLikGen = SEIRModel.PathLik(DataGen,SEIRModel,Parameters);
ind = 1;
Res = ResGibbs;
subplot(5,1,1)
plot(Res.Thetas(ind,:))
% hold on
% plot(ResGibbs.PropTransfThetas,'g')
% hold off
xlabel('Iterations')
ylabel('\sigma traceplot')
subplot(5,1,2)
plot(Res.LogLiks)
xlabel('Iterations')
ylabel('p(y|\sigma) traceplot')
subplot(5,1,3)
[F,XI]=KSDENSITY(Res.Thetas(ind,:)) ;
plot(XI,F)
xlabel('\sigma')
ylabel('p(\sigma|y) posterior')
subplot(5,1,4)
plot(mean(squeeze(Res.Paths(:,1,:))))
hold on
plot(quantile(squeeze(Res.Paths(:,1,:)),0.025),'r')
plot(quantile(squeeze(Res.Paths(:,1,:)),0.975),'r')
plot((betas),'g')
xlabel('time')
ylabel('\beta_t')
hold off
subplot(5,1,5)
plot(Res.Coalescence/Res.Parameters.NbParticules)


subplot(4,2,1)
plot(ResPMMH.Thetas(1,:),ResPMMH.Paths(:,1,20),'.')
xlabel('\sigma')
ylabel('\beta(T/2)')
subplot(4,2,2)
plot(ResPMMH.Thetas(2,:),ResPMMH.Paths(:,1,20),'.')
xlabel('noise')
ylabel('\beta(T/2)')
subplot(4,2,3)
plot(ResPMMH.Thetas(1,:),ResPMMH.Thetas(2,:),'.')
xlabel('\sigma')
ylabel('noise')
subplot(4,2,4)
ESSs = [];
for i = 2:length(Data.Instants)
    ESSs(i) = 100/(1+2*sum(autocorrelation(ResPMMH.Paths(:,1,Data.Instants(i)),100)));
end
plot(ESSs(2:end))
xlabel('time')
ylabel('\beta ESS')
subplot(4,2,5)
plot(ResPMMH.Thetas(1,:))
tmp = 100/(1+2*sum(autocorrelation(ResPMMH.Thetas(1,:),100)));
title(['\sigma ESS= ' num2str(tmp,3) '%'])
xlabel('iterations')
ylabel('\sigma')
subplot(4,2,7)
[F,XI]=KSDENSITY(ResPMMH.Thetas(1,:)) ;
plot(XI,F)
hold on
plot(0.1, 0,'og')
hold off
xlabel('\sigma')
subplot(4,2,6)
plot(ResPMMH.Thetas(2,:))
tmp = 100/(1+2*sum(autocorrelation(ResPMMH.Thetas(2,:),100)));
title(['noise ESS= ' num2str(tmp,3) '%'])
xlabel('iterations')
ylabel('noise')
subplot(4,2,8)
[F,XI]=KSDENSITY(ResPMMH.Thetas(2,:)) ;
plot(XI,F)
hold on
plot(0.1, 0,'og')
hold off
xlabel('noise')


