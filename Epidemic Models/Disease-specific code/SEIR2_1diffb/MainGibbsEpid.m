
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SEIRModel = struct();

SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),coeff*Data.Observations(5,IndTime),coeff*Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
SEIRModel.PathLik = @SEIR_PathLik;

clf
Parameters.ComputationTStep = 1;
% betas =  exp(randn(1,1)*0.2+1.2+ cumsum(randn(1,30*7+1)*0.1*sqrt(Parameters.ComputationTStep)));
RealAmpl = 0.1;
RealBaseline = 0.4;
T = 30*7;
step = 1;
RealInflpt = 0.4*T/step;
RealSteepness = T/20;
xis = 0:step:T;
betas =  exp(RealBaseline-RealAmpl*Sigmoid((xis-RealInflpt*step)/RealSteepness));
Parameters.betainit.Value = betas(1);
Parameters.betainit.TransfValue = log(Parameters.betainit.Value);

Data.Instants = [0:30]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 1*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];

Parameters.ObsNoise = 0.05;
Parameters.SigmaObs.Value = Parameters.ObsNoise;
DataGen = SEIR_CreateData(betas,Parameters,Data,SEIRModel);

subplot(2,1,1)
plot(DataGen.Observations(5,:))
subplot(2,1,2)
plot(betas)

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.GibbsSigma = 0.045;
Parameters.km1.Estimated = 1;
Parameters.km1.GibbsSigma = 0.018;
Parameters.gammam1.Estimated = 1;
Parameters.gammam1.GibbsSigma = 0.0028;
Parameters.betainit.Estimated = 1;
Parameters.betainit.GibbsSigma = 0.0026;
Parameters.EInitProp.Estimated = 1;
Parameters.EInitProp.GibbsSigma = 0.12;
Parameters.IInitProp.Estimated = 1;
Parameters.IInitProp.GibbsSigma = 0.11;
Parameters.RInitProp.Estimated = 1;
Parameters.RInitProp.GibbsSigma = 0.015;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


Parameters.Betas = betas;
Parameters.NoPaths = 0;
Parameters.GibbsAdaptC = 0.995;
Parameters.MCMCType = 'Rand';
Parameters.ModelType = 'SMC';
Parameters.NbParticules = 1000;
Parameters.GMeth = 'frfr';
Parameters.PathsToKeep = [1:6]';
Parameters.PMCMC = 'Gibbs';
Parameters.GMeth = 'cst given';
Parameters.G = 1;
Parameters.SigmaRW.Value = 0.1;
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.Epsil = 1;
TempPar = ProposeInitialParameter(DataGen, SEIRModel, Parameters);
TempPar.SigmaRW.TransfValue = log(0.08);
Parameters.SigmaRW.Value = 0.08;
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.Model = 'SEIR';
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.IndBetaVar = 6;
Parameters.PMCMC = 'Gibbs';
Parameters.GMeth = 'cst given';
Parameters.G = 45^(-1);
Parameters.AdaptC = 0.999;
ResGibbs = RunEstimationMethod(DataGen, SEIRModel,Parameters,TempPar,10000);
Parameters.PMCMC = 'Gibbs';
Parameters.GMeth = 'frfr';
ResGibbs = RunEstimationMethod(DataGen, SEIRModel,Parameters,TempPar,20000);
% Parameters.PMCMC = 'GibbsRepar';
% Parameters.GMeth = 'frfr';
% Parameters.SigmaRW.GibbsSigma = 3;
% ResGibbsRepar = RunEstimationMethod(DataGen, SEIRModel,Parameters,TempPar,20000);
% 

Ress.PMMH = ResPMMH;
Ress.Gibbs = ResGibbs;
Ress.GibbsRepar = ResGibbsRepar;

SavePath = 'S:\Results\';
save([SavePath 'GibbsSEIRTests20NoiseVeryShort.mat'],'Ress')

% ResGibbs2  = SimpleGibbs(betas-betas(1),DataGen.Observations-betas(1),Parameters,100);

Names = ResGibbs.Parameters.Names.Estimated;
k = ceil(sqrt(length(Names)));
for i = 1:length(Names)
    subplot(k,k,i)
    plot(ResGibbs.Thetas(i,:))
    title(Names{i})
end


clf
subplot(3,1,1)
plot(ResGibbsRepar.TransfThetas)
subplot(3,1,2)
plot(ResGibbs.TransfThetas,'g')
subplot(3,1,3)
plot(ResPMMH.TransfThetas,'r')

plot(ResGibbs.TransfThetas,'g')
hold off
subplot(2,1,2)
plot(ResGibbsRepar.Paths(:,6,ceil(length(betas)/2)))
hold on
plot(ResGibbs.Paths(:,6,ceil(length(betas)/2)),'g')
plot(ResPMMH.Paths(:,6,ceil(length(betas)/2)),'r')
hold off

clf
subplot(5,2,1:2)
n = size(ResGibbs.Paths,3);
plot((0:n-1)/7,mean(squeeze(ResPMMH.Paths(:,6,:))))
hold on
plot((0:n-1)/7,quantile(squeeze(ResPMMH.Paths(:,6,:)),0.025),'r')
plot((0:n-1)/7,quantile(squeeze(ResPMMH.Paths(:,6,:)),0.975),'r')
plot((0:n-1)/7,log(betas),'g')
xlabel('time')
ylabel('\beta_t')
hold off
title('22 diffusion steps observed through 3 epidemic observations with high uncertainty (25/%)','FontWeight','bold')
subplot(5,2,3)
plot(autocorrelation(ResGibbsRepar.TransfThetas,300))
hold on
plot(autocorrelation(ResGibbs.TransfThetas,300),'g')
plot(autocorrelation(ResPMMH.TransfThetas,300),'r')
hold off
title('\sigma: samples autocorrelation','FontWeight','bold')
subplot(5,2,4)
[F,XI]=KSDENSITY(ResGibbsRepar.TransfThetas(ind,:)) ;
plot(XI,(F),'b')
hold on
[F,XI]=KSDENSITY(ResGibbs.TransfThetas(ind,:)) ;
plot(XI,(F),'g')
[F,XI]=KSDENSITY(ResPMMH.TransfThetas(ind,:)) ;
plot(XI,(F),'r')
hold off
title('\sigma: posterior distribution estimate','FontWeight','bold')
subplot(5,2,5)
plot(autocorrelation(ResGibbsRepar.Paths(:,6,ceil(size(ResGibbs.Paths,3)/3)),30))
hold on
plot(autocorrelation(ResGibbs.Paths(:,6,ceil(size(ResGibbs.Paths,3)/3)),30),'g')
plot(autocorrelation(ResPMMH.Paths(:,6,ceil(size(ResGibbs.Paths,3)/3)),30),'r')
hold off
title('\beta(T/3): samples autocorrelation','FontWeight','bold')
subplot(5,2,6)
[F,XI]=KSDENSITY(ResGibbsRepar.Paths(:,6,ceil(size(ResGibbs.Paths,3)/3))) ;
plot(XI,F,'b')
hold on
[F,XI]=KSDENSITY(ResGibbs.Paths(:,6,ceil(size(ResGibbs.Paths,3)/3))) ;
plot(XI,F,'g')
[F,XI]=KSDENSITY(ResPMMH.Paths(:,6,ceil(size(ResGibbs.Paths,3)/3))) ;
plot(XI,F,'r')
hold off
title('\beta(T/3): posterior distribution estimate','FontWeight','bold')
subplot(5,2,7)
plot(autocorrelation(ResGibbsRepar.Paths(:,6,ceil(size(ResGibbs.Paths,3)*2/3)),30))
hold on
plot(autocorrelation(ResGibbs.Paths(:,6,ceil(size(ResGibbs.Paths,3)*2/3)),30),'g')
plot(autocorrelation(ResPMMH.Paths(:,6,ceil(size(ResGibbs.Paths,3)*2/3)),30),'r')
hold off
title('\beta(2T/3): samples autocorrelation','FontWeight','bold')
subplot(5,2,8)
[F,XI]=KSDENSITY(ResGibbsRepar.Paths(:,6,ceil(size(ResGibbs.Paths,3)*2/3))) ;
plot(XI,F,'b')
hold on
[F,XI]=KSDENSITY(ResGibbs.Paths(:,6,ceil(size(ResGibbs.Paths,3)*2/3))) ;
plot(XI,F,'g')
[F,XI]=KSDENSITY(ResPMMH.Paths(:,6,ceil(size(ResGibbs.Paths,3)*2/3))) ;
plot(XI,F,'r')
hold off
title('\beta(2T/3): posterior distribution estimate','FontWeight','bold')
subplot(5,2,9)
plot(autocorrelation(ResGibbsRepar.Paths(:,6,ceil(size(ResGibbs.Paths,3)*3/3)),30))
hold on
plot(autocorrelation(ResGibbs.Paths(:,6,ceil(size(ResGibbs.Paths,3)*3/3)),30),'g')
plot(autocorrelation(ResPMMH.Paths(:,6,ceil(size(ResGibbs.Paths,3)*3/3)),30),'r')
hold off
title('\beta(T): samples autocorrelation','FontWeight','bold')
subplot(5,2,10)
[F,XI]=KSDENSITY(ResPMMH.Paths(:,6,ceil(size(ResGibbs.Paths,3)*3/3))) ;
plot(XI,F,'r')
hold on
[F,XI]=KSDENSITY(ResGibbs.Paths(:,6,ceil(size(ResGibbs.Paths,3)*3/3))) ;
plot(XI,F,'g')
[F,XI]=KSDENSITY(ResGibbsRepar.Paths(:,6,ceil(size(ResGibbs.Paths,3)*3/3))) ;
plot(XI,F,'b')
hold off
title('\beta(T): posterior distribution estimate','FontWeight','bold')
legend('PMMH','Gibbs','Gibbs with repar')

plot(ResGibbs.Thetas)
% hold on
% plot(ResGibbs2.SigmasRecord,'g')
% plot(SigmasRecord,'r')
% hold off

clf
[F,XI]=KSDENSITY(log(gamrnd(0.001,0.001^(-1),100,1).^(-1))) ;
plot(XI,F)


plot(ResGibbs.Thetas(ind,:),ResGibbs.Paths(:,6,ceil(length(betas)/2)),'.')

empirsigs = [];
Res = ResGibbsRepar;
for i = 1:2000
    Beta = Res.Paths(i,6,:);
    empirsigs(i) = sqrt(mean(diff(Beta).^2));
end
plot(empirsigs,Res.Thetas,'.')
hold on
plot(min(empirsigs):0.001:max(empirsigs),min(empirsigs):0.001:max(empirsigs),'g')
hold off

% DataGen = SEIR_CreateData(betas,Parameters,Data,SEIRModel);
% PathLikGen = SEIRModel.PathLik(DataGen,SEIRModel,Parameters);
ind = 1;
Res = ResPMMH;
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
plot(mean(squeeze(Res.Paths(:,6,:))))
hold on
plot(quantile(squeeze(Res.Paths(:,6,:)),0.025),'r')
plot(quantile(squeeze(Res.Paths(:,6,:)),0.975),'r')
plot(log(betas),'g')
xlabel('time')
ylabel('\beta_t')
hold off
subplot(5,1,5)
plot(Res.Coalescence/Res.Parameters.NbParticules)


