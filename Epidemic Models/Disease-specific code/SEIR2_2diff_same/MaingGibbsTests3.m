
SEIRModel = struct();
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),coeff*Data.Observations(5,IndTime),coeff*Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);


clf
Parameters.ComputationTStep = 0.1;
betas =  (0.8 + cumsum(randn(1,20+1)*0.1*sqrt(Parameters.ComputationTStep)));
Parameters.betainit.Value = 0.8;% betas(1);
Parameters.betainit.TransfValue = log(0.8);

Data.Instants = [0:50]/Parameters.ComputationTStep;
Data.ObservedVariables = 1*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];

Parameters.ObsNoise =  0.1 ;
DataGen = SEIR_CreateData(betas,Parameters,Data,SEIRModel);
Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
plot(DataGen.)

Parameters.Betas = betas;
Parameters.NoPaths = 0;
Parameters.GibbsAdaptC = 0.99;
Parameters.MCMCType = 'Rand';
Parameters.ModelType = 'SMC';
Parameters.NbParticules = 500;
Parameters.GMeth = 'frfr';
Parameters.PathsToKeep = [1 2]';
Parameters.PMCMC = 'Gibbs';
Parameters.GMeth = 'cst given';
Parameters.G = 1;
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
Parameters.PMCMC = 'Gibbs';
Parameters.GMeth = 'frfr';
Parameters.IndBetaVar = 1;
Parameters.SigmaRW.GibbsSigma = 0.3;
ResGibbs = RunEstimationMethod(DataGen, SEIRModel,Parameters,TempPar,2000);





% DataGen = SEIR_CreateData(betas,Parameters,Data,SEIRModel);
% PathLikGen = SEIRModel.PathLik(DataGen,SEIRModel,Parameters);
ind = 1;
subplot(5,1,1)
plot(ResGibbs.Thetas(ind,:))
% hold on
% plot(ResGibbs.PropTransfThetas,'g')
% hold off
xlabel('Iterations')
ylabel('\sigma traceplot')
subplot(5,1,2)
plot(ResGibbs.LogLiks)
xlabel('Iterations')
ylabel('p(y|\sigma) traceplot')
subplot(5,1,3)
[F,XI]=KSDENSITY(ResGibbs.Thetas(ind,:)) ;
plot(XI,F)
hold on
plot(sqrt(mean(diff((DataGen.Observations(1,2:end))).^2/Parameters.ComputationTStep)),0,'og')
hold off
xlabel('\sigma')
ylabel('p(\sigma|y) posterior')
subplot(5,1,4)
plot(mean(squeeze(ResGibbs.Paths(:,1,:))))
hold on
plot(quantile(squeeze(ResGibbs.Paths(:,1,:)),0.025),'r')
plot(quantile(squeeze(ResGibbs.Paths(:,1,:)),0.975),'r')
plot((betas),'g')
xlabel('time')
ylabel('\beta_t')
hold off
subplot(5,1,5)
plot(ResGibbs.Coalescence/50)

