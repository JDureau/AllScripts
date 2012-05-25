function [] = FullSEIRinference_knownpars(Data,DiffType,ObsType,Name)


SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc';
tmp = load([SavePath '/ParametersSEIR.mat']);
Parameters = tmp.Parameters;

Data.Instants = [1:size(Data.Observations,2)]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];

if strcmp(ObsType,'Fixed')
    Parameters.SigmaObs.Value = 0.1;
elseif strcmp(ObsType,'Estimated')
    Parameters.SigmaObs.Value = sqrt(log(1.02));
    Parameters.SigmaObs.Min = -10^14;
    Parameters.SigmaObs.Max =  10^14;
%     Parameters.SigmaObs.MinLim = sqrt(log(1.0001));
%     Parameters.SigmaObs.MaxLim = sqrt(log(1.16));
    Parameters.SigmaObs.Estimated = 1;
    Parameters.SigmaObs.TransfType = 'Log';
    Parameters.SigmaObs.Sample = 1;
    Parameters.SigmaObs.Init = 0;
elseif strcmp(ObsType,'CoeffStudy')
    Parameters.SigmaObs.Value = 0.1;
    Parameters.MultCoeff.Value = Data.Coeff;
    Parameters.MultCoeff.Min =  -10^14;
    Parameters.MultCoeff.Max =   10^14;
    Parameters.MultCoeff.MinLim = 5;
    Parameters.MultCoeff.MaxLim = 50;
    Parameters.MultCoeff.Estimated = 0;
    Parameters.MultCoeff.TransfType = 'Logit';
    Parameters.MultCoeff.Sample = 1;
end
Parameters.Correction = 1;


Parameters.betainit.TransfValue = 0.7;
Parameters = UpdateParsTransfToNoTransf(Parameters);


Parameters.DiffusionType = DiffType;
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';

% 
% 
% Parameters.EInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% % Parameters.EInitProp.Min = 0;%0.0000000002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% % Parameters.EInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% Parameters.IInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% Parameters.IInitProp.Min = 0;%0.0000000002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% Parameters.IInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)

    

SEIRModel = struct();
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
disp('BE CAREFUL! Lik function changed for Gibbs')
SEIRModel.LikFunction = 'normpdf(log(Variables(:,5)),log(coeff*Data.Observations(5,IndTime)),Parameters.SigmaObs.Value)';
% SEIRModel.LikFunction = 'normpdf(log(Variables(:,5)),log(coeff*Data.Observations(5,IndTime))-0.5*Parameters.SigmaObs.Value^2,Parameters.SigmaObs.Value)';
% SEIRModel.LikFunction = 'normpdf(Variables(:,5),coeff*Data.Observations(5,IndTime),coeff*Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);

% temp = zeros(1,7);
% temp(1,5) = 1;
% SEIRModel.ObservationJacobian = {};
% SEIRModel.ObservationMeasurementNoise = {};
% for i = 1:length(Data.Instants)
%     SEIRModel.ObservationJacobian{i} = temp;
%     SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs.Value*Data.Observations(5,i))^2;
% end


Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
% Parameters.km1.Estimated = 1;
% Parameters.gammam1.Estimated = 1;
% Parameters.EInitProp.Estimated = 1;
% Parameters.IInitProp.Estimated = 1;
% Parameters.RInitProp.Estimated = 1;
% Parameters.betainit.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
if strcmp(ObsType,'Estimated')
    Parameters.SigmaObs.Estimated = 1;
end

Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Test = 0;
try
    Temp = EstimationEKFGen(Data, SEIRModel, Parameters);    
    if (Temp.LogLik>-600)
        Test = 1;
    end
end
NbIts = 0;
while not(Test)
    Parameters = SampleParameters(Parameters);
    Parameters.SigmaRW.Value = rand(1,1)*2;
    Parameters = UpdateParsNoTransfToTransf(Parameters);
    Parameters.InitialCov = rand(1,1)*0.4;

    try
        Temp = EstimationEKFGen(Data, SEIRModel, Parameters);
        if (Temp.LogLik>-400)
            Test = 1;
        end
    end
    NbIts = NbIts + 1;
    if NbIts>20000
        disp('Can''t initialize IBM')
        die
    end
end

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.SigmaRW.Estimated = 1;
% Parameters.betainit.Estimated = 1;
% Parameters.EInitPropNoise.Estimated = 1;
% Parameters.IInitPropNoise.Estimated = 1;
% Parameters.RInitPropNoise.Estimated = 1;
if strcmp(Parameters.DiffusionType,'IBM')
    Parameters.betaderinit.Estimated = 1;
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);



Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
% Parameters.km1.Estimated = 1;
% Parameters.gammam1.Estimated = 1;
% Parameters.EInitProp.Estimated = 1;
% Parameters.IInitProp.Estimated = 1;
% Parameters.RInitProp.Estimated = 1;
% Parameters.betainit.Estimated = 1;
% Parameters.EInitPropNoise.Estimated = 1;
% Parameters.IInitPropNoise.Estimated = 1;
% Parameters.RInitPropNoise.Estimated = 1;
% Parameters.betainitNoise.Estimated = 1;
% Parameters.InitialCov = Parameters.InitialCov/2;
if strcmp(ObsType,'Estimated')
    Parameters.SigmaObs.Estimated = 1 ;
elseif strcmp(ObsType,'CoeffEst')
    Parameters.MultCoeff.Estimated = 1 ;
end
if strcmp(Parameters.DiffusionType,'IBM')
    Parameters.betaderinit.Estimated = 0;
    %     Parameters.betaderinitNoise.Estimated = 1;
end
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);


% Parameters.InitialCov = 0;
% Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);


Test = 0;

NbIts = 0;

while not(Test)
    
    Parameters = KalmOpt(Parameters,Data,SEIRModel,1500);

    try
        KalHess = Parameters.KalHess;
        Test = Parameters.KalmMaxReached;
    end
    
    NbIts = NbIts + 1;
    disp(NbIts)
    if NbIts>50
        return
    end
end


Cov = (-KalHess)^-1;

Parameters.Correction = 0;
Parameters = KalmOpt(Parameters,Data,SEIRModel,1500);
Parameters.Correction = 1;


%% SMC Optimization

disp('Phase 2')



Parameters.NoPaths = 1;
Parameters.PathsToKeep = [1:7]';
if strcmp(Parameters.DiffusionType,'IBM')
    Parameters.NbParticules = 8000;
else
    Parameters.NbParticules = 4000;
end
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;

Initialization = [];
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
Parameters.Correction = 1;
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizewithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',20,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  

% 
% % defining NbParts 
% nbparts = [  100 1000 2000 3000 4000  5000 ];
% % nbparts = [  5000 6000 7000 8000 10000];
% NbTests = 50;
% LogLikrecords = [];
% for i = 1:length(nbparts)
%     Parameters.NbParticules = nbparts(i);
%     for j = 1:NbTests
%         disp([i j])
%         Temp = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);
%         LogLikrecords(i,j) = Temp.LogLik;
%     end
% end
% clf
% hold on
% for i = 1:length(nbparts)
%     plot(nbparts(i),std(LogLikrecords(i,:)),'o')
% end
% hold off
% clf
% hold on
% for i = 1:length(nbparts)
%     plot(nbparts(i),mean(LogLikrecords(i,:)),'.')
%     plot(nbparts(i),mean(LogLikrecords(i,:))+std(LogLikrecords(i,:)),'.r')
%     plot(nbparts(i),mean(LogLikrecords(i,:))-std(LogLikrecords(i,:)),'.r')
% end
% hold off
% 
% 
% Parameters.NbParticules = 4000;
% 
% % defining TStep
% DataTests = Data;
% 
% res = [2 3 10].^-1;
% logliks = [];
% paths = [];
% for i = 1:length(res)
%     disp(i)
%     for j = 1:20
%         Parameters.ComputationTStep = res(i);
%         DataTests.Instants = [1:size(Data.Observations,2)]*7/Parameters.ComputationTStep;
%         DataTests.ObservedVariables = 5*ones(1,length(DataTests.Instants));
%         DataTests.NbComputingSteps = [0 diff(DataTests.Instants)];
%         Temp = EstimationSMCsmoothGen(DataTests, SEIRModel, Parameters);
%         logliks(i,j) = Temp.LogLik;
%         paths(i,j,:) = exp(Temp.PosteriorMeansRecord(6,:));
%     end
% end
% logliks'
% clf
% hold on
% for i = 1:length(res)
%     plot(res(i),mean(logliks(i,:)),'.')
%     plot(res(i),mean(logliks(i,:))+std(logliks(i,:)),'.r')
%     plot(res(i),mean(logliks(i,:))-std(logliks(i,:)),'.r')
% end
% hold off
% 
% cols = rand(length(res),3);
% clf
% hold on
% for i = 1:length(res)
%     plot(mean(squeeze(paths(i,:,:))),'col',cols(i,:),'LineWidth',2)
%     plot(quantile(squeeze(paths(i,:)),0.95),'col',cols(i,:))
%     plot(quantile(squeeze(paths(i,:)),0.05),'col',cols(i,:))
% end
% legend
% hold off


% if strcmp(Parameters.PMCMC,'Gibbs')
%     tmpCov = (-KalHess)^-1;
%     for i = 1:length(Names)
%         ind = Parameters.(Names{i}).Index;
%         Parameters.(Names{i}).GibbsSigma = sqrt(2.38^2*tmpCov(ind,ind));
%         Parameters.(Names{i}).GibbsSampler = @GibbsRWsampler;
%     end
% end



if strcmp(Parameters.DiffusionType,'IBM')
    Parameters.ComputationTStep = 1/3;
    Parameters.NbParticules = 10000;
else
    Parameters.NbParticules = 1000;
    Parameters.ComputationTStep = 1/3;
end
    

Data.Instants = [0:size(Data.Observations,2)-1]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];


Parameters.ModelType = 'SMC';
dim = length(Parameters.Names.Estimated);
Cov = 2.38^2/dim*(-KalHess)^-1;
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
% if strcmp(Parameters.PMCMC,'Gibbs')
%     Parameters.NoPaths = 0;
% else
%     Parameters.NoPaths = 1;
% end
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
Parameters.AdaptC = 0.99;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,1000);

TempRes = Res;

dim = length(Parameters.Names.Estimated);
Cov = 2.38^2/dim*cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
TempPar = TempRes.TempPar;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Res2 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,1000);

SavePath = 'S:\Results\';
save([Name '.mat'],'Res2')


TempRes = Res2;

Cov = 2.38^2/dim*cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
TempPar = TempRes.TempPar;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);

% if strcmp(Parameters.DiffusionType,'IBM')
%     Parameters.NoPaths = 0;
% else
%     Parameters.NoPaths = 0;
% end
% Parameters.PathsToKeep = [1:7]';
Res3 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,10000);
Res3.Description = 'This one has unif prior on days for latent and inf periods. It also uses the actual averaged data and not only >65 as FirstEst';

save(Name,'Res3')

'uniform priors?'

