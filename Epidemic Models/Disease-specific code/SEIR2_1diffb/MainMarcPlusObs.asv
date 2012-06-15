% HPA data

DataPath = 'H:\My Documents\PhD Work\Data\HPA';

A = load([DataPath '\andre_estimates_31_01.txt']);

Data.Dates = {};
Data.NewCases = {};
for i = 1:size(A,1)
    Data.Dates{i} = i;
    for j = 1:7
        Data.NewCases{i}{j} = A(i,j);
    end
end

InitialDate = struct();
InitialDate.Month = 6;
InitialDate.Day = 1;
InitialDate.Year = 2009;
Data.Dates = ApproxWeeklyDates(InitialDate,35);



PopWeigths = [667600,2461800,5904100,6862500,14417400,12847800,7929300];
PopProps = PopWeigths/sum(PopWeigths)*100;
Weigthed = sum((A(:,1:7)*diag(PopWeigths.^-1)*diag(PopProps/100)*100000)');

plot(A*diag(PopWeigths.^-1)*100000)
hold on
plot(Weigthed,'k','LineWidth',2)
hold off
legend('<1','1-4','5-14','15-24','25-44','45-64','65+')




%% Inference
cd('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Filtering')
addpath('H:\My Documents\PhD Work\Matlab Scripts\General Tools')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes\Resampling\pf_resampling')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Joint Sampling')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\MIF')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Model Selection')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Optimization Approach')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\Models')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\SEIR')

Parameters = struct();

Parameters.Problem = 'MarcFluPlusObs';
Parameters.NbVariables = 7;
Parameters.SigmaObs = 0.1;
Parameters.DiffusionType = 'Add';
Parameters.ObservationLength = 7*35;
Parameters.ComputationTStep = 0.1;
Parameters.TotalPopulation = 100000;

Data.Observations = zeros(6,35);
Data.Observations(5,:) = Weigthed;
Data.Instants = [1:35]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];



Parameters.km1.Value = 2;
Parameters.km1.Min = -10^14;
Parameters.km1.Max = 10^14;
Parameters.km1.MinLim = 0.1;
Parameters.km1.MaxLim = 2.1;
Parameters.km1.Estimated = 1;
Parameters.km1.TransfType = 'Logit';
Parameters.gammam1.Value = 2;
Parameters.gammam1.Min = -10^14;
Parameters.gammam1.Max = 10^14;
Parameters.gammam1.MinLim = 0.4;
Parameters.gammam1.MaxLim = 2.1;
Parameters.gammam1.Estimated = 1;
Parameters.gammam1.TransfType = 'Logit';
Parameters.betainit.Value = 0.25;
Parameters.betainit.Min = -10^14;
Parameters.betainit.Max = 10^14;
Parameters.betainit.Estimated = 1;
Parameters.betainit.TransfType = 'Log';
Parameters.betainit.Init = 1;
Parameters.ObsPropinit.Value = 0.10;
Parameters.ObsPropinit.Min = 0.08;
Parameters.ObsPropinit.Max = 0.12;
Parameters.ObsPropinit.MinLim = 0;
Parameters.ObsPropinit.MaxLim = 1;
Parameters.ObsPropinit.Estimated = 1;
Parameters.ObsPropinit.TransfType = 'Logit';
Parameters.ObsPropinit.Init = 1;
Parameters.EInitProp.Value = max(1,10*Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Min = -10^14;%0.2*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Max = 10^14;%5*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.MinLim = 0;
Parameters.EInitProp.MaxLim = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.EInitProp.TransfType = 'Logit';
Parameters.EInitProp.Init = 1;
Parameters.IInitProp.Value = max(1,10*Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Min = -10^14;%0.2*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Max = 10^14;%5*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.MinLim = 0;
Parameters.IInitProp.MaxLim = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.IInitProp.TransfType = 'Logit';
Parameters.IInitProp.Init = 1;
Parameters.RInitProp.Value = 0.2;
Parameters.RInitProp.Min = 0;
Parameters.RInitProp.Max = 0.30;
Parameters.RInitProp.MinLim = 0;
Parameters.RInitProp.MaxLim = 0.6;
Parameters.RInitProp.Estimated = 1;
Parameters.RInitProp.TransfType = 'Logit';
Parameters.RInitProp.Init = 1;
Parameters.SigmaRW.Value = exp(-0.6);
Parameters.SigmaRW.Min = -10^14;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Log';
Parameters.SigmaRWObs.Value = exp(-5.6);
Parameters.SigmaRWObs.Min = -10^14;
Parameters.SigmaRWObs.Max = 10^14;
Parameters.SigmaRWObs.Estimated = 1;
Parameters.SigmaRWObs.TransfType = 'Log';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)

SavePath = 'S:\Results\';
Res = load([SavePath 'Marc_FirstEst_WithPaths.mat']);
ResPost = Res.Res3;
ParametersPost = ResPost.Parameters;
NamesPost = ParametersPost.Names.Estimated;
for i = 1:length(NamesPost)
    Parameters.(NamesPost{i}).Value = mean(ResPost.Thetas(i,:));
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

    

SEIRPlusObsModel = struct();
SEIRPlusObsModel.EKF_projection = @SEIRPlusObs_EKF_projection;
SEIRPlusObsModel.InitializeParameters = @SEIRPlusObsInitialize;
SEIRPlusObsModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs)';
SEIRPlusObsModel.SMC_projection = @SEIRPlusObs_SMC_projection;
Parameters = SEIRPlusObsModel.InitializeParameters(Parameters);

% EKF Optimization
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];



temp = zeros(1,7);
temp(1,5) = 1;
SEIRPlusObsModel.ObservationJacobian = {};
SEIRPlusObsModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.Instants)
    SEIRPlusObsModel.ObservationJacobian{i} = temp;
    SEIRPlusObsModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs*Data.Observations(5,i))^2;
end


Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.k.Estimated = 1;
Parameters.gamma.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.ObsPropinit.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRWObs.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRPlusObsModel,Parameters),Initialization,optimset('MaxIter',100,'TolX',1e-8,'TolFun',1e-5));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  

 Temp = EstimationEKFGen(Data, SEIRPlusObsModel, Parameters);


Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).SamplStd = 0.01*Parameters.(Names{i}).Value;
end
ResKal = KalmanNumericDerivativesWithPrior(Data,SEIRPlusObsModel,Parameters);

Test = mean(eig(-ResKal.Hess)>0)==1;
disp(Test)
Cov = (-ResKal.Hess)^-1;

%% SMC Optimization

Parameters.ComputationTStep = 0.1;

Data.Observations = zeros(6,35);
Data.Observations(5,:) = A(:,7);
Data.Instants = [1:35]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];


Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.ObsPropinit.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRWObs.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRPlusObsModel,Parameters),Initialization,optimset('MaxIter',100,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  


% defining NbParts 
nbparts = [500 800 1000 1300 1600];
NbTests = 50;
LogLikrecords = [];
for i = 1:length(nbparts)
    Parameters.NbParticules = nbparts(i);
    for j = 1:NbTests
        disp([i j])
        Temp = EstimationSMCsmoothGen(Data, SEIRPlusObsModel, Parameters);
        LogLikrecords(i,j) = Temp.LogLik;
    end
end
clf
hold on
for i = 1:length(nbparts)
    inds = find(not(isinf(LogLikrecords(i,:))));
    plot(nbparts(i),std(LogLikrecords(i,inds)),'o')
end
hold off


Parameters.NbParticules = 1000;


% defining TStep
DataTests = Data;

res = [1:10].^-1;
logliks = [];
paths = [];
for i = 1:length(res)
    disp(i)
    for j = 1:10
        Parameters.ComputationTStep = res(i);
        DataTests.Observations = zeros(6,35);
        DataTests.Observations(5,:) = A(:,7);
        DataTests.Instants = [1:35]*7/Parameters.ComputationTStep;
        DataTests.ObservedVariables = 5*ones(1,length(DataTests.Instants));
        DataTests.NbComputingSteps = [0 diff(DataTests.Instants)];
        Temp = EstimationSMCsmoothGen(DataTests, SEIRPlusObsModel, Parameters);
        logliks(i,j) = Temp.LogLik;
        paths(i,j,:) = exp(Temp.PosteriorMeansRecord(6,:));
    end
end
clf
hold on
for i = 1:length(res)
    plot(res(i),mean(logliks(i,:)),'.')
    plot(res(i),mean(logliks(i,:))+std(logliks(i,:)),'.r')
    plot(res(i),mean(logliks(i,:))-std(logliks(i,:)),'.r')
end
hold off

cols = rand(length(res),3);
clf
hold on
for i = 1:length(res)
    plot(mean(squeeze(paths(i,:,:))),'col',cols(i,:),'LineWidth',2)
    plot(quantile(squeeze(paths(i,:)),0.95),'col',cols(i,:))
    plot(quantile(squeeze(paths(i,:)),0.05),'col',cols(i,:))
end
legend
hold off


Parameters.ComputationTStep = 1/3;

Data.Observations = zeros(6,35);
Data.Observations(5,:) = Weigthed;
Data.Instants = [0:34]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];

Names = Parameters.Names.Estimated;     
Cov = zeros(length(Names));
for i = 1:length(Names)
    Cov(i,i) = 0.1*abs(Parameters.(Names{i}).TransfValue);
end
    

Parameters.ModelType = 'SMC';
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.1;
TempPar = ProposeInitialParameter(Data, SEIRPlusObsModel, Parameters);
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Res = RunEstimationMethod(Data, SEIRPlusObsModel,Parameters,TempPar,5000);

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
Res2 = RunEstimationMethod(Data, SEIRPlusObsModel,Parameters,TempPar,5000);

SavePath = 'S:\Results\';
save([SavePath 'MarcPlusObs_FirstEst_NoPaths.mat'],'Res2')


TempRes = Res2;

Cov = 2.38^2/dim*cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.6;
TempPar = TempRes.TempPar;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:7]';
Res3 = RunEstimationMethod(Data, SEIRPlusObsModel,Parameters,TempPar,10000);

SavePath = 'S:\Results\';
save([SavePath 'MarcPlusObs_FirstRes_WithPaths.mat'],'Res3')

Res = load([SavePath 'MarcPlusObs_FirstEst_WithPaths.mat'])











