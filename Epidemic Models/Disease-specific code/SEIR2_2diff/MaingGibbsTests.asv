% Main Gibbs tests


%% load everything, set up dataset


cd('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Disease-specific code\SEIR')
addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\General Tools')
addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Generic PMCMC tools')
addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Disease-specific code\DirectObs')
addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Toolboxes')


DataPath = 'H:\My Documents\PhD Work\Data\HPA';

A = load([DataPath '\andre_estimates_31_01.txt']);

Data.Dates = {};
Data.NewCases = {};
for i = 1:size(A,1)
    Data.Dates{i} = i;
    for j = 1:7
    Data.NewCases{i}{j} = A(i,j)*10;
    end
end

InitialDate = struct();
InitialDate.Month = 6;
InitialDate.Day = 1;
InitialDate.Year = 2009;
Data.Dates = ApproxWeeklyDates(InitialDate,35);

plot(A)
legend('0-4','5-14','15-24','25-44','45-64','65+')

PopWeigths = [667600,2461800,5904100,6862500,14417400,12847800,7929300];
PopProps = PopWeigths/sum(PopWeigths)*100;
Weigthed = sum((A(:,1:7)*diag(PopWeigths.^-1)*diag(PopProps/100)*100000)');

plot(A*diag(PopWeigths.^-1)*100000)
hold on
plot(Weigthed,'k','LineWidth',2)
hold off
legend('<1','1-4','5-14','15-24','25-44','45-64','65+')


SavePath = 'S:\Results\';
tmp = load([SavePath 'ParametersSEIR_MarcAdd.mat']);
Parameters = tmp.Parameters;

TStep = Parameters.ComputationTStep;


Data.Observations = zeros(7,35);
Data.Observations(5,:) = Weigthed*10;

Difftype = 'Add';



Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:7]';
Parameters.TotalPopulation = 100000;





% FullSEIRinference(Data,'Add','Fixed','test')



SavePath = 'S:\Results\';
tmp = load([SavePath 'ParametersSEIR.mat']);
Parameters = tmp.Parameters;

Data.Instants = [1:size(Data.Observations,2)]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];

ObsType = 'Fixed';
if strcmp(ObsType,'Fixed')
    Parameters.SigmaObs.Value = 0.1;
elseif strcmp(ObsType,'Estimated')
    Parameters.SigmaObs.Value = 0.1;
    Parameters.SigmaObs.Min = 0.05;
    Parameters.SigmaObs.Max = 0.15;
    Parameters.SigmaObs.MinLim = 0.01;
    Parameters.SigmaObs.MaxLim = 0.19;
    Parameters.SigmaObs.Estimated = 1;
    Parameters.SigmaObs.TransfType = 'Log';
    Parameters.SigmaObs.Sample = 1;
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

Parameters.DiffusionType = Difftype;
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';



Parameters.EInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% Parameters.EInitProp.Min = 0;%0.0000000002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% Parameters.EInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
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
% SEIRModel.LikFunction = 'normpdf(log(Variables(:,5)),log(coeff*Data.Observations(5,IndTime)),sqrt(log(1+coeff*Parameters.SigmaObs.Value)))';
SEIRModel.LikFunction = 'normpdf(Variables(:,5),coeff*Data.Observations(5,IndTime),coeff*Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
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
Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
% Parameters.betainit.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;

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
Parameters.betainit.Estimated = 1;
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
Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.betainit.Estimated = 1;
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
Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);


Test = 0;

NbIts = 0;

while not(Test)
    
    Parameters = KalmOpt(Parameters,Data,SEIRModel,15000);

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
Parameters = KalmOpt(Parameters,Data,SEIRModel,15000);
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
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizewithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',50,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  




SavePath = 'S:\Results\ResultsBiostats\';
Res = load([SavePath 'MarcData_Add.mat']);
Res = Res.Res3;

EmpirCov = cov(Res.TransfThetas');


Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.betainit.Estimated = 1;
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


Parameters.PMCMC = 'Gibbs';
SEIRModel.PathLik = @SEIR_PathLik;

if strcmp(Parameters.PMCMC,'Gibbs')
%     tmpCov = (-KalHess)^-1;
    for i = 1:length(Names)
        ind = Parameters.(Names{i}).Index;
        Parameters.(Names{i}).GibbsSigma = sqrt(2.38^2/7*EmpirCov(ind,ind));
        Parameters.(Names{i}).GibbsSampler = @GibbsRWsampler;
    end
end
if strcmp(Parameters.DiffusionType,'IBM')
    Parameters.ComputationTStep = 1/3;
    Parameters.NbParticules = 10000;
else
    Parameters.NbParticules = 2000;
    Parameters.ComputationTStep = 1/3;
end
Data.Instants = [0:size(Data.Observations,2)-1]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];


Parameters.NbParticules = 2000;
Parameters.ModelType = 'SMC';
dim = length(Parameters.Names.Estimated);
Cov = 2.38^2/dim*EmpirCov;
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
if strcmp(Parameters.PMCMC,'Gibbs')
    Parameters.NoPaths = 0;
else
    Parameters.NoPaths = 1;
end
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);

Parameters.MCMCType =  'Rand';
Parameters.PMCMC = 'Gibbs';
Parameters.GibbsRepar = 0;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
ResGibbs = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,5000);
SavePath = 'S:\Results\';
save([SavePath 'MarcData_Add_GibbsTest_NoRepar.mat'],'ResGibbs');


Parameters.NoPaths = 0;
Parameters.PMCMC = 'Gibbs';
Parameters.GibbsRepar = 1;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
ResPMMH = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,10);
SavePath = 'S:\Results\';
save([SavePath 'MarcData_Add_GibbsTest.mat'],'ResPMMH');




load([SavePath 'MarcData_Add_GibbsTest_NoRepar.mat'])
load([SavePath 'MarcData_Add_PMMHTest.mat'])
figure(1)
PlotMarc(ResGibbs,8)
figure(2)
PlotMarc(ResPMMH,8)

clf
Res = ResGibbs
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
   subplot(4,2,i)
%    plot(autocorrelation(ResPMMH.Thetas(i,:),100),'k')
%    hold on
%    plot(autocorrelation(ResGibbs.Thetas(i,:),100),'g')
%    hold off
   plot(Res.TransfThetas(i,:));%,Res.Paths(:,6,400),'.')
   title([Names{i} ' (Acc=' num2str(ResGibbs.GibbsAccRates(i)) ')'])
end

%% Add SigmaObs





Parameters.SigmaObs.Estimated = 1;
Parameters.SigmaObs.GibbsSigma = 0.5;
Parameters.SigmaObs.GibbsSampler = @GibbsRWsampler;
Parameters.SigmaObs.Init = 1;
Parameters.SigmaObs.TransfType = 'Log';
Parameters.SigmaObs.Min = -10^14;
Parameters.SigmaObs.Max = 10^14;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);




% PMMH
tmp = zeros(8,8);
tmp(2:8,2:8) = 2.38^2/dim*EmpirCov;
tmp(1,1) = 2.38^2*0.2^2;
Parameters.G = tmp^-1;
Parameters.Epsil = 0.8;
% Parameters.G = 0.8^-2;
Parameters.PMCMC = 'PMMH';
Parameters.NoPaths = 1;
Parameters.Adapt = 0;
Parameters.AdaptC = 0.99;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
ResPMMH = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,5000);


NewEmpirCov = 2.38^2/8*cov(ResPMMH.TransfThetas');
Parameters.Epsil = 0.8;
Parameters.G = NewEmpirCov^-1;
Parameters.PMCMC = 'PMMH';
Parameters.NoPaths = 0;
Parameters.Adapt = 0;
Parameters.AdaptC = 0.99;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
ResPMMH = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,1000);


SavePath = 'S:\Results\';
Res = load([SavePath 'MarcData_Add.mat']);
save([SavePath 'MarcData_Add_PMMHTest_WithSigmaObs.mat'],'ResPMMH');

% Gibbs


tmpCov = NewEmpirCov;
for i = 1:length(Names)
    ind = Parameters.(Names{i}).Index;
    Parameters.(Names{i}).GibbsSigma = sqrt(2.38^2*NewEmpirCov(ind,ind));
    Parameters.(Names{i}).GibbsSampler = @GibbsRWsampler;
end


SEIRModel.PathLik = @SEIR_PathLik;
Parameters.PMCMC = 'Gibbs';
Parameters.NoPaths = 0;
Parameters.MCMCType = 'Rand';
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
ResGibbs = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,1000);


SavePath = 'S:\Results\';
Res = load([SavePath 'MarcData_Add.mat']);
save([SavePath 'MarcData_Add_GibbsTest_WithSigmaObs.mat'],'ResGibbs');


load([SavePath 'MarcData_Add_PMMHTest_WithSigmaObs.mat'])
load([SavePath 'MarcData_Add_GibbsTest_WithSigmaObs.mat'])
figure(1)
PlotMarc(ResGibbs,8)
figure(2)
PlotMarc(ResPMMH,8)

clf
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
   subplot(4,2,i)
%    plot(autocorrelation(ResPMMH.Thetas(i,:),100),'k')
%    hold on
%    plot(autocorrelation(ResGibbs.Thetas(i,:),100),'g')
%    hold off
   [N,X] = hist(ResPMMH.Thetas(i,:));
   hist(ResPMMH.Thetas(i,:));
   hold on
   try
       yis = normpdf(X,(Parameters.(Names{i}).Max+Parameters.(Names{i}).Min)/2,(Parameters.(Names{i}).Max-Parameters.(Names{i}).Min)/4)
       plot(X,yis/mean(yis)*mean(N),'g')
   end
   title(Names{i})
end
hold off

Res = Res.Res3;
clf
Names = Res.Parameters.Names.Estimated;
ind = 0;
for i = 1:length(Names)-1
    for j = i+1:length(Names)
       ind = ind+1;
       subplot(5,5,ind)
%    plot(autocorrelation(ResPMMH.Thetas(i,:),100),'k')
%    hold on
%    plot(autocorrelation(ResGibbs.Thetas(i,:),100),'g')
%    hold off
       plot(Res.TransfThetas(i,:),Res.TransfThetas(j,:),'.')
       title([Names{i} ' ' Names{j}])
    end
end




%% Add SigmaObs and rho




Parameters.MultCoeff.Value = 10;
Parameters.MultCoeff.Estimated = 1;
Parameters.MultCoeff.Init = 0;
Parameters.MultCoeff.TransfType = 'Log';
Parameters.MultCoeff.Min = -10^14;
Parameters.MultCoeff.Max = 10^14;
Parameters.SigmaObs.Estimated = 1;
Parameters.SigmaObs.GibbsSigma = 0.5;
Parameters.SigmaObs.GibbsSampler = @GibbsRWsampler;
Parameters.SigmaObs.Init = 0;
Parameters.SigmaObs.TransfType = 'Log';
Parameters.SigmaObs.Min = -10^14;
Parameters.SigmaObs.Max = 10^14;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);




% PMMH
tmp = zeros(9,9);
tmp(2:8,2:8) = 2.38^2/9*EmpirCov;
tmp(1,1) = 2.38^2*0.2^2;
tmp(9,9) = 2.38^2*0.2^2;
Parameters.G = tmp^-1;
Parameters.Epsil = 0.8;
% Parameters.G = 0.8^-2;
Parameters.PMCMC = 'PMMH';
Parameters.MCMCType = 'Rand'
Parameters.NoPaths = 1;
Parameters.Adapt = 0;
Parameters.AdaptC = 0.99;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
ResPMMH = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,500);


NewEmpirCov = 2.38^2/8*cov(ResPMMH.TransfThetas');
Parameters.Epsil = 0.8;
Parameters.G = NewEmpirCov^-1;
Parameters.PMCMC = 'PMMH';
Parameters.NoPaths = 0;
Parameters.Adapt = 0;
Parameters.AdaptC = 0.99;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
ResPMMH = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,1000);


SavePath = 'S:\Results\';
Res = load([SavePath 'MarcData_Add.mat']);
save([SavePath 'MarcData_Add_PMMHTest_WithSigmaObs.mat'],'ResPMMH');

% Gibbs


tmpCov = NewEmpirCov;
for i = 1:length(Names)
    ind = Parameters.(Names{i}).Index;
    Parameters.(Names{i}).GibbsSigma = sqrt(2.38^2*NewEmpirCov(ind,ind));
    Parameters.(Names{i}).GibbsSampler = @GibbsRWsampler;
end


SEIRModel.PathLik = @SEIR_PathLik;
Parameters.PMCMC = 'Gibbs';
Parameters.NoPaths = 0;
Parameters.MCMCType = 'Rand';
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
ResGibbs = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,1000);


SavePath = 'S:\Results\';
Res = load([SavePath 'MarcData_Add.mat']);
save([SavePath 'MarcData_Add_GibbsTest_WithSigmaObs.mat'],'ResGibbs');


load([SavePath 'MarcData_Add_PMMHTest_WithSigmaObs.mat'])
load([SavePath 'MarcData_Add_GibbsTest_WithSigmaObs.mat'])
figure(1)
PlotMarc(ResGibbs,8)
figure(2)
PlotMarc(ResPMMH,8)


Names = Parameters.Names.Estimated;
for i = 1:length(Names)
   subplot(3,3,i)
%    plot(autocorrelation(ResPMMH.Thetas(i,:),100),'k')
%    hold on
%    plot(autocorrelation(ResGibbs.Thetas(i,:),100),'g')
%    hold off
   plot(ResPMMH.Thetas(i,:))
   title(Names{i})
end


Res = Res.Res3;
clf
Names = Res.Parameters.Names.Estimated;
ind = 0;
for i = 1:length(Names)-1
    for j = i+1:length(Names)
       ind = ind+1;
       subplot(5,5,ind)
%    plot(autocorrelation(ResPMMH.Thetas(i,:),100),'k')
%    hold on
%    plot(autocorrelation(ResGibbs.Thetas(i,:),100),'g')
%    hold off
       plot(Res.TransfThetas(i,:),Res.TransfThetas(j,:),'.')
       title([Names{i} ' ' Names{j}])
    end
end


%% On sims

Parameters.ComputationTStep = 0.1;
Parameters.ObsLength = 10;
Parameters.Sigma.RealValue = 0.1;
Parameters.SigmaObs = 0.05;

n = Parameters.ObsLength/Parameters.ComputationTStep;
Data.Beta = cumsum(randn(1,n)*Parameters.Sigma.RealValue*sqrt(Parameters.ComputationTStep));
Data.NbComputingSteps = [0 ones(1,Parameters.ObsLength)/Parameters.ComputationTStep];
Data.Observations = [0 Data.Beta(cumsum(ones(1,Parameters.ObsLength)/Parameters.ComputationTStep))+Parameters.SigmaObs*randn(1,Parameters.ObsLength)];



NbIts = 2000;
NbParts = 500;

BetasRecord = zeros(NbIts,sum(Data.NbComputingSteps)+1);
SigmasRecord = Parameters.Sigma.RealValue;
Parameters.Sigma.Value = 0.9;%Parameters.Sigma.RealValue;
Parameters.Sigma.TransfValue = log(Parameters.Sigma.Value);
Accepted = [];
lambda = 0.0172;
C = 0.995;
% beta given sigma
BetasSamples = zeros(NbParts,sum(Data.NbComputingSteps));
ind = 1;
for IndStep = 2:Parameters.ObsLength+1
    for IndCStep = 1:1/Parameters.ComputationTStep
        BetasSamples(:,ind+1) = BetasSamples(:,ind)+randn(NbParts,1)*Parameters.Sigma.Value*sqrt(Parameters.ComputationTStep);
        ind = ind+1;
    end
    Weigths = normpdf(BetasSamples(:,ind),Data.Observations(IndStep),Parameters.SigmaObs);
    Weigths = Weigths/sum(Weigths);
    u = rand(1,1)/NbParts;
    s = 0;
    KeptInds = [];
    resind = 1;
    for ipart = 1:NbParts
        k = 0;
        s = s+Weigths(ipart);
        while s>u
            k=k+1;
            u = u+1/NbParts;
            KeptInds(resind) = ipart;
            resind = resind+1;
        end
    end
    BetasSamples = BetasSamples(KeptInds,:);
end
plot(Data.Beta)
hold on
plot(cumsum(Data.NbComputingSteps),Data.Observations,'g')
RandInd = ceil(rand(1,1)*NbParts);
Beta = BetasSamples(RandInd,:);
plot(Beta,'r')
hold off
for IndIt = 2:NbIts
    % Sigma given Beta
    RandInd = ceil(rand(1,1)*NbParts);
    Beta = BetasSamples(RandInd,:);
    BetasRecord(IndIt,:) = Beta;
    Parameters.SigmaStar.TransfValue = Parameters.Sigma.TransfValue + lambda*randn(1,1);
    Parameters.SigmaStar.Value = exp(Parameters.SigmaStar.TransfValue);
    BetaStar = Beta/Parameters.Sigma.Value*Parameters.SigmaStar.Value;
    LogLik = sum(log(normpdf(Beta(cumsum(Data.NbComputingSteps(2:end))+1),Data.Observations(2:end),Parameters.SigmaObs)));
    LogLikStar = sum(log(normpdf(BetaStar(cumsum(Data.NbComputingSteps(2:end))+1),Data.Observations(2:end),Parameters.SigmaObs)));
    LogCorr = log(exp(-Parameters.Sigma.TransfValue));
    LogCorrStar = log(exp(-Parameters.SigmaStar.TransfValue));
    AccRate = LogLikStar-LogCorrStar-LogLik+LogCorr;
    if log(rand(1,1))<AccRate
        Parameters.Sigma.Value = Parameters.SigmaStar.Value;
        Parameters.Sigma.TransfValue = Parameters.SigmaStar.TransfValue;
        Beta = BetaStar;
        Accepted(IndIt-1) = 1;
    else
        Accepted(IndIt-1) = 0;
    end
    if IndIt>10
        lambda = exp(log(lambda)-C^IndIt*(0.23-mean(Accepted)));
    end
    SigmasRecord(IndIt) = Parameters.Sigma.Value;
    
    % beta given sigma
    BetasSamples = zeros(NbParts,sum(Data.NbComputingSteps));
    ind = 1;
    for IndStep = 2:Parameters.ObsLength+1
        for IndCStep = 1:1/Parameters.ComputationTStep
            BetasSamples(:,ind+1) = BetasSamples(:,ind)+randn(NbParts,1)*Parameters.Sigma.Value*sqrt(Parameters.ComputationTStep);
            BetasSamples(1,ind+1) = Beta(ind+1);
            ind = ind+1;
        end
        Weigths = normpdf(BetasSamples(:,ind),Data.Observations(IndStep),Parameters.SigmaObs);
        Weigths = Weigths/sum(Weigths);
        u = rand(1,1)/NbParts;
        s = 0;
        KeptInds = [];
        resind = 1;
        for ipart = 1:NbParts
            k = 0;
            s = s+Weigths(ipart);
            while s>u
                k=k+1;
                u = u+1/NbParts;
                KeptInds(resind) = ipart;
                resind = resind+1;
            end
        end
        KeptInds(1) = 1;  
        BetasSamples = BetasSamples(KeptInds,:);
    end
    disp([num2str(IndIt) '  ' num2str(mean(Accepted)*100,3) '%   ' num2str(lambda,3) ])
end


subplot(3,1,1)
plot(SigmasRecord)
subplot(3,1,2)

plot(Data.Beta)
hold on
plot(cumsum((Data.NbComputingSteps)),Data.Observations,'g')
plot(mean(BetasRecord),'k')
plot(quantile(BetasRecord,0.025),'r')
plot(quantile(BetasRecord,0.975),'r')
hold off

subplot(3,1,3)
hist(SigmasRecord)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DObs = struct();

DObs.InitializeParameters = @DObsInitialize;
DObs.LikFunction = 'normpdf(Data.Observations(IndTime),Variables(:,1),Parameters.ObsNoise)';
DObs.SMC_projection = @DObs_SMC_projection;
DObs.PathLik = @DObs_PathLik;

clf
Parameters.ComputationTStep = 1;
betas =  (0 + cumsum(randn(1,20+1)*0.1*sqrt(Parameters.ComputationTStep)));
% betas = 1+0.1*sqrt(Parameters.ComputationTStep)*(1:50);
Parameters.betainit.Value = 0;% betas(1);
Parameters.betainit.TransfValue = 0;

Data.Instants = [0:20]/Parameters.ComputationTStep;
Data.ObservedVariables = 1*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];

Parameters.ObsNoise =  0.01 ;
DataGen = DObs_CreateData(betas,Parameters,Data,DObs);
% DataGen = Data;
% DataGen.Observations = [0 Data.Observations];


Parameters.TrueB = (log(betas)-Parameters.betainit.TransfValue)/Parameters.SigmaRW.Value;
Parameters.SigmaRW.GibbsSampler = @GibbsRWsampler;

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
% Parameters.betainit.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

% Parameters.SigmaRW.Value = 0.1;
% Parameters = UpdateParsNoTransfToTransf(Parameters);
% eps = [1,0.5,0.25,0.2,0.1,0.05,0.02];
% LogLiks = [];
% for i = 1:length(eps)
%     Data.Instants = [1:20]/eps(i);
%     Data.ObservedVariables = 5*ones(1,length(Data.Instants));
%     Data.NbComputingSteps = [0 diff(Data.Instants)];
%     for j = 1:10
%         [i,j]
%         Res = EstimationSMCsmoothGen(DataGen,DObs,Parameters);
%         LogLiks(i,j) = Res.LogLik;
%     end
% end

Parameters.Betas = betas;
Parameters.NoPaths = 0;
Parameters.GibbsAdaptC = 0.99;

Parameters.MCMCType = 'Rand';
Parameters.ModelType = 'SMC';
Parameters.NbParticules = 500;
Parameters.GMeth = 'frfr';
Parameters.PathsToKeep = [1 2]';
Parameters.betainit.GibbsSigma = 0.5;
Parameters.Epsil = 0;
Parameters.SigmaRW.Value = 0.1;
Parameters.PMCMC = 'Gibbs';
Parameters.PMCMC = 'PMMH';
Parameters.GMeth = 'cst given';
Parameters.G = 1;
Parameters = UpdateParsNoTransfToTransf(Parameters);
TempPar = ProposeInitialParameter(DataGen, DObs, Parameters);
clf
plot(betas)
hold on
plot(DataGen.Instants, DataGen.Observations(1,:),'g')
plot(TempPar.Paths(1,:),'r')
hold off
TempPar.SigmaRW.TransfValue = log(0.1);
Parameters.SigmaRW.Value = 0.1;
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.GibbsReparType = 'classic';
Parameters.Model = 'DObs';
Parameters.SigmaRW.Min = -10^14;
Parameters.SigmaRW.Max = 10^14;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
% Parameters.PMCMC = 'PMMH';
% Parameters.GMeth = 'cst given';
Parameters.PMCMC = 'Gibbs';
Parameters.GMeth = 'frfr';
Parameters.G = 1;
Parameters.Epsil = 1;
Parameters.IndBetaVar = 1;
Parameters.SigmaRW.GibbsSigma = 0.1;
TempPar.Paths = [Beta;Beta];
DataGen.Observations = Data.Observations;
ResGibbs = RunEstimationMethod(DataGen, DObs,Parameters,TempPar,1000);





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
% subplot(4,1,3)
% plot(ResGibbs.MargLogLiksStars-ResGibbs.MargLogLiks)
% xlabel('Iterations')
% ylabel('p(y|B,\sigma)')
subplot(5,1,3)
[F,XI]=KSDENSITY(ResGibbs.Thetas(ind,:)) ;
plot(XI,F)
hold on
plot(sqrt(mean(diff((DataGen.Observations(1,2:end))).^2/Parameters.ComputationTStep)),0,'og')
hold off
xlabel('\sigma')
ylabel('p(\sigma|y) posterior')
subplot(5,1,4)
% plot(mean(squeeze(ResGibbs.Paths(:,6,:))))
% hold on
% plot(quantile(squeeze(ResGibbs.Paths(:,6,:)),0.025),'r')
% plot(quantile(squeeze(ResGibbs.Paths(:,6,:)),0.975),'r')
% plot(log(betas),'g')
% xlabel('time')
% hold off
plot(mean(squeeze(ResGibbs.Paths(:,1,:))))
hold on
plot(quantile(squeeze(ResGibbs.Paths(:,1,:)),0.025),'r')
plot(quantile(squeeze(ResGibbs.Paths(:,1,:)),0.975),'r')
plot((betas),'g')
xlabel('time')
ylabel('\beta_t')
hold off
subplot(5,1,5)
plot(ResGibbs.Coalescence/500)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


Parameters.ObsLength = 20;
Parameters.Sigma.RealValue = 0.1;
Parameters.SigmaObs = 0.1;

n =  Parameters.ObsLength +1;
Data.Beta = [0 cumsum(randn(1,n-1)*Parameters.Sigma.RealValue)];
Data.Observations = [Data.Beta+Parameters.SigmaObs*randn(1,n)];


NbIts = 1000;
NbParts = 500;

BetasRecord = zeros(NbIts,n);
SigmasRecord = Parameters.Sigma.RealValue;
Parameters.Sigma.Value = Parameters.Sigma.RealValue;

Accepted = [];
lambda = 0.35
C = 0.995;
% beta given sigma
BetasSamples = zeros(NbParts,n);
for IndStep = 2:n
    BetasSamples(:,IndStep) = BetasSamples(:,IndStep-1)+randn(NbParts,1)*Parameters.Sigma.Value;
    Weigths = normpdf(BetasSamples(:,IndStep),Data.Observations(IndStep),Parameters.SigmaObs);
    Weigths = Weigths/sum(Weigths);
    u = rand(1,1)/NbParts;
    s = 0;
    KeptInds = [];
    resind = 1;
    for ipart = 1:NbParts
        k = 0;
        s = s+Weigths(ipart);
        while s>u
            k=k+1;
            u = u+1/NbParts;
            KeptInds(resind) = ipart;
            resind = resind+1;
        end
    end
    BetasSamples = BetasSamples(KeptInds,:);
end
plot(Data.Beta)
hold on
plot(Data.Observations,'g')
RandInd = ceil(rand(1,1)*NbParts);
Beta = BetasSamples(RandInd,:);
plot(Beta,'r')
hold off
alpha = 0.0001;
beta = 0.0001;
Accepted = [];
C = 0.99;

for IndIt = 2:NbIts
    disp(IndIt)
    % Sigma given Beta
    RandInd = ceil(rand(1,1)*NbParts);
    Beta = BetasSamples(RandInd,:);
    BetasRecord(IndIt,:) = Beta;
    
%     oldvalue = Parameters.Sigma.Value;
%     Parameters.Sigma.Value = sqrt(gamrnd((n-1)/2+alpha,(beta+sum(diff(Beta).^2)/2)^-1)^-1);
%     Beta = (Beta/oldvalue)*Parameters.SigmaStar.Value;
    
    Parameters.SigmaStar.Value = Parameters.Sigma.Value + randn(1,1)*lambda;
    LogLikStar = sum(log(normpdf(Data.Observations,(Beta/Parameters.Sigma.Value)*Parameters.SigmaStar.Value,Parameters.SigmaObs)));
    LogLik = sum(log(normpdf(Data.Observations,Beta,Parameters.SigmaObs)));
    AccRate = LogLikStar-LogLik;
    if and(Parameters.SigmaStar.Value>0,log(rand(1,1))<AccRate)
        Beta = (Beta/Parameters.Sigma.Value)*Parameters.SigmaStar.Value;
        Parameters.Sigma.Value = Parameters.SigmaStar.Value;
        Accepted(end+1) = 1;
    else
        Accepted(end+1) = 0;
    end
    if IndIt>100
        lambda = exp(log(lambda)-C^IndIt*(0.23-mean(Accepted)));
    end
    [mean(Accepted) lambda]
    
    SigmasRecord(IndIt) = Parameters.Sigma.Value;
    
    % beta given sigma
    BetasSamples = zeros(NbParts,n);
    for IndStep = 2:n
        BetasSamples(:,IndStep) = BetasSamples(:,IndStep-1)+randn(NbParts,1)*Parameters.Sigma.Value;
        BetasSamples(1,IndStep) = Beta(IndStep);
        Weigths = normpdf(BetasSamples(:,IndStep),Data.Observations(IndStep),Parameters.SigmaObs);
        Weigths = Weigths/sum(Weigths);
        u = rand(1,1)/NbParts;
        s = 0;
        KeptInds = [];
        resind = 1;
        for ipart = 1:NbParts
            k = 0;
            s = s+Weigths(ipart);
            while s>u
                k=k+1;
                u = u+1/NbParts;
                KeptInds(resind) = ipart;
                resind = resind+1;
            end
        end
        KeptInds(1) = 1; 
        BetasSamples = BetasSamples(KeptInds,:);
    end
end
mean(SigmasRecord.^2)
var(SigmasRecord.^2)

subplot(3,1,1)
plot(SigmasRecord)
subplot(3,1,2)
plot(Data.Beta)
hold on
plot(Data.Observations,'g')
plot(mean(BetasRecord),'k')
plot(quantile(BetasRecord,0.025),'r')
plot(quantile(BetasRecord,0.975),'r')
hold off

subplot(3,1,3)
plot(autocorrelation(SigmasRecord,300))



