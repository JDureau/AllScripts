% Google Flu Trends


%% Load Data
DataPath = 'H:\PhD Work\Data\GoogleFluTrends';

%France
fid = fopen([DataPath '\France.txt'],'r');
A = textscan(fid,'%s');
fclose(fid);
A = A{:};

regions = {'France','Alsace','Aquitaine','Auvergne','Burgundy','Brittany','Centre','Champagne-Ardenne','Franche-Comte','Ile-de-France','Languedoc-Roussillon','Lorraine','Midi-Pyrenees','Nord-Pas-de-Calais','Normandy - Lower','Normandy - Upper','Pays de la Loire','Picardie','Poitou-Charentes','Provence-Alpes-Cote d Azur','Rhône-Alpes'};
DataFr.Dates = {};
DataFr.NewCases = {};
for i = 2:length(A)
    Data.Dates{i-1} = A{i}(1:10);
    inds = [regexp(A{i}, ',') length(A{i})+1];
    for j = 1:length(inds)-1
        if inds(j)+1<inds(j+1)
            DataFr.NewCases{i-1}(j) = str2num(A{i}(inds(j)+1:inds(j+1)-1));
        else
            DataFr.NewCases{i-1}(j) = 0;
        end
    end
end
DataFr.regions = regions

%US
fid = fopen([DataPath '\US.txt'],'r');
A = textscan(fid,'%s');
fclose(fid);
A = A{:};

regions = {'United States','Alabama','Alaska','Arizona','Arkansas','California'};
DataUS.Dates = {};
DataUS.NewCases = {};
for i = 2:length(A)
    DataUS.Dates{i-1} = A{i}(1:10);
    inds = [regexp(A{i}, ',') length(A{i})+1];
    for j = 1:length(inds)-1
        if inds(j)+1<inds(j+1)
            DataUS.NewCases{i-1}(j) = str2num(A{i}(inds(j)+1:inds(j+1)-1));
        else
            DataUS.NewCases{i-1}(j) = 0;
        end
    end
end
DataUS.regions = regions;

ToPlot = {'United States'};
temp = [];
for j = 1:length(ToPlot)
    ind = find(strcmp(regions,ToPlot{j}));
    for i = 1:length(Data.NewCases)
        temp(j,i) = Data.NewCases{i}(ind);
    end
end
plot(temp')
legend(ToPlot)

%% Inference
cd('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Filtering')
addpath('H:\PhD Work\Matlab Scripts\General Tools')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Parameter Estimation')
addpath('H:\PhD Work\Matlab Scripts\Toolboxes\Resampling\pf_resampling')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Joint Sampling')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\MIF')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Model Selection')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Optimization Approach')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Generic Parameter Estimation')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Generic Parameter Estimation\Models')
addpath('H:\PhD Work\Matlab Scripts\Toolboxes')

Parameters = struct();

Parameters.NbVariables = 6;
Parameters.SigmaObs = 0.1;
Parameters.DiffusionType = 'Add';
Parameters.ObservationLength = 7*50;
% Parameters.ObservationLength = 7*3*4;
Parameters.ComputationTStep = 0.1;
Parameters.TotalPopulation = 1;


Data = ExtractGoogleTSeries(DataFr,Parameters,'France',2009);
plot(Data.Observations(5,:))

Parameters.k.Value = 1/3.5;
Parameters.k.Min = 1/4;
Parameters.k.Max = 1/3;
Parameters.k.Estimated = 1;
Parameters.k.TransfType = 'Log';
Parameters.gamma.Value = 1/7.5;
Parameters.gamma.Min = 1/8;
Parameters.gamma.Max = 1/3;
Parameters.gamma.Estimated = 1;
Parameters.gamma.TransfType = 'Log';
Parameters.betainit.Value = 0.25;
Parameters.betainit.Min = -10;
Parameters.betainit.Max = 10;
Parameters.betainit.Estimated = 1;
Parameters.betainit.TransfType = 'Log';
Parameters.EInitProp.Value = Data.Observations(5,1)/Parameters.TotalPopulation;
Parameters.EInitProp.Min = 0.5*Data.Observations(5,1)/Parameters.TotalPopulation;
Parameters.EInitProp.Max = 1.5*Data.Observations(5,1)/Parameters.TotalPopulation;
Parameters.EInitProp.Estimated = 1;
Parameters.EInitProp.TransfType = 'Log';
Parameters.IInitProp.Value = Data.Observations(5,1)/Parameters.TotalPopulation;
Parameters.IInitProp.Min = 0.5*Data.Observations(5,1)/Parameters.TotalPopulation;
Parameters.IInitProp.Max = 1.5*Data.Observations(5,1)/Parameters.TotalPopulation;
Parameters.IInitProp.Estimated = 1;
Parameters.IInitProp.TransfType = 'Log';
Parameters.RInitProp.Value = 0.5;
Parameters.RInitProp.Min = 0;
Parameters.RInitProp.Max = 0.9;
Parameters.RInitProp.Estimated = 1;
Parameters.RInitProp.TransfType = 'Logit';
Parameters.SigmaRW.Value = exp(-0.6);
Parameters.SigmaRW.Min = 0;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.Estimated = 0;
Parameters.SigmaRW.TransfType = 'Log';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)



SEIRModel = struct();
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);





temp = zeros(1,6);
temp(1,5) = 1;
SEIRModel.ObservationJacobian = {};
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.Instants)
    SEIRModel.ObservationJacobian{i} = temp;
    SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs*Data.Observations(5,i))^2;
end

% EKF Optimization
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];


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
Parameters.SigmaRW.Estimated = 1;
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
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  





% SMC Optimization
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];

NbTrials = 20;
Pars = {};
LogLiks = {};

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.k.Estimated = 1;
Parameters.gamma.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
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
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  










ToEstimate = {'SigmaRW','RInitProp','k','gamma'};
Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
for j = 1:length(ToEstimate)
    Parameters.(ToEstimate{j}).Estimated = 1;
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Parameters.MCMCType = 'Lang';
Parameters.GMeth = 'cst given';
Parameters.aim = 0.23;
Parameters.Epsil = 6;





% Compute Hessian.

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters) 

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).SamplStd = 0.01*Parameters.(Names{i}).Value;
end
Res = KalmanNumericDerivativesWithPrior(Data,SEIRModel,Parameters);

Test = mean(eig(-Res.Hess)>0)==1;
disp(Test)

Parameters.ComputationTStep = 0.5;
Data.Instants = [1:52]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];
Parameters.NoPaths = 1;

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',2000,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    ParametersC.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  


Cov = (-Res.Hess)^-1;
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.5;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
[Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
ResGoog = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,30000);


TempPar = ResGoog2.TempPar;
TransfThetasSamples = ResGoog2.TransfThetas;
Parameters = DefineScalingPars(TransfThetasSamples,Parameters);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,Parameters);
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'GMM';
Parameters.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',Parameters,1);
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);

Parameters.Epsil = 1;
[Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Parameters.PathsToKeep = [1:6]';
Parameters.NoPaths = 0;
ResGoog2 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,30000);


SavePath = 'C:\Users\dureau\Desktop\Results';
save([SavePath '\ResGoog_Alsace_RW.mat'],'ResGoog2')








%% Inference On Simulated data
cd('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Filtering')
addpath('H:\PhD Work\Matlab Scripts\General Tools')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Parameter Estimation')
addpath('H:\PhD Work\Matlab Scripts\Toolboxes\Resampling\pf_resampling')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Joint Sampling')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\MIF')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Model Selection')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Optimization Approach')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Generic Parameter Estimation')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Generic Parameter Estimation\Models')
addpath('H:\PhD Work\Matlab Scripts\Toolboxes')

Parameters = ResGoog2.Parameters;

Data = ResGoog2.Data;
ResultSMC = EstimationSMCsmoothGen(Data, SEIRModel, Parameters)

t = 0:1:Data.Instants(end)/Parameters.ComputationTStep;
plot(mean(squeeze(ResGoog2.Paths(:,6,:))))

FtPath = 0.6 - 1/6*cos((t+200)*(2*pi)/1300); 
plot(mean(squeeze(exp(ResGoog2.Paths(:,6,:)))))
hold on
plot(t*Parameters.ComputationTStep,FtPath)
hold off

record = [];
Parameters = SEIRInitialize(Parameters);
mpred = Parameters.InitialState;
TStep = Parameters.ComputationTStep;
ind = 0;
Obs = [];
RtPath = [];
for IndTime = 1:length(Data.Instants)
    mpred(5) = 0;
    for IndDiscr = 1:Data.NbComputingSteps(IndTime)
        TotPop = Parameters.TotalPopulation;
        mtemp = mpred;
        ind = ind+1;
        mtemp(6) = FtPath(ind);
        mpred(1) = mpred(1) + (-mtemp(6)*mtemp(1)*mtemp(3)/TotPop)*TStep;
        mpred(2) = mpred(2) + ( mtemp(6)*mtemp(1)*mtemp(3)/TotPop- Parameters.k.Value*mtemp(2))*TStep;
        mpred(3) = mpred(3) + ( Parameters.k.Value*mtemp(2) - Parameters.gamma.Value*mtemp(3))*TStep;
        mpred(4) = mpred(4) + ( Parameters.gamma.Value*mtemp(3))*TStep;
        mpred(5) = mpred(5) + ( Parameters.k.Value*mtemp(2))*TStep;
        mpred(6) = mtemp(6);
        record(:,end+1) = mpred;
        RtPath(end+1) = mpred(1)/TotPop*mpred(6)/Parameters.gamma.Value;
    end
    Obs(end+1) = mpred(5);
end
plot(record(5,:))


Data.Observations(5,:) = Obs;
ResultSMC = EstimationSMCsmoothGen(Data, SEIRModel, Parameters)
PlotResGoogle(ResultSMC,8)



SEIRModel = struct();
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);



temp = zeros(1,6);
temp(1,5) = 1;
SEIRModel.ObservationJacobian = {};
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.Instants)
    SEIRModel.ObservationJacobian{i} = temp;
    SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs*Data.Observations(5,i))^2;
end

% EKF Optimization
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];


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
Parameters.SigmaRW.Estimated = 1;
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
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  





% SMC Optimization
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];

NbTrials = 20;
Pars = {};
LogLiks = {};

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.k.Estimated = 1;
Parameters.gamma.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
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
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-8));
% [x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  


ToEstimate = {'SigmaRW','RInitProp','k','gamma'};
Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
for j = 1:length(ToEstimate)
    Parameters.(ToEstimate{j}).Estimated = 1;
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Parameters.MCMCType = 'Lang';
Parameters.GMeth = 'cst given';
Parameters.aim = 0.23;
Parameters.Epsil = 6;



% Compute Hessian.

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters) 

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).SamplStd = 0.01*Parameters.(Names{i}).Value;
end
Res = KalmanNumericDerivativesWithPrior(Data,SEIRModel,Parameters);

Test = mean(eig(-Res.Hess)>0)==1;
disp(Test)

Parameters.ComputationTStep = 0.5;
Data.Instants = [1:50]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];
Parameters.MCMCType = ',ll';
Parameters.GMeth = 'cst given';

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:6]';
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',2000,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    ParametersC.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  


Parameters.NoPaths = 1;
Parameters.PathsToKeep = [1:6]';
Cov = (-Res.Hess)^-1;
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.5;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
[Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
ResGoog = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,30000);


TempPar = ResGoog2.TempPar;
TransfThetasSamples = ResGoog2.TransfThetas;
Parameters = DefineScalingPars(TransfThetasSamples,Parameters);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,Parameters);
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'GMM';
Parameters.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',Parameters,1);
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);

Parameters.Epsil = 1;
[Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Parameters.PathsToKeep = [1:6]';
Parameters.NoPaths = 0;
ResGoog2 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,30000);


SavePath = 'C:\Users\dureau\Desktop\Results';
save([SavePath '\ResGoog_Alsace_RW_Sim.mat'],'ResGoog2')


%%%%%%%%%%%%%%%%%







ToEstimate = {'SigmaRW','RInitProp','k','gamma','betainit','EInitProp','IInitProp'};
% ToEstimate = {'SigmaRW','RInitProp'};
Parameters.NbParticules = 1000;
for i = 1:length(ToEstimate)
    Names = Parameters.Names.Estimated;
    for j = 1:length(Names)
        Parameters.(Names{j}).Estimated = 0;
    end
    Parameters.(ToEstimate{i}).Estimated = 1;
    Parameters = DefineEstimatedParametersIndexes(Parameters);
    Parameters = DefineTransfFunctions(Parameters);
    Parameters = DefinePriors(Parameters);
    Parameters = UpdateParsNoTransfToTransf(Parameters);
    Initialization = Parameters.(ToEstimate{i}).TransfValue ;
    [x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',2000,'TolX',1e-8,'TolFun',1e-8));
    Parameters.(ToEstimate{i}).TransfValue = (x);
    Parameters = UpdateParsNoTransfToTransf(Parameters);
    Cov = (Parameters.(ToEstimate{i}).TransfValue*0.1)^2;
    Parameters.Epsil = 0.1;
    Parameters.ComputeRWsamplCov = 0;
    Parameters.G = Cov^-1;
    Parameters.NoPaths = 1;
    Parameters.MCMCType = 'Rand';
    Parameters.GMeth = 'cst given';
    Parameters.aim = 0.23;
    [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
    Parameters.Epsil = 1;
    TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
    TempRess{i} = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,1000);
end

Names = Parameters.Names.Estimated;
Cov = zeros(length(Names),length(Names));
for i = 1:length(ToEstimate)
    ind = Parameters.(ToEstimate{i}).Index;
    Parameters.(ToEstimate{i}).Estimated = 1;
    Parameters.(ToEstimate{i}).TransfValue = TempRess{i}.TransfThetas(end);
    Cov(ind,ind) = std(TempRess{i}.TransfThetas)^2;
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsTransfToNoTransf(Parameters);

Parameters.Epsil = 1;
Parameters.ComputeRWsamplCov = 0;
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.aim = 0.23;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
[Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);

TempPar = ResGoog_RW_Calib3.TempPar;
TransfThetasSamples = ResGoog_RW_Calib3.TransfThetas;
Parameters = DefineScalingPars(TransfThetasSamples,Parameters);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,Parameters);
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'GMM';
Parameters.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',Parameters,1);
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);

Parameters.Epsil = 2.38/sqrt(7);
ResGoog_RW_Calib3 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,10000);

TempPar = ResGoog_RW_Calib3.TempPar;
TransfThetasSamples = ResGoog_RW_Calib3.TransfThetas;
Parameters = DefineScalingPars(TransfThetasSamples,Parameters);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,Parameters);
Parameters.MCMCType = 'GMM';
Parameters.GMeth = 'GMM';
Parameters.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',Parameters,1);
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
Parameters.Epsil = 1;
ResGoog_RW_Calib4 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,500);




SavePath = 'C:\Users\dureau\Desktop\Results';
save([SavePath '\ResGoog_AlsaceOct_RW_Calib.mat'],'ResGoog_RW_Calib')


TempPar = ResGoog_RW_Calib.TempPar;
Parameters.NoPaths = 0;
Parameters.PathsToKeep = 1:6;
TransfThetasSamples = ResGoog_RW_Calib.TransfThetas;
Parameters = DefineScalingPars(TransfThetasSamples,Parameters);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,Parameters);
Parameters.MCMCType = 'GMM';
Parameters.GMeth = 'GMM';
Parameters.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',Parameters,3);
Parameters.NoPaths = 1;
Parameters.aim = 0.50;
Parameters.Epsil = 1.3;
[Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);

Res = ResGoog_RW;
autocor = autocorrelation(Res.Thetas(2,:),100);
nESS = length(Res.Thetas(2,:))/(1+2*sum(autocor));
disp(nESS/length(Res.Thetas(2,:)))



SavePath = 'C:\Users\dureau\Desktop\Results';
save([SavePath '\ResGoog_AlsaceOct_RW.mat'],'ResGoog_RW')


Cov = zeros(3,3)
Cov(1,1) = 0.2^2;
Cov(2,2) = 10^2;
Cov(3,3) = 0.25;


Parameters.NbParticules = 500;
Parameters.PathsToKeep = [1:6]';
Parameters.MCMCType = 'Add';
ResultSMC = EstimationSMCsmoothGen(Data, SEIRModel, Parameters)

subplot(2,1,1)
plot(Data.Instants,ResultSMC.PosteriorMeansRecord(3,:))
hold on
plot(Data.Instants,Data.Observations(1,:),'g')
hold off
set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
subplot(2,1,2)
plot(squeeze(mean(exp(ResultSMC.CompletePaths(:,4,:)))))
hold on
plot(cumsum(Data.NbComputingSteps),0*ones(size(Data.NbComputingSteps)),'.k')
set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
title(ResultSMC.LogLik)
pause(0.01)


Res = ResGoog_RW;
% plot ResRun
figure(1)
for i = 1:5
    subplot(6,1,i)
    plot(mean(squeeze(Res.Paths(:,i,:))))
    hold on
%     plot(quantile(squeeze(Res.Paths(:,i,:)),0.025),'r')
%     plot(quantile(squeeze(Res.Paths(:,i,:)),0.975),'r')
    hold off
end
subplot(6,1,6)
plot(mean(squeeze(exp(Res.Paths(:,6,:)))))
hold on
plot(quantile(squeeze(exp(Res.Paths(:,6,:))),0.025),'r')
plot(quantile(squeeze(exp(Res.Paths(:,6,:))),0.975),'r')
hold off    
        
figure(2)
Names = Parameters.Names.Estimated;
k = length(Names);
for i = 1:k
    subplot(ceil(sqrt(k)),ceil(sqrt(k)),i)
    plot(Res.TransfThetas(i,:))   
%     plot(Res.Thetas(i,:))
    title(Names{i})
end
    
        
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).Estimated = 0;
end
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
        
        
Parameters.NbParticules = 10000;
LogLiks = [];
for i = 1:50
    i
    ResultSMC = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);
    LogLiks(i) = ResultSMC.LogLik;
end
hist(LogLiks)

Temp = ResultSMC
subplot(2,1,1)
plot(Data.Instants,Temp.PosteriorMeansRecord(3,:))
hold on
plot(Data.Instants,Data.Observations(5,:),'g')
hold off
set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
subplot(2,1,2)
plot(Data.Instants,exp(Temp.PosteriorMeansRecord(6,:)))
set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
title(Temp.LogLik)
pause(0.01)
        
        
        
        
        
