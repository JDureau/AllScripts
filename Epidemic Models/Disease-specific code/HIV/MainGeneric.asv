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

%% SIR 
InitialState = [10000 100 0 2*10^-4]';


Parameters = struct();
%%% Parameters definition.
Parameters.Names.All = {'Beta', 'Gamma','SigmaObs', 'SigmaRW','SigmaOU','KappaOU','MuOU'};
for i = 1:length(Parameters.Names.All)
    Parameters.(Parameters.Names.All{i}).Estimated = 0;
end

% Estimated parameters
Parameters.Beta.Value = 1;
Parameters.Beta.Transf = @log;
Parameters.Beta.InvTransf = @exp;
Parameters.Beta.Prior = 'Parameters.Beta.Value*normpdf(Parameters.Beta.Value,0,10^12)';
Parameters.Beta.Estimated = 1;
Parameters.Gamma.Value = 1;
Parameters.Gamma.Transf = @log;
Parameters.Gamma.InvTransf = @exp;
Parameters.Gamma.Prior = 'Parameters.Gamma.Value*normpdf(Parameters.Gamma.Value,0,10^12)';
Parameters.Gamma.Estimated = 1;
Parameters.SigmaRW.Value = 10^-5;
Parameters.SigmaRW.Transf = @log;
Parameters.SigmaRW.InvTransf = @exp;
Parameters.SigmaRW.Prior = 'Parameters.SigmaRW.Value*normpdf(Parameters.SigmaRW.Value,0,10^12)';
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
for i = 1:length(Parameters.Names.Estimated)
    Parameters.(Parameters.Names.Estimated{i}).SamplStd = 0.1*Parameters.(Parameters.Names.Estimated{i}).Value;
end

% Other parameters
Parameters.NbVariables = 4;
Parameters.SigmaObs = 0.01;
Parameters.ObsNoise = 0.0;
Parameters.InitialState = InitialState;
Parameters.DiffusionType = 'Add';
Parameters.ObservedVariables = 2;

%%% Model
SIRModel = struct();
SIRModel.EKF_projection = @SIR_EKF_projection;
temp = zeros(1,4);
temp(2) = 1;
SIRModel.InitializeParameters = @SIR_Initialize;
SIRModel.ObservationJacobian = temp;
SIRModel.ObservationMeasurementNoise = Parameters.SigmaObs^2;
SIRModel.SMC_projection = @SIR_SMC_projection;
SIRModel.Lik = @SIR_Lik;

%%% Data
Parameters.ObservationLength = 15;
Parameters.ObservationTStep = 0.1;
Parameters.ComputationTStep = 0.01;

time = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
RealLambdas =1*(2-Sigmoid(time-7))/InitialState(1);
RealTrajectory = ComputeTrajectory(InitialState,RealLambdas,Parameters);

Data.Observations = RealTrajectory.Observations;
Data.Instants = RealTrajectory.ObsTime;
Data.NbComputingSteps = [0 round(diff(Data.Instants)/Parameters.ComputationTStep)];

%%% MLE
Initialization = log([1 1 10^-5]);
Temp = EstimationEKFGen(Data, SIRModel, Parameters);
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimize(x,Data,SIRModel,Parameters),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
TempParameters = Parameters;
TempParameters.Beta.TransfValue = (x(1));
TempParameters.Gamma.TransfValue = (x(2));
TempParameters.SigmaRW.TransfValue = (x(3));
TempParameters = UpdateParsTransfToNoTransf(TempParameters);
Res = KalmanNumericDerivatives(Data,SIRModel,TempParameters);
Parameters.G = -Res.Hess;
Parameters.RWsamplCov = -Res.Hess;
eig(-Res.Hess)>0

% PMCMC
Parameters.Epsil = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.NbParticules = 1000;
Parameters.NoPaths = 1;
TempPar = ProposeInitialParameter(Data, SIRModel, Parameters);
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
[Parameters, TempPar] = CalibrateMethod( Data, SIRModel, Parameters, TempPar);
ResRWML = RunEstimationMethod(Data, SIRModel,Parameters,TempPar,10000);

TransfThetasSamples = ResRWML.TransfThetas;
Parameters = DefineScalingPars(TransfThetasSamples,Parameters);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,Parameters);
Parameters.GMeth = 'GMM';
Parameters.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',Parameters,1);
Parameters.aim = 0.23;
[Parameters, TempPar] = CalibrateMethod( Data, SIRModel, Parameters, TempPar);
ResRW1 = RunEstimationMethod(Data, SIRModel,Parameters,TempPar,10000);


%% HIV Imperial Mysore



Parameters = struct();

% Initializing all model parameters
Parameters.TotalFSW.Value = 2144;
Parameters.TotalFSW.Min = 804;
Parameters.TotalFSW.Max = 3752;
Parameters.TotalFSW.Estimates = 0;
Parameters.TotalFSW.TransfType = 'Log';
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialFt.Value = 0.008;
Parameters.InitialFt.Min = 0;
Parameters.InitialFt.Max = 0.15;
Parameters.InitialFt.Estimated = 1;
Parameters.InitialFt.TransfType = 'Logit';
Parameters.SecondFt.Value = 0.02;
Parameters.SecondFt.Min = 0;
Parameters.SecondFt.Max = 0.9;
Parameters.SecondFt.Estimated = 1;
Parameters.SecondFt.TransfType = 'Logit';
Parameters.ThirdFt.Value = 0.01;
Parameters.ThirdFt.Min = 0;
Parameters.ThirdFt.Max = 0.95;
Parameters.ThirdFt.Estimated = 1;
Parameters.ThirdFt.TransfType = 'Logit';
Parameters.FourthFt.Value = 0.7;
Parameters.FourthFt.Min = 0;
Parameters.FourthFt.Max = 0.95;
Parameters.FourthFt.Estimated = 1;
Parameters.FourthFt.TransfType = 'Logit';
Parameters.InitialSF1 = 0.98*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = 0.02*Parameters.TotF1.Value;
Parameters.InitialSF2 = 0.98*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = 0.02*Parameters.TotF2.Value;
Parameters.TotMFactor.Value = 12;
Parameters.TotMFactor.Min = 7;
Parameters.TotMFactor.Max = 19;
Parameters.TotMFactor.Estimated = 1;
Parameters.TotMFactor.TransfType = 'Log';
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = 0.998*Parameters.TotM.Value;
Parameters.InitialHIVM = 0.002*Parameters.TotM.Value;
Parameters.Alpha.Value = 1/110;
Parameters.Alpha.Min = 1/138.5;
Parameters.Alpha.Max = 1/92;
Parameters.Alpha.Estimated = 1;
Parameters.Alpha.TransfType = 'Log';
Parameters.MuF.Value = 1/60;
Parameters.MuF.Min = 1/66;
Parameters.MuF.Max = 1/53;
Parameters.MuF.Estimated = 1;
Parameters.MuF.TransfType = 'Log';
Parameters.MuM.Value = 1/150;
Parameters.MuM.Min = 1/240;
Parameters.MuM.Max = 1/84;
Parameters.MuM.Estimated = 1;
Parameters.MuM.TransfType = 'Log';
Parameters.BetaMFPerAct.Value = 0.0029;
Parameters.BetaMFPerAct.Min = 0.0006*2;
Parameters.BetaMFPerAct.Max = 0.0011*3;
Parameters.BetaMFPerAct.Estimated = 1;
Parameters.BetaMFPerAct.TransfType = 'Log';
Parameters.BetaFMPerAct.Value = 0.0023;
Parameters.BetaFMPerAct.Min = 0.0001*2;
Parameters.BetaFMPerAct.Max = 0.0014*3;
Parameters.BetaFMPerAct.Estimated = 1;
Parameters.BetaFMPerAct.TransfType = 'Log';
Parameters.NumberActsPerClient.Value = 2;
Parameters.NumberActsPerClient.Min = 1;
Parameters.NumberActsPerClient.Max = 3;
Parameters.NumberActsPerClient.Estimated = 1;
Parameters.NumberActsPerClient.TransfType = 'Log';
Parameters.BetaFM.Value = 1-(1-Parameters.BetaFMPerAct.Value)^Parameters.NumberActsPerClient.Value;
Parameters.BetaMF.Value = 1-(1-Parameters.BetaMFPerAct.Value)^Parameters.NumberActsPerClient.Value;
Parameters.eHIV.Value = 0.82;
Parameters.eHIV.Min = 0.8;
Parameters.eHIV.Max = 0.9;
Parameters.eHIV.Estimated = 1;
Parameters.eHIV.TransfType = 'Log';
Parameters.CF1.Value = 16.26;
Parameters.CF1.Min = 15.22;
Parameters.CF1.Max = 17.3;
Parameters.CF1.Estimated = 1;
Parameters.CF1.TransfType = 'Log';
Parameters.CF2.Value = 51.99;
Parameters.CF2.Min = 47.37;
Parameters.CF2.Max = 56.6;
Parameters.CF2.Estimated = 1;
Parameters.CF2.TransfType = 'Log';
Parameters.SigmaRW.Value = 10;
Parameters.SigmaRW.Min = 0;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Log';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.StableCUseConstraint = 0;
TellParsValues(Parameters)

% Temp = EstimationEKFGen(Data, HIVModel, Parameters);
% PlotHIV(Temp,Data)

% Other parameters
Parameters.NbVariables = 9;
Parameters.SigmaObs = 0.1;
Parameters.DiffusionType = 'Affine';
Parameters.ObservationLength = 25*12;
Parameters.ComputationTStep = 0.5;


% t0 = jan 85.
Data.Observations = zeros(9,5);
Data.Observations(7,2) = 26.11;
Data.Observations(7,3) = 24.24;
Data.Observations(8,4) = 5.4;
Data.Observations(7,5) = 11.10;
Data.Instants = round([0 236 264 286 292]/(Parameters.ComputationTStep));
Data.ObservedVariables = [ 0 7 7 8 7];
Data.NbComputingSteps = [0 diff(Data.Instants)];
% Data.Observations = zeros(9,4);
% Data.Observations(7,2) = 26.11;
% Data.Observations(7,3) = 24.24;
% Data.Observations(8,4) = 5.4;
% Data.Instants = round([0 236 264 286]/(Parameters.ComputationTStep));
% Data.ObservedVariables = [ 0 7 7 8];
% Data.NbComputingSteps = [0 diff(Data.Instants)];


%%% Model
HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIV_Initialize;
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*Parameters.SigmaObs).*(Res.WentOutOrNot)';
HIVModel.SMC_projection = @HIV_SMC_projection;

temp7 = zeros(1,9);
temp8 = zeros(1,9);
temp7(1,7) = 1;
temp8(1,8) = 1;
HIVModel.ObservationJacobian = {};
HIVModel.ObservationJacobian{2} = temp7;
HIVModel.ObservationJacobian{3} = temp7;
HIVModel.ObservationJacobian{4} = temp8;
HIVModel.ObservationJacobian{5} = temp7;
HIVModel.ObservationMeasurementNoise = {};
HIVModel.ObservationMeasurementNoise{2} = (Parameters.SigmaObs*Data.Observations(7,2))^2;
HIVModel.ObservationMeasurementNoise{3} = (Parameters.SigmaObs*Data.Observations(7,3))^2;
HIVModel.ObservationMeasurementNoise{4} = (Parameters.SigmaObs*Data.Observations(8,4))^2;
HIVModel.ObservationMeasurementNoise{5} = (Parameters.SigmaObs*Data.Observations(7,5))^2;


% Kal Opt to find MLE parameters (no Cov so no update / rough variations of Ft)
Names = Parameters.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Parameters.(Names{i}).Estimated = 1;
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters.DiffusionType = 'Affine';

[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,Parameters),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
Initialization = x;
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)

%%% Let's start PMCMC
ParametersPMCMC = Parameters;

% first, model parameters are fixed to MLE, we just estimate the diffusion
% parameters
ParametersPMCMC.NbParticules = 1000;
ParametersPMCMC.StableCUseConstraint = 0;
ParametersPMCMC.NoPaths = 0;
ParametersPMCMC.PathsToKeep = [7 8 9];
ParametersPMCMC.DiffusionType = 'Add';
ParametersPMCMC.MCMCType = 'jiji';
ParametersPMCMC.GMeth = 'cst given';
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end

ParametersPMCMC.SigmaRW.Value = 0.1;
ParametersPMCMC.SigmaRW.Min = 0;
ParametersPMCMC.SigmaRW.Max = 10^14;
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC.SigmaRW.TransfType = 'Log';
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);
ResultSMC = EstimationSMCsmoothGen(Data, HIVModel, ParametersPMCMC)
PlotresHIV(ResultSMC)

% optimize RW parameters with SMC
ParametersPMCMC.NoPaths = 0;
Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)




% optimize all parameters.
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC.TotMFactor.Estimated = 1;
ParametersPMCMC.Alpha.Estimated = 1;
ParametersPMCMC.MuF.Estimated = 1;
ParametersPMCMC.MuM.Estimated = 1;
ParametersPMCMC.BetaMFPerAct.Estimated = 1;
ParametersPMCMC.BetaFMPerAct.Estimated = 1;
ParametersPMCMC.NumberActsPerClient.Estimated = 1;
ParametersPMCMC.eHIV.Estimated = 1;
ParametersPMCMC.CF1.Estimated = 1;
ParametersPMCMC.CF2.Estimated = 1;
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);
Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',2000));
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)

% Real PMCMC
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end
ParametersPMCMC.SigmaRW.Estimated = 1;
% ParametersPMCMC.NumberActsPerClient.Estimated = 1;
% ParametersPMCMC.eHIV.Estimated = 1;
% ParametersPMCMC.Alpha.Estimated = 1;
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);


Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC.DiffusionType = 'Affine';

[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
Initialization = x;
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).SamplStd = 0.01*ParametersPMCMC.(Names{i}).Value;
end
Res = KalmanNumericDerivativesWithPrior(Data,HIVModel,ParametersPMCMC);
Test = mean(eig(-Res.Hess)>0)==1;
Test


% Define Cov
Cov = [];
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    Cov(ParametersPMCMC.(Names{i}).Index,ParametersPMCMC.(Names{i}).Index) =  (0.1*ParametersPMCMC.(Names{i}).TransfValue)^2;
end

% Estimating only RW parameters:
ParametersPMCMC.Epsil = 1;
ParametersPMCMC.ComputeRWsamplCov = 0;
ParametersPMCMC.G = -Res.Hess;%Cov^-1;
ParametersPMCMC.NoPaths = 1;
ParametersPMCMC.MCMCType = 'Rand';
ParametersPMCMC.GMeth = 'cst given';
ParametersPMCMC.DiffusionType = 'Add';
ParametersPMCMC.aim = 0.23;
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
% [Parameters, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ParametersPMCMC.Epsil = 1;
ResRW1_RW_Calib_Rds_1_2_3 = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

TempPar = ResRW1_RW_Calib_Rds_1_2_3.TempPar;
ParametersPMCMC.NoPaths = 1;
ParametersPMCMC.PathsToKeep = [7 8 9];
TransfThetasSamples = ResRW1_RW_Calib_Rds_1_2_3.TransfThetas;
ParametersPMCMC = DefineScalingPars(TransfThetasSamples,ParametersPMCMC);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,ParametersPMCMC);
ParametersPMCMC.GMeth = 'GMM';
ParametersPMCMC.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',ParametersPMCMC,1);
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
ParametersPMCMC.Epsil = 0.2;
ResRW1_RW_Calib_Rds_1_2_3 = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);


SavePath = 'C:\Users\dureau\Desktop\Results';
% SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_Mysore_withPars_Calib.mat'],'ResRW1_RW_Calib_Rds_1_2_3')


TempPar = ResRW1_RW_Calib_Rds_1_2_3.TempPar;
ParametersPMCMC.NoPaths = 0;
ParametersPMCMC.PathsToKeep = [1:9];
TransfThetasSamples = ResRW1_RW_Calib_Rds_1_2_3.TransfThetas;
ParametersPMCMC = DefineScalingPars(TransfThetasSamples,ParametersPMCMC);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,ParametersPMCMC);
ParametersPMCMC.GMeth = 'GMM';
ParametersPMCMC.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',ParametersPMCMC,1);
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
ParametersPMCMC.Epsil = 1;
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_RW_Rds_1_2 = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);


SavePath = 'C:\Users\dureau\Desktop\Results';
% SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_Mysore_sigmaonly_forprojections.mat'],'ResRW1_RW_Rds_1_2')

%% End copypast

% Kal Opt to find MLE parameters (no Cov so no update / rough variations of Ft)
Names = Parameters.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end

[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,Parameters),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
Initialization = x;
Initialization = x;
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)

%%% Let's start PMCMC
ParametersPMCMC = Parameters;

% first, model parameters are fixed to MLE, we just estimate the diffusion
% parameters
ParametersPMCMC.NbParticules = 1000;
ParametersPMCMC.StableCUseConstraint = 0;
ParametersPMCMC.NoPaths = 0;
ParametersPMCMC.PathsToKeep = [7 8 9];
ParametersPMCMC.DiffusionType = 'Add';
ParametersPMCMC.MCMCType = 'jiji';
ParametersPMCMC.GMeth = 'cst given';
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end

ParametersPMCMC.SigmaRW.Value = 0.1;
ParametersPMCMC.SigmaRW.Min = 0;
ParametersPMCMC.SigmaRW.Max = 10^14;
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC.SigmaRW.TransfType = 'Log';
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);
ResultSMC = EstimationSMCsmoothGen(Data, HIVModel, ParametersPMCMC)
PlotresHIV(ResultSMC)

% optimize RW parameters with SMC
ParametersPMCMC.NoPaths = 0;
Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)

% optimize all parameters.
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC.TotMFactor.Estimated = 1;
ParametersPMCMC.Alpha.Estimated = 1;
ParametersPMCMC.MuF.Estimated = 1;
ParametersPMCMC.MuM.Estimated = 1;
ParametersPMCMC.BetaMFPerAct.Estimated = 1;
ParametersPMCMC.BetaFMPerAct.Estimated = 1;
ParametersPMCMC.NumberActsPerClient.Estimated = 1;
ParametersPMCMC.eHIV.Estimated = 1;
ParametersPMCMC.CF1.Estimated = 1;
ParametersPMCMC.CF2.Estimated = 1;
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);
Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',2000));
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)

% Real PMCMC
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);

% Define Cov
Cov = [];
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    Cov(ParametersPMCMC.(Names{i}).Index,ParametersPMCMC.(Names{i}).Index) =  (0.1*ParametersPMCMC.(Names{i}).TransfValue)^2;
end

% Estimating only RW parameters:
ParametersPMCMC.Epsil = 1;
ParametersPMCMC.ComputeRWsamplCov = 0;
ParametersPMCMC.G = Cov^-1;
ParametersPMCMC.NoPaths = 1;
ParametersPMCMC.MCMCType = 'Rand';
ParametersPMCMC.GMeth = 'cst given';
ParametersPMCMC.aim = 0.23;
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
[Parameters, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_RW_Calib_Rds_1_2_3 = RunEstimationMethod(Data, HIVModel,Parameters,TempPar,5000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_ResRW1_RW_Calib.mat'],'ResRW1_RW_Calib_Rds_1_2_3')


TempPar = ResRW1_RW_Calib_Rds_1_2_3.TempPar;
ParametersPMCMC.NoPaths = 0;
ParametersPMCMC.PathsToKeep = [7 8 9];
TransfThetasSamples = ResRW1_RW_Calib_Rds_1_2_3.TransfThetas;
ParametersPMCMC = DefineScalingPars(TransfThetasSamples,ParametersPMCMC);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,ParametersPMCMC);
ParametersPMCMC.GMeth = 'GMM';
ParametersPMCMC.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',ParametersPMCMC,1);
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_RW_Rds_1_2 = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_ResRW1_RW.mat'],'ResRW1_RW_Rds_1_2_3')

% Now introduce some of the model parameters
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end

ParametersPMCMC.DiffusionType = 'Affine';
ParametersPMCMC.InitialFt.Estimated = 1;
ParametersPMCMC.SecondFt.Estimated = 1;
ParametersPMCMC.ThirdFt.Estimated = 1;
ParametersPMCMC.FourthFt.Estimated = 1;
ParametersPMCMC.NumberActsPerClient.Estimated = 1;
ParametersPMCMC.NumberActsPerClient.TransfType = 'Log';
ParametersPMCMC.InitialFt.TransfType = 'Log';
ParametersPMCMC.InitialFt.Estimated = 1;
ParametersPMCMC.eHIV.Estimated = 1;
ParametersPMCMC.eHIV.TransfType = 'Log';
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);

Names = ParametersPMCMC.Names.Estimated;
Cov = [];
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).SamplStd = 0.1*ParametersPMCMC.(Names{i}).Value;
end

Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
Initialization = x;
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)


Res = KalmanNumericDerivativesWithPrior(Data,HIVModel,ParametersPMCMC);

tempCov = (-Res.Hess)^-1;
tempCov(end+1,end+1) = (std(ResRW1_RW.Thetas))^2;

Parameters.G = tempCov^-1;
Parameters.RWsamplCov = tempCov^-1;
eig(tempCov^-1)>0

ParametersPMCMC.DiffusionType = 'Add';
ParametersPMCMC.InitialFt.Estimated = 0;
ParametersPMCMC.SecondFt.Estimated = 0;
ParametersPMCMC.ThirdFt.Estimated = 0;
ParametersPMCMC.FourthFt.Estimated = 0;
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);


ParametersPMCMC.Epsil = 1;
ParametersPMCMC.ComputeRWsamplCov = 0;
ParametersPMCMC.NoPaths = 1;
ParametersPMCMC.MCMCType = 'Rand';
ParametersPMCMC.GMeth = 'cst given';
ParametersPMCMC.aim = 0.23;
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_RW_andPars_Calib = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_ResRW1_RW_andPars_Calib.mat'],'ResRW1_RW_andPars_Calib')


TempPar = ResRW1_RW_andPars_Calib.TempPar;
ParametersPMCMC.NoPaths = 0;
TransfThetasSamples = ResRW1_RW_andPars_Calib.TransfThetas;
ParametersPMCMC = DefineScalingPars(TransfThetasSamples,ParametersPMCMC);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,ParametersPMCMC);
ParametersPMCMC.GMeth = 'GMM';
ParametersPMCMC.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',ParametersPMCMC,1);
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_andPars_RW = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_ResRW1_RW_andPars.mat'],'ResRW1_RW_andPars')




clf
Names = Parameters.Names.Estimated;
k = length(Names);
for i = 1:k
    subplot(ceil(sqrt(k)),ceil(sqrt(k)),i)
%     plot(ResRWML.Thetas(i,:),ResRWML.Paths(:,8,end),'.')
    plot(ResRW1_RW_andPars_Calib.Thetas(i,:))
    hold on
    title(Names{i})
end




%% HIV Imperial Belgaum



Parameters = struct();

% Initializing all model parameters
Parameters.TotalFSW.Value = 1742;
Parameters.TotalFSW.Min = 804;
Parameters.TotalFSW.Max = 3752;
Parameters.TotalFSW.Estimates = 0;
Parameters.TotalFSW.TransfType = 'Log';
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialFt.Value = 0.15;
Parameters.InitialFt.Min = 0;
Parameters.InitialFt.Max = 0.9;
Parameters.InitialFt.Estimated = 1;
Parameters.InitialFt.TransfType = 'Logit';
Parameters.SecondFt.Value = 0.02;
Parameters.SecondFt.Min = 0;
Parameters.SecondFt.Max = 0.9;
Parameters.SecondFt.Estimated = 1;
Parameters.SecondFt.TransfType = 'Logit';
Parameters.ThirdFt.Value = 0.01;
Parameters.ThirdFt.Min = 0;
Parameters.ThirdFt.Max = 0.95;
Parameters.ThirdFt.Estimated = 1;
Parameters.ThirdFt.TransfType = 'Logit';
Parameters.FourthFt.Value = 0.7;
Parameters.FourthFt.Min = 0;
Parameters.FourthFt.Max = 0.95;
Parameters.FourthFt.Estimated = 1;
Parameters.FourthFt.TransfType = 'Logit';
Parameters.InitialSF1 = 0.98*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = 0.02*Parameters.TotF1.Value;
Parameters.InitialSF2 = 0.98*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = 0.02*Parameters.TotF2.Value;
Parameters.TotMFactor.Value = 21;
Parameters.TotMFactor.Min = 7;
Parameters.TotMFactor.Max = 35;
Parameters.TotMFactor.Estimated = 1;
Parameters.TotMFactor.TransfType = 'Log';
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = 0.998*Parameters.TotM.Value;
Parameters.InitialHIVM = 0.002*Parameters.TotM.Value;
Parameters.Alpha.Value = 1/110;
Parameters.Alpha.Min = 1/138.5;
Parameters.Alpha.Max = 1/92;
Parameters.Alpha.Estimated = 1;
Parameters.Alpha.TransfType = 'Log';
Parameters.MuF.Value = 1/120;
Parameters.MuF.Min = 1/250;
Parameters.MuF.Max = 1/94;
Parameters.MuF.Estimated = 1;
Parameters.MuF.TransfType = 'Log';
Parameters.MuM.Value = 1/120;
Parameters.MuM.Min = 1/143;
Parameters.MuM.Max = 1/83;
Parameters.MuM.Estimated = 1;
Parameters.MuM.TransfType = 'Log';
Parameters.BetaMFPerAct.Value = 0.0029;
Parameters.BetaMFPerAct.Min = 0.0006*2;
Parameters.BetaMFPerAct.Max = 0.0011*3;
Parameters.BetaMFPerAct.Estimated = 1;
Parameters.BetaMFPerAct.TransfType = 'Log';
Parameters.BetaFMPerAct.Value = 0.0023;
Parameters.BetaFMPerAct.Min = 0.0001*2;
Parameters.BetaFMPerAct.Max = 0.0014*3;
Parameters.BetaFMPerAct.Estimated = 1;
Parameters.BetaFMPerAct.TransfType = 'Log';
Parameters.NumberActsPerClient.Value = 2;
Parameters.NumberActsPerClient.Min = 1;
Parameters.NumberActsPerClient.Max = 3;
Parameters.NumberActsPerClient.Estimated = 1;
Parameters.NumberActsPerClient.TransfType = 'Log';
Parameters.BetaFM.Value = 1-(1-Parameters.BetaFMPerAct.Value)^Parameters.NumberActsPerClient.Value;
Parameters.BetaMF.Value = 1-(1-Parameters.BetaMFPerAct.Value)^Parameters.NumberActsPerClient.Value;
Parameters.eHIV.Value = 0.82;
Parameters.eHIV.Min = 0.8;
Parameters.eHIV.Max = 0.9;
Parameters.eHIV.Estimated = 1;
Parameters.eHIV.TransfType = 'Log';
Parameters.CF1.Value = 23;
Parameters.CF1.Min = 22.1;
Parameters.CF1.Max = 25.1;
Parameters.CF1.Estimated = 1;
Parameters.CF1.TransfType = 'Log';
Parameters.CF2.Value = 95;
Parameters.CF2.Min = 84.8;
Parameters.CF2.Max = 103.8;
Parameters.CF2.Estimated = 1;
Parameters.CF2.TransfType = 'Log';
Parameters.SigmaRW.Value = 10;
Parameters.SigmaRW.Min = 0;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.Estimated = 0;
Parameters.SigmaRW.TransfType = 'Log';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.StableCUseConstraint = 0;
TellParsValues(Parameters)

% Temp = EstimationEKFGen(Data, HIVModel, Parameters);
% PlotHIV(Temp,Data)

% Other parameters
Parameters.NbVariables = 9;
Parameters.SigmaObs = 0.1;
Parameters.DiffusionType = 'Affine';
Parameters.ObservationLength = 26*12;
Parameters.ComputationTStep = 0.5;


% t0 = jan 85.
Data.Observations = zeros(9,5);
Data.Observations(7,2) = 33.9;
Data.Observations(7,3) = 27.3;
Data.Observations(8,4) = 6.2;
Data.Observations(7,5) = 22.3;
Data.Instants = round([0 250 283 286 309]/(Parameters.ComputationTStep));
Data.ObservedVariables = [ 0 7 7 8 7];
Data.NbComputingSteps = [0 diff(Data.Instants)];
% Data.Observations = zeros(9,4);
% Data.Observations(7,2) = 33.9;
% Data.Observations(7,3) = 27.3;
% Data.Observations(8,4) = 6.2;
% Data.Instants = round([0 250 283 286]/(Parameters.ComputationTStep));
% Data.ObservedVariables = [ 0 7 7 8];
% Data.NbComputingSteps = [0 diff(Data.Instants)];


%%% Model
HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIVInitialize;
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*Parameters.SigmaObs).*(Res.WentOutOrNot)';
HIVModel.SMC_projection = @HIV_SMC_projection;

temp7 = zeros(1,9);
temp8 = zeros(1,9);
temp7(1,7) = 1;
temp8(1,8) = 1;
HIVModel.ObservationJacobian = {};
HIVModel.ObservationJacobian{2} = temp7;
HIVModel.ObservationJacobian{3} = temp7;
HIVModel.ObservationJacobian{4} = temp8;
HIVModel.ObservationJacobian{5} = temp7;
HIVModel.ObservationMeasurementNoise = {};
HIVModel.ObservationMeasurementNoise{2} = (Parameters.SigmaObs*Data.Observations(7,2))^2;
HIVModel.ObservationMeasurementNoise{3} = (Parameters.SigmaObs*Data.Observations(7,3))^2;
HIVModel.ObservationMeasurementNoise{4} = (Parameters.SigmaObs*Data.Observations(8,4))^2;
HIVModel.ObservationMeasurementNoise{5} = (Parameters.SigmaObs*Data.Observations(7,5))^2;


% Kal Opt to find MLE parameters (no Cov so no update / rough variations of Ft)
Names = Parameters.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end

[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,Parameters),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
Initialization = x;
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)

%%% Let's start PMCMC
ParametersPMCMC = Parameters;

% first, model parameters are fixed to MLE, we just estimate the diffusion
% parameters
ParametersPMCMC.NbParticules = 1000;
ParametersPMCMC.StableCUseConstraint = 0;
ParametersPMCMC.NoPaths = 0;
ParametersPMCMC.PathsToKeep = [7 8 9];
ParametersPMCMC.DiffusionType = 'Add';
ParametersPMCMC.MCMCType = 'jiji';
ParametersPMCMC.GMeth = 'cst given';
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end

ParametersPMCMC.SigmaRW.Value = 0.1;
ParametersPMCMC.SigmaRW.Min = 0;
ParametersPMCMC.SigmaRW.Max = 10^14;
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC.SigmaRW.TransfType = 'Log';
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);
ResultSMC = EstimationSMCsmoothGen(Data, HIVModel, ParametersPMCMC)
PlotresHIV(ResultSMC)

% optimize RW parameters with SMC
ParametersPMCMC.NoPaths = 0;
Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)

% optimize all parameters.
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC.TotMFactor.Estimated = 1;
ParametersPMCMC.Alpha.Estimated = 1;
ParametersPMCMC.MuF.Estimated = 1;
ParametersPMCMC.MuM.Estimated = 1;
ParametersPMCMC.BetaMFPerAct.Estimated = 1;
ParametersPMCMC.BetaFMPerAct.Estimated = 1;
ParametersPMCMC.NumberActsPerClient.Estimated = 1;
ParametersPMCMC.InitialFt.Estimated = 1;
ParametersPMCMC.eHIV.Estimated = 1;
ParametersPMCMC.CF1.Estimated = 1;
ParametersPMCMC.CF2.Estimated = 1;
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);
Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',2000));
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)

% Real PMCMC
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);

% Define Cov
Cov = [];
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    Cov(ParametersPMCMC.(Names{i}).Index,ParametersPMCMC.(Names{i}).Index) =  (0.1*ParametersPMCMC.(Names{i}).TransfValue)^2;
end

% Estimating only RW parameters:
ParametersPMCMC.Epsil = 1;
ParametersPMCMC.ComputeRWsamplCov = 0;
ParametersPMCMC.G = Cov^-1;
ParametersPMCMC.NoPaths = 1;
ParametersPMCMC.MCMCType = 'Rand';
ParametersPMCMC.GMeth = 'cst given';
ParametersPMCMC.aim = 0.23;
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_RW_Calib_Rds_1_2_3 = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_Belgaum_Calib_ResRW1_RW_123.mat'],'ResRW1_RW_Calib_Rds_1_2_3')


TempPar = ResRW1_RW_Calib_Rds_1_2_3.TempPar;
ParametersPMCMC.NoPaths = 0;
ParametersPMCMC.PathsToKeep = [7 8 9];
TransfThetasSamples = ResRW1_RW_Calib_Rds_1_2_3.TransfThetas;
ParametersPMCMC = DefineScalingPars(TransfThetasSamples,ParametersPMCMC);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,ParametersPMCMC);
ParametersPMCMC.GMeth = 'GMM';
ParametersPMCMC.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',ParametersPMCMC,1);
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_RW_Rds_1_2_3 = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,10000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_Belgaum_ResRW1_RW_123.mat'],'ResRW1_RW_Rds_1_2_3')

% Now introduce some of the model parameters
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end

ParametersPMCMC.DiffusionType = 'Affine';
ParametersPMCMC.InitialFt.Estimated = 1;
ParametersPMCMC.SecondFt.Estimated = 1;
ParametersPMCMC.ThirdFt.Estimated = 1;
ParametersPMCMC.FourthFt.Estimated = 1;
ParametersPMCMC.NumberActsPerClient.Estimated = 1;
ParametersPMCMC.NumberActsPerClient.TransfType = 'Log';
ParametersPMCMC.InitialFt.TransfType = 'Log';
ParametersPMCMC.InitialFt.Estimated = 1;
ParametersPMCMC.eHIV.Estimated = 1;
ParametersPMCMC.eHIV.TransfType = 'Log';
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);

Names = ParametersPMCMC.Names.Estimated;
Cov = [];
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).SamplStd = 0.1*ParametersPMCMC.(Names{i}).Value;
end

Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
Initialization = x;
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)


Res = KalmanNumericDerivativesWithPrior(Data,HIVModel,ParametersPMCMC);

tempCov = (-Res.Hess)^-1;
tempCov(end+1,end+1) = (std(ResRW1_RW.Thetas))^2;

Parameters.G = tempCov^-1;
Parameters.RWsamplCov = tempCov^-1;
eig(tempCov^-1)>0

ParametersPMCMC.DiffusionType = 'Add';
ParametersPMCMC.InitialFt.Estimated = 0;
ParametersPMCMC.SecondFt.Estimated = 0;
ParametersPMCMC.ThirdFt.Estimated = 0;
ParametersPMCMC.FourthFt.Estimated = 0;
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);


ParametersPMCMC.Epsil = 1;
ParametersPMCMC.ComputeRWsamplCov = 0;
ParametersPMCMC.NoPaths = 1;
ParametersPMCMC.MCMCType = 'Rand';
ParametersPMCMC.GMeth = 'cst given';
ParametersPMCMC.aim = 0.23;
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_RW_andPars_Calib = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_ResRW1_RW_andPars_Calib.mat'],'ResRW1_RW_andPars_Calib')


TempPar = ResRW1_RW_andPars_Calib.TempPar;
ParametersPMCMC.NoPaths = 0;
TransfThetasSamples = ResRW1_RW_andPars_Calib.TransfThetas;
ParametersPMCMC = DefineScalingPars(TransfThetasSamples,ParametersPMCMC);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,ParametersPMCMC);
ParametersPMCMC.GMeth = 'GMM';
ParametersPMCMC.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',ParametersPMCMC,1);
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_andPars_RW = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_ResRW1_RW_andPars.mat'],'ResRW1_RW_andPars')




clf
Names = Parameters.Names.Estimated;
k = length(Names);
for i = 1:k
    subplot(ceil(sqrt(k)),ceil(sqrt(k)),i)
%     plot(ResRWML.Thetas(i,:),ResRWML.Paths(:,8,end),'.')
    plot(ResRW1_RW_andPars_Calib.Thetas(i,:))
    hold on
    title(Names{i})
end


%% HIV Imperial Tests
SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
load([SavePath '\HIV_Mysore_ResRW1_RW_andPars.mat'])


Parameters = ResRW1_RW_Rds_1_2_3.Parameters;
Data = ResRW1_RW_Rds_1_2_3.Data;

% t = 0:1:Data.Instants(end);
% FtPath = 0.05 + 0.8*sigmoid((t-3*t(end)/4)/50); Test Mysore 1
t = 0:1:1300;
Parameters.ComputationTStep = 1;
FtPath = 0.3*ones(size(t)); 
Parameters.InitialFt.Value = 0.3;
plot(t,FtPath)
record = [];
Parameters = HIV_Initialize(Parameters);
mpred = Parameters.InitialState;
TStep = Parameters.ComputationTStep;
for IndDiscr = 1:length(t);
    mtemp = mpred;
    TotF1 = mtemp(1) + mtemp(2);
    TotF2 = mtemp(3) + mtemp(4);
    TotM  = mtemp(5) + mtemp(6);
    mtemp(9) = FtPath(IndDiscr);
    mpred(1) = mpred(1) + ( Parameters.MuF.Value*(TotF1) + Parameters.Alpha.Value*mtemp(2) - Parameters.BetaMF.Value*Parameters.CF1.Value*(1-Parameters.eHIV.Value*mtemp(9))*mtemp(1).*mtemp(6)./TotM-Parameters.MuF.Value*mtemp(1))*TStep;
    mpred(2) = mpred(2) + ( Parameters.BetaMF.Value*Parameters.CF1.Value*(1-Parameters.eHIV.Value*mtemp(9))*mtemp(1).*mtemp(6)./TotM - (Parameters.MuF.Value + Parameters.Alpha.Value)*mtemp(2))*TStep;
    mpred(3) = mpred(3) + ( Parameters.MuF.Value*(TotF2) + Parameters.Alpha.Value*mtemp(4) - Parameters.BetaMF.Value*Parameters.CF2.Value*(1-Parameters.eHIV.Value*mtemp(9))*mtemp(3).*mtemp(6)./TotM-Parameters.MuF.Value*mtemp(3))*TStep;
    mpred(4) = mpred(4) + ( Parameters.BetaMF.Value*Parameters.CF2.Value*(1-Parameters.eHIV.Value*mtemp(9))*mtemp(3).*mtemp(6)./TotM - (Parameters.MuF.Value + Parameters.Alpha.Value)*mtemp(4))*TStep;
    mpred(5) = mpred(5) + ( Parameters.MuM.Value*(TotM) + Parameters.Alpha.Value*mtemp(6) - Parameters.BetaFM.Value*(1-Parameters.eHIV.Value*mtemp(9))*mtemp(5).*(Parameters.CF1.Value*TotF1./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*mtemp(2)./TotF1+Parameters.CF2.Value*TotF2./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*mtemp(4)./TotF2)-Parameters.MuM.Value*mtemp(5))*TStep;
    mpred(6) = mpred(6) + ( Parameters.BetaFM.Value*(1-Parameters.eHIV.Value*mtemp(9))*mtemp(5).*(Parameters.CF1.Value*TotF1./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*mtemp(2)./TotF1+Parameters.CF2.Value*TotF2./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*mtemp(4)./TotF2) - (Parameters.MuM.Value+Parameters.Alpha.Value)*mtemp(6))*TStep;
    mpred(7) = (mpred(2) + mpred(4))/(TotF1+TotF2)*100; 
    mpred(8) = (mpred(6))/(TotM)*100; 
    record(:,IndDiscr) = mpred;
end
plot(record(7,:))

% obs 2: late
% obs 1: early
Data.Instants = round([0 1000 1100 1200 1300]/(Parameters.ComputationTStep)); %Mysore30 Obs 2
% Data.Instants = round([0 600 700 750 800]/(Parameters.ComputationTStep)); %Mysore5 Obs 2
% Data.Instants = round([0 200 400 500 600]/(Parameters.ComputationTStep)); %Mysore5 Obs 1
% Data.Instants = round([0 120 264 286 292]/(Parameters.ComputationTStep)); %Mysore1 Obs 1
% Data.Instants = round([0 236 264 286 292]/(Parameters.ComputationTStep));
% %Mysore1 Obs 2
% t0 = jan 85.
Data.Observations = zeros(9,5);
Data.Observations(7,2) = record(7,Data.Instants(2));
Data.Observations(7,3) = record(7,Data.Instants(3));
Data.Observations(8,4) = record(8,Data.Instants(4));
Data.Observations(7,5) = record(7,Data.Instants(5));
Data.ObservedVariables = [ 0 7 7 8 7];
Data.NbComputingSteps = [0 diff(Data.Instants)];


%%% Model
HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIV_Initialize;
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*Parameters.SigmaObs).*(Res.WentOutOrNot)';
HIVModel.SMC_projection = @HIV_SMC_projection;

temp7 = zeros(1,9);
temp8 = zeros(1,9);
temp7(1,7) = 1;
temp8(1,8) = 1;
HIVModel.ObservationJacobian = {};
HIVModel.ObservationJacobian{2} = temp7;
HIVModel.ObservationJacobian{3} = temp7;
HIVModel.ObservationJacobian{4} = temp8;
HIVModel.ObservationJacobian{5} = temp7;
HIVModel.ObservationMeasurementNoise = {};
HIVModel.ObservationMeasurementNoise{2} = (Parameters.SigmaObs*Data.Observations(7,2))^2;
HIVModel.ObservationMeasurementNoise{3} = (Parameters.SigmaObs*Data.Observations(7,3))^2;
HIVModel.ObservationMeasurementNoise{4} = (Parameters.SigmaObs*Data.Observations(8,4))^2;
HIVModel.ObservationMeasurementNoise{5} = (Parameters.SigmaObs*Data.Observations(7,5))^2;


% Kal Opt to find MLE parameters (no Cov so no update / rough variations of Ft)
% Names = Parameters.Names.All;
% Initialization = [];
% for i = 1:length(Names)
%     Parameters.(Names{i}).Estimated = 1;
%     Initialization(i) = Parameters.(Names{i}).TransfValue ;
% end
% Parameters = DefineEstimatedParametersIndexes(Parameters);
% Parameters.DiffusionType = 'Affine';
% 
% [x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,Parameters),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
% Initialization = x;
% Initialization = x;
% Names = Parameters.Names.Estimated;
% for i = 1:length(Names)
%     Parameters.(Names{i}).TransfValue = (x(i));
% end
% Parameters = UpdateParsTransfToNoTransf(Parameters);
% TellParsValues(Parameters)

%%% Let's start PMCMC
ParametersPMCMC = Parameters;

% first, model parameters are fixed to MLE, we just estimate the diffusion
% parameters
ParametersPMCMC.NbParticules = 1000;
ParametersPMCMC.StableCUseConstraint = 0;
ParametersPMCMC.NoPaths = 0;
ParametersPMCMC.PathsToKeep = [7 8 9];
ParametersPMCMC.DiffusionType = 'Add';
ParametersPMCMC.MCMCType = 'jiji';
ParametersPMCMC.GMeth = 'cst given';
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end
ParametersPMCMC.SigmaRW.Value = 0.1;
ParametersPMCMC.SigmaRW.Min = 0;
ParametersPMCMC.SigmaRW.Max = 10^14;
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC.SigmaRW.TransfType = 'Log';
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);
ResultSMC = EstimationSMCsmoothGen(Data, HIVModel, ParametersPMCMC)
PlotresHIV(ResultSMC)

% optimize RW parameters with SMC
ParametersPMCMC.NoPaths = 0;
Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)

% optimize all parameters.
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC.TotMFactor.Estimated = 1;
ParametersPMCMC.Alpha.Estimated = 1;
ParametersPMCMC.MuF.Estimated = 1;
ParametersPMCMC.MuM.Estimated = 1;
ParametersPMCMC.BetaMFPerAct.Estimated = 1;
ParametersPMCMC.BetaFMPerAct.Estimated = 1;
ParametersPMCMC.NumberActsPerClient.Estimated = 1;
ParametersPMCMC.eHIV.Estimated = 1;
ParametersPMCMC.CF1.Estimated = 1;
ParametersPMCMC.CF2.Estimated = 1;
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);
Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',2000));
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)

% Real PMCMC
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC.NumberActsPerClient.Estimated = 1;
ParametersPMCMC.eHIV.Estimated = 1;
ParametersPMCMC.Alpha.Estimated = 1;
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);

% Define Cov
Cov = [];
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    Cov(ParametersPMCMC.(Names{i}).Index,ParametersPMCMC.(Names{i}).Index) =  (0.1*ParametersPMCMC.(Names{i}).TransfValue)^2;
end

% Estimating only RW parameters:
ParametersPMCMC.Epsil = 1;
ParametersPMCMC.ComputeRWsamplCov = 0;
ParametersPMCMC.G = Cov^-1;
ParametersPMCMC.NoPaths = 1;
ParametersPMCMC.MCMCType = 'Rand';
ParametersPMCMC.GMeth = 'cst given';
ParametersPMCMC.aim = 0.23;
ParametersPMCMC.NoPaths = 1;
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
[Parameters, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ParametersPMCMC.Epsil = 8;
ResRW1_RW_Calib_Rds_1_2_3 = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_Tests_Mysore05_Obs2.mat'],'ResRW1_RW_Calib_Rds_1_2_3')


TempPar = ResRW1_RW_Calib_Rds_1_2_3.TempPar;
ParametersPMCMC.NoPaths = 0;
ParametersPMCMC.PathsToKeep = [7 8 9];
TransfThetasSamples = ResRW1_RW_Calib_Rds_1_2_3.TransfThetas;
ParametersPMCMC = DefineScalingPars(TransfThetasSamples,ParametersPMCMC);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,ParametersPMCMC);
ParametersPMCMC.GMeth = 'GMM';
ParametersPMCMC.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',ParametersPMCMC,1);
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_RW_Rds_1_2 = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

ResRW1_RW_Rds_1_2.FtPath = FtPath;
SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_Tests_Mysore05_Obs2.mat'],'ResRW1_RW_Rds_1_2')

% Now introduce some of the model parameters
Names = ParametersPMCMC.Names.All;
for i = 1:length(Names);
    ParametersPMCMC.(Names{i}).Estimated = 0;
end

ParametersPMCMC.DiffusionType = 'Affine';
ParametersPMCMC.InitialFt.Estimated = 1;
ParametersPMCMC.SecondFt.Estimated = 1;
ParametersPMCMC.ThirdFt.Estimated = 1;
ParametersPMCMC.FourthFt.Estimated = 1;
ParametersPMCMC.NumberActsPerClient.Estimated = 1;
ParametersPMCMC.NumberActsPerClient.TransfType = 'Log';
ParametersPMCMC.InitialFt.TransfType = 'Log';
ParametersPMCMC.InitialFt.Estimated = 1;
ParametersPMCMC.eHIV.Estimated = 1;
ParametersPMCMC.eHIV.TransfType = 'Log';
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);

Names = ParametersPMCMC.Names.Estimated;
Cov = [];
for i = 1:length(Names)
    ParametersPMCMC.(Names{i}).SamplStd = 0.1*ParametersPMCMC.(Names{i}).Value;
end

Names = ParametersPMCMC.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersPMCMC.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,ParametersPMCMC),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
Initialization = x;
Names = ParametersPMCMC.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
ParametersPMCMC = UpdateParsTransfToNoTransf(ParametersPMCMC);
TellParsValues(ParametersPMCMC)


Res = KalmanNumericDerivativesWithPrior(Data,HIVModel,ParametersPMCMC);

tempCov = (-Res.Hess)^-1;
tempCov(end+1,end+1) = (std(ResRW1_RW.Thetas))^2;

Parameters.G = tempCov^-1;
Parameters.RWsamplCov = tempCov^-1;
eig(tempCov^-1)>0

ParametersPMCMC.DiffusionType = 'Add';
ParametersPMCMC.InitialFt.Estimated = 0;
ParametersPMCMC.SecondFt.Estimated = 0;
ParametersPMCMC.ThirdFt.Estimated = 0;
ParametersPMCMC.FourthFt.Estimated = 0;
ParametersPMCMC.SigmaRW.Estimated = 1;
ParametersPMCMC = DefineEstimatedParametersIndexes(ParametersPMCMC);
ParametersPMCMC = DefineTransfFunctions(ParametersPMCMC);
ParametersPMCMC = UpdateParsNoTransfToTransf(ParametersPMCMC);
ParametersPMCMC = DefinePriors(ParametersPMCMC);


ParametersPMCMC.Epsil = 1;
ParametersPMCMC.ComputeRWsamplCov = 0;
ParametersPMCMC.NoPaths = 1;
ParametersPMCMC.MCMCType = 'Rand';
ParametersPMCMC.GMeth = 'cst given';
ParametersPMCMC.aim = 0.23;
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_RW_andPars_Calib = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_ResRW1_RW_andPars_Calib.mat'],'ResRW1_RW_andPars_Calib')


TempPar = ResRW1_RW_andPars_Calib.TempPar;
ParametersPMCMC.NoPaths = 0;
TransfThetasSamples = ResRW1_RW_andPars_Calib.TransfThetas;
ParametersPMCMC = DefineScalingPars(TransfThetasSamples,ParametersPMCMC);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,ParametersPMCMC);
ParametersPMCMC.GMeth = 'GMM';
ParametersPMCMC.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',ParametersPMCMC,1);
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_andPars_RW = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_ResRW1_RW_andPars.mat'],'ResRW1_RW_andPars')




clf
Names = Parameters.Names.Estimated;
k = length(Names);
for i = 1:k
    subplot(ceil(sqrt(k)),ceil(sqrt(k)),i)
%     plot(ResRWML.Thetas(i,:),ResRWML.Paths(:,8,end),'.')
    plot(ResRW1_RW_andPars_Calib.Thetas(i,:))
    hold on
    title(Names{i})
end



