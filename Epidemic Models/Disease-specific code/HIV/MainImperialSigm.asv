%% HIV Imperial Tests With Pars: Sigmoid

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

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
load([SavePath '\HIV_Mysore_ResRW1_RW_andPars.mat'])


Parameters = ResRW1_RW_Rds_1_2_3.Parameters;
Data = ResRW1_RW_Rds_1_2_3.Data;

%%% Model
HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIV_Initialize;
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*Parameters.SigmaObs).*(Res.WentOutOrNot)';
HIVModel.SMC_projection = @HIV_SMC_projection;


t = 0:1:Data.Instants(end);
FtPath = 0.05 + 0.8*sigmoid((t-3*t(end)/4)/50); %Test Mysore 1plot(t,FtPath)
plot(FtPath)

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
% Data.Instants = round([0 100 200 250 275]/(Parameters.ComputationTStep)); %Mysore5 Obs 1
Data.Instants = round([0 200 240 260 280]/(Parameters.ComputationTStep)); %MysoreSigm Obs2
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
Names = Parameters.Names.All;
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
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
[Parameters, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
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
ParametersPMCMC.Epsil = 1.2;
ResRW1_RW_Calib_Rds_1_2_3 = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

SavePath = 'C:\Users\dureau\Desktop\Results';
% SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_Tests_MysoreSigm_Obs2_withPars.mat'],'ResRW1_RW_Calib_Rds_1_2_3')


TempPar = ResRW1_RW_Calib_Rds_1_2_3.TempPar;
ParametersPMCMC.NoPaths = 0;
ParametersPMCMC.PathsToKeep = [7 8 9];
TransfThetasSamples = ResRW1_RW_Calib_Rds_1_2_3.TransfThetas;
ParametersPMCMC = DefineScalingPars(TransfThetasSamples,ParametersPMCMC);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,ParametersPMCMC);
ParametersPMCMC.GMeth = 'GMM';
ParametersPMCMC.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',ParametersPMCMC,1);
TempPar = ProposeInitialParameter(Data, HIVModel, ParametersPMCMC);
ParametersPMCMC.Epsil = 1;
[ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
ResRW1_RW_Rds_1_2 = RunEstimationMethod(Data, HIVModel,ParametersPMCMC,TempPar,5000);

ResRW1_RW_Rds_1_2.FtPath = FtPath;

SavePath = 'C:\Users\dureau\Desktop\Results';
% SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\HIV_Tests_MysoreSig_Obs2_withPars.mat'],'ResRW1_RW_Rds_1_2')

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



