function Res = HIVapplyInference(Data,Parameters)

HIVModel = struct();
HIVModel.EKF_projection = @HIV3_EKF_projection;
HIVModel.InitializeParameters = @HIV2_Initialize;
HIVModel.SMC_projection = @HIV3_SMC_projection;

try if Parameters.RealData

        Parameters.NbVariables = 9;
        Parameters.Problem = 'ImperialHIV';
        Parameters.ObservationLength = 25*12;
        Parameters.ComputationTStep = 0.5;
        Parameters.TypeWork = 'Normal';

        ObsVars = Parameters.ObsVars;
        ObsMin = Parameters.ObsMin;
        ObsMax = Parameters.ObsMax;
        ObsYears = Parameters.ObsYears;
        
        Data.Observations = zeros(10,length(ObsVars)+1);
        for i = 1:length(ObsVars)
            Data.Observations(ObsVars(i),1+i) = (ObsMin(i)+ObsMax(i))/2*100;
            Data.ObsSigmas(i+1) = ((ObsMax(i)-ObsMin(i))*100/4);
        end
        Instants = round((ObsYears-1985)*12);
        Data.Instants = round([0 Instants]/(Parameters.ComputationTStep));
        Data.ObservedVariables = [ 0 ObsVars];
        Data.NbComputingSteps = [0 diff(Data.Instants)];

        HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),(Parameters.ObsMax(IndTime-1)+Parameters.ObsMin(IndTime-1))*100/2,(Parameters.ObsMax(IndTime-1)-Parameters.ObsMin(IndTime-1))*100/4)';

        temp7 = zeros(1,10);
        temp8 = zeros(1,10);
        temp7(1,7) = 1;
        temp8(1,8) = 1;
        temps{7} = temp7;
        temps{8} = temp8;
        HIVModel.ObservationJacobian = {};
        for i = 1:length(ObsVars)
            HIVModel.ObservationJacobian{i+1} = temps{ObsVars(i)};
        end
        HIVModel.ObservationMeasurementNoise = {};
        for i = 1:length(ObsVars)
            HIVModel.ObservationMeasurementNoise{i+1} = ((Parameters.ObsMax(i)-Parameters.ObsMin(i))*100/4)^2;%(Data.Observations(ObsVars(i),i+1)*(100-Data.Observations(ObsVars(i),i+1))/400);
        end
    end
    
    NbItsPMCMC = 20000;
catch

    temp7 = zeros(1,10);
    temp8 = zeros(1,10);
    temp7(1,7) = 1;
    temp8(1,8) = 1;
    temps{7} = temp7;
    temps{8} = temp8;
    HIVModel.ObservationJacobian = {};
    ObsVars = Data.ObservedVariables;
    for i = 1:length(ObsVars)-1
        HIVModel.ObservationJacobian{i+1} = temps{ObsVars(i+1)};
    end
    HIVModel.ObservationMeasurementNoise = {};
    for i = 1:length(ObsVars)-1
    %     HIVModel.ObservationMeasurementNoise{i+1} = 0.01*Data.Observations(ObsVars(i+1),i+1);%sqrt(Data.Observations(ObsVars(i+1),i+1)*(100-Data.Observations(ObsVars(i+1),i+1))/400);
        HIVModel.ObservationMeasurementNoise{i+1} = sqrt(Data.Observations(ObsVars(i+1),i+1)*(100-Data.Observations(ObsVars(i+1),i+1))/400);
    end
    
    
    HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),Parameters.MultNoise*Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))';
    % HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';
    NbItsPMCMC = 3000;
end
    
Parameters.Problem = 'ImperialHIV2';
% Parameters.DiffusionType = 'Add';
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.Value = 0.05;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters = DefinePriors(Parameters);
% 
% Parameters.DiffusionType = 'OUD';
% Parameters.SigmaOU.Estimated = 1;
% Parameters.KappaOU.Estimated = 1;
% Parameters.MuOU.Estimated = 1;
% Parameters.SigmaOU.Value = 0.05;
% Parameters.KappaOU.Value = 0.0001;
% Parameters.MuOU.Value = 0.5;
% Parameters = DefineEstimatedParametersIndexes(Parameters);
% Parameters = DefineTransfFunctions(Parameters);
% Parameters = UpdateParsNoTransfToTransf(Parameters);
% Parameters = DefinePriors(Parameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('Initialize maxim. alg.')
% Names = Parameters.Names.All;
% for i = 1:length(Names);
%     Parameters.(Names{i}).Estimated = 0;
% end
% Parameters.InitialFt.Estimated = 1;
% Parameters.InitialIPropF.Estimated = 1;
% Parameters.InitialIPropM.Estimated = 1;
% Parameters.TotMFactor.Estimated = 1;
% Parameters.Alpham1.Estimated = 1;
% Parameters.MuFm1.Estimated = 1;
% Parameters.MuMm1.Estimated = 1;
% Parameters.BetaMFPerAct.Estimated = 1;
% Parameters.BetaFMPerAct.Estimated = 1;
% Parameters.NumberActsPerClient.Estimated = 1;
% Parameters.eHIV.Estimated = 1;
% Parameters.CF1.Estimated = 1;
% Parameters.CF2.Estimated = 1;
% Parameters.SigmaRW.Estimated = 1;
% Parameters = DefineEstimatedParametersIndexes(Parameters);
% Parameters = DefineTransfFunctions(Parameters);
% Parameters = UpdateParsNoTransfToTransf(Parameters);
% Parameters = DefinePriors(Parameters);
% Names = Parameters.Names.Estimated;
% 
% test = 0;
% for k = 1:length(Names)
%     Parameters.(Names{k}).Sample = 1;
% end
% NbIts = 0;
% Temp = EstimationEKFGen(Data, HIVModel, Parameters);
% while not(test)
%     Parameters = SampleParameters(Parameters);
%     try
%         Temp = EstimationEKFGen(Data, HIVModel, Parameters);
%         if Temp.LogLik>-1000
%             test = 1;
%         end
%     end
%     NbIts = NbIts+1;
%     if NbIts == 10000
%         test = 1;
%     end
% end

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('EKF Maxim. alg.')
% Parameters.Correction = 0;
% Initialization = [];
% for i = 1:length(Names)
%     Initialization(i) = Parameters.(Names{i}).TransfValue ;
% end
% Parameters.RWinEKF = 0;
% [x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,Parameters),Initialization,optimset('MaxIter',5000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',5000));
% Names = Parameters.Names.Estimated;
% for i = 1:length(Names)
%     Parameters.(Names{i}).TransfValue = (x(i));
% end
% Parameters = UpdateParsTransfToNoTransf(Parameters);
% TellParsValues(Parameters)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('EKF Maxim. alg. all pars')
% ParametersKalman = Parameters;
% Names = ParametersKalman.Names.All;
% for i = 1:length(Names);
%     ParametersKalman.(Names{i}).Estimated = 0;
% end
% ParametersKalman.InitialFt.Estimated = 1;
% % ParametersKalman.SecondFt.Estimated = 1;
% % ParametersKalman.ThirdFt.Estimated = 1;
% % ParametersKalman.FourthFt.Estimated = 1;
% % ParametersKalman.FifthFt.Estimated = 1;
% % ParametersKalman.SixthFt.Estimated = 1;
% ParametersKalman.InitialIPropF.Estimated = 1;
% ParametersKalman.InitialIPropM.Estimated = 1;
% ParametersKalman.TotMFactor.Estimated = 1;
% ParametersKalman.Alpham1.Estimated = 1;
% ParametersKalman.MuFm1.Estimated = 1;
% ParametersKalman.MuMm1.Estimated = 1;
% ParametersKalman.BetaMFPerAct.Estimated = 1;
% ParametersKalman.BetaFMPerAct.Estimated = 1;
% ParametersKalman.NumberActsPerClient.Estimated = 1;
% ParametersKalman.eHIV.Estimated = 1;
% ParametersKalman.CF1.Estimated = 1;
% ParametersKalman.CF2.Estimated = 1;
% if strcmp(Parameters.DiffusionType,'Add')
%     ParametersKalman.SigmaRW.Estimated = 1;
% elseif strcmp(Parameters.DiffusionType,'OUD')
%     ParametersKalman.SigmaOU.Estimated = 1;
%     ParametersKalman.KappaOU.Estimated = 1;
%     ParametersKalman.MuOU.Estimated = 1;
% end
% 
% ParametersKalman = DefineEstimatedParametersIndexes(ParametersKalman);
% ParametersKalman = DefineTransfFunctions(ParametersKalman);
% ParametersKalman = UpdateParsNoTransfToTransf(ParametersKalman);
% ParametersKalman = DefinePriors(ParametersKalman);
% Names = ParametersKalman.Names.Estimated;
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('EKF Maxim. alg.')
% ParametersKalman.Correction = 0;
% Initialization = [];
% for i = 1:length(Names)
%     Initialization(i) = ParametersKalman.(Names{i}).TransfValue ;
% end
% % ParametersKalman.RWinEKF = 1;
% ParametersKalman.Problem = 'ImperialHIV2';
% [x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,ParametersKalman),Initialization,optimset('MaxIter',5000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',2000));
% Names = ParametersKalman.Names.Estimated;
% for i = 1:length(Names)
%     ParametersKalman.(Names{i}).TransfValue = (x(i));
% end
% ParametersKalman = UpdateParsTransfToNoTransf(ParametersKalman);
% TellParsValues(ParametersKalman)
% 
% Parameters = ParametersKalman; % we take those because they maximize the posterior

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParametersKalman = Parameters;

disp('EKF Maxim. alg. all pars for Hess')
Names = ParametersKalman.Names.All;
for i = 1:length(Names);
    ParametersKalman.(Names{i}).Estimated = 0;
end
if strcmp(ParametersKalman.DiffusionType,'Logistic')
    % ParametersKalman.SecondFt.Estimated = 1;
    % ParametersKalman.ThirdFt.Estimated = 1;
    % ParametersKalman.FourthFt.Estimated = 1;
    % ParametersKalman.FifthFt.Estimated = 1;
    % ParametersKalman.SixthFt.Estimated = 1;
    ParametersKalman.InitialIPropF.Estimated = 1;
    ParametersKalman.InitialIPropM.Estimated = 1;
    ParametersKalman.TotMFactor.Estimated = 1;
    ParametersKalman.Alpham1.Estimated = 1;
    ParametersKalman.MuFm1.Estimated = 1;
    ParametersKalman.MuMm1.Estimated = 1;
    ParametersKalman.BetaMFPerAct.Estimated = 1;
    ParametersKalman.BetaFMPerAct.Estimated = 1;
    ParametersKalman.NumberActsPerClient.Estimated = 1;
    ParametersKalman.eHIV.Estimated = 1;
    ParametersKalman.CF1.Estimated = 1;
    ParametersKalman.CF2.Estimated = 1;
    ParametersKalman.CUinit.Estimated = 1;
    ParametersKalman.CUdelta.Estimated = 1;
    ParametersKalman.k.Estimated = 1;
    ParametersKalman.SigmaRW.Estimated = 1;
else
    ParametersKalman.InitialFt.Estimated = 1;
    % ParametersKalman.SecondFt.Estimated = 1;
    % ParametersKalman.ThirdFt.Estimated = 1;
    % ParametersKalman.FourthFt.Estimated = 1;
    % ParametersKalman.FifthFt.Estimated = 1;
    % ParametersKalman.SixthFt.Estimated = 1;
    ParametersKalman.InitialIPropF.Estimated = 1;
    ParametersKalman.InitialIPropM.Estimated = 1;
    ParametersKalman.TotMFactor.Estimated = 1;
    ParametersKalman.Alpham1.Estimated = 1;
    ParametersKalman.MuFm1.Estimated = 1;
    ParametersKalman.MuMm1.Estimated = 1;
    ParametersKalman.BetaMFPerAct.Estimated = 1;
    ParametersKalman.BetaFMPerAct.Estimated = 1;
    ParametersKalman.NumberActsPerClient.Estimated = 1;
    ParametersKalman.eHIV.Estimated = 1;
    ParametersKalman.CF1.Estimated = 1;
    ParametersKalman.CF2.Estimated = 1;
    ParametersKalman.SigmaRW.Estimated = 1;
end
% ParametersKalman.SigmaRW.Value = 0.01;
ParametersKalman = DefineEstimatedParametersIndexes(ParametersKalman);
ParametersKalman = DefineTransfFunctions(ParametersKalman);
ParametersKalman = UpdateParsNoTransfToTransf(ParametersKalman);
ParametersKalman = DefinePriors(ParametersKalman);
Names = ParametersKalman.Names.Estimated;


% ParametersKalman.SigmaRW.Value = 0.05;
% ParametersKalman = DefineEstimatedParametersIndexes(ParametersKalman);
% ParametersKalman = DefineTransfFunctions(ParametersKalman);
% ParametersKalman = UpdateParsNoTransfToTransf(ParametersKalman);
% ParametersKalman = DefinePriors(ParametersKalman);


disp('EKF Maxim. alg. for Hess')
ParametersKalman.Correction = 1;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersKalman.(Names{i}).TransfValue ;
end
% ParametersKalman.RWinEKF = 1;
% ParametersKalman.DiffusionType = 'Logistic';
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,ParametersKalman),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',2000));
Names = ParametersKalman.Names.Estimated;
for i = 1:length(Names)
    ParametersKalman.(Names{i}).TransfValue = (x(i));
end
ParametersKalman = UpdateParsTransfToNoTransf(ParametersKalman);
TellParsValues(ParametersKalman)

ParametersKalman.Finalx = x;
ResKal = KalmanNumericDerivativesWithPrior(Data,HIVModel,ParametersKalman);
Test = mean(eig(-ResKal.Hess)>0)==1;
Cov = (-ResKal.Hess)^-1;
Test


while not(Test)
    ParametersKalman.Correction = 1;
    Initialization = [];
    for i = 1:length(Names)
        Initialization(i) = ParametersKalman.(Names{i}).TransfValue ;
    end
    % ParametersKalman.RWinEKF = 1;
    ParametersKalman.DiffusionType = 'Add';
    [x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,ParametersKalman),Initialization,optimset('MaxIter',5000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',2000));
    Names = ParametersKalman.Names.Estimated;
    for i = 1:length(Names)
        ParametersKalman.(Names{i}).TransfValue = (x(i));
    end
    ParametersKalman = UpdateParsTransfToNoTransf(ParametersKalman);
    TellParsValues(ParametersKalman)

    ParametersKalman.Finalx = x;
    ResKal = KalmanNumericDerivativesWithPrior(Data,HIVModel,ParametersKalman);
    Test = mean(eig(-ResKal.Hess)>0)==1;
    Cov = (-ResKal.Hess)^-1;
    Test
end

Parameters = ParametersKalman;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('SMC Maxim. alg. only sigma')
Parameters.NbParticules = 1000;
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [7 8 9];
% try
%     if Parameters.SwitchToIBM == 1
%         Parameters.DiffusionType = 'IBM';
%         Parameters.SigmaRW.Value = 0.005;
%     else
%         Parameters.DiffusionType = 'Add';
%         Parameters.SigmaRW.Value = 0.05;
%     end
% catch    
%     Parameters.DiffusionType = 'Add';
%     Parameters.SigmaRW.Value = 0.15;
% end
Parameters.MCMCType = 'jiji';
Parameters.GMeth = 'cst given';
Names = Parameters.Names.All;
for i = 1:length(Names);
    Parameters.(Names{i}).Estimated = 0;
end


Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters = DefinePriors(Parameters);

% optimize RW parameters with SMC
Parameters.NoPaths = 0;
Names = Parameters.Names.Estimated;
Parameters.Correction = 1;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizewithPrior(x,Data,HIVModel,Parameters),Initialization,optimset('MaxIter',50,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)

% Result = EstimationEKFGen(Data, HIVModel, Parameters)
% ResultSMC = EstimationSMCsmoothGen(Data, HIVModel, Parameters)
% Result.LogLik
% ResultSMC.LogLik

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('SMC Maxim. alg. all pars')
% Names = Parameters.Names.All;
% for i = 1:length(Names);
%     Parameters.(Names{i}).Estimated = 0;
% end
% Parameters.InitialFt.Estimated = 1;
% Parameters.InitialIPropF.Estimated = 1;
% Parameters.InitialIPropM.Estimated = 1;
% Parameters.TotMFactor.Estimated = 1;
% Parameters.Alpham1.Estimated = 1;
% Parameters.MuFm1.Estimated = 1;
% Parameters.MuMm1.Estimated = 1;
% Parameters.BetaMFPerAct.Estimated = 1;
% Parameters.BetaFMPerAct.Estimated = 1;
% Parameters.NumberActsPerClient.Estimated = 1;
% Parameters.eHIV.Estimated = 1;
% Parameters.CF1.Estimated = 1;
% Parameters.CF2.Estimated = 1;
% Parameters.SigmaRW.Estimated = 1;
% Parameters = DefineEstimatedParametersIndexes(Parameters);
% Parameters = DefineTransfFunctions(Parameters);
% Parameters = UpdateParsNoTransfToTransf(Parameters);
% Parameters = DefinePriors(Parameters);
% Names = Parameters.Names.Estimated;
% 
% % optimize RW parameters with SMC
% Parameters.NoPaths = 0;
% Names = Parameters.Names.Estimated;
% Initialization = [];
% for i = 1:length(Names)
%     Initialization(i) = Parameters.(Names{i}).TransfValue ;
% end
% [x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizewithPrior(x,Data,HIVModel,Parameters),Initialization,optimset('MaxIter',80,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',5000));
% Names = Parameters.Names.Estimated;
% for i = 1:length(Names)
%     Parameters.(Names{i}).TransfValue = (x(i));
% end
% Parameters = UpdateParsTransfToNoTransf(Parameters);
% TellParsValues(Parameters)
% 
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Quick PMCMC')
for i = 1:length(Names);
    Parameters.(Names{i}).Estimated = 0;
end
Parameters.InitialFt.Estimated = 1;
Parameters.InitialIPropF.Estimated = 1;
Parameters.InitialIPropM.Estimated = 1;
Parameters.TotMFactor.Estimated = 1;
Parameters.Alpham1.Estimated = 1;
Parameters.MuFm1.Estimated = 1;
Parameters.MuMm1.Estimated = 1;
Parameters.BetaMFPerAct.Estimated = 1;
Parameters.BetaFMPerAct.Estimated = 1;
Parameters.NumberActsPerClient.Estimated = 1;
Parameters.eHIV.Estimated = 1;
Parameters.CF1.Estimated = 1;
Parameters.CF2.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters = DefinePriors(Parameters);
Names = Parameters.Names.Estimated;



Cov =  2.38^2/length(Names)*(-ResKal.Hess)^-1;%Parameters.CovInit;
Parameters.G = Cov^-1;
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [7 8 9];
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
try
    if Parameters.SwitchToIBM == 1
        Parameters.DiffusionType = 'IBM';
    end
catch    
    
end
Parameters.aim = 0.23;
Parameters.Epsil = 1;
Parameters.ModelType = 'SMC';
Parameters.Correction = 1.5;
TempPar = ProposeInitialParameter(Data, HIVModel, Parameters);
Parameters.KeepAll = 1;
Parameters.AdaptC = 0.999;
% [ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
Res = RunEstimationMethod(Data, HIVModel,Parameters,TempPar,2000);
Res.Parameters = Parameters;

try
    SavePath = 'S:\Results';
    save([SavePath Parameters.NameToSave '.mat'],'Res')
end

Cov =  cov(Res.TransfThetas');
Parameters.G = Cov^-1;
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [7 8 9];
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
try
    if Parameters.SwitchToIBM == 1
        Parameters.DiffusionType = 'IBM';
    end
catch    
    Parameters.DiffusionType = 'Add';
end
Parameters.aim = 0.23;
Parameters.Epsil = 1;
Parameters.ModelType = 'SMC';
Parameters.Correction = 1;
TempPar = ProposeInitialParameter(Data, HIVModel, Parameters);
Parameters.KeepAll = 1;
Parameters.AdaptC = 0.999;
% [ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
Res = RunEstimationMethod(Data, HIVModel,Parameters,TempPar,NbItsPMCMC);
Res.Parameters = Parameters;

% 


% 
% disp('Step 3: Sigma Opt')
% % generate samples from posterior
% Parameters.NbParticules = 5000;
% Parameters.DiffusionType = 'Add';
% Res = EstimationSMCfiltGen(Data, HIVModel, Parameters);
% Res.LogLik
% 
% NbIts = 100;
% Paths = zeros(NbIts,length(Parameters.PathsToKeep),sum(Data.NbComputingSteps));
% LogLiks = [];
% for i = 1:NbIts
%     disp(i)
%     Temp = EstimationSMCsmoothGen(Data, HIVModel, Parameters);
%     ind = ceil(rand(1,1)*Parameters.NbParticules);
%     Paths(i,:,:) = Temp.CompletePaths(ind,:,:);
%     LogLiks(i) = Temp.LogLik;
% end
% Res.Paths = Paths;