function Res = DengueapplyInference(Data,Parameters)

DengueModel = struct();
DengueModel.EKF_projection = @Dengue_EKF_projection;
DengueModel.InitializeParameters = @Dengue_Initialize;
DengueModel.SMC_projection = @Dengue_SMC_projection;



DengueModel.LikFunction = 'normpdf(Parameters.rep1.Value*Parameters.rep2.Value*Variables(:,11),Data.Observations(11,IndTime),sqrt(Parameters.rep2.Value*(1.0-Parameters.rep2.Value)*Parameters.rep1.Value*(Data.Observations(11,IndTime)) + (Parameters.rep2.Value*Parameters.phi.Value*Parameters.rep1.Value*(Data.Observations(11,IndTime)))^2))';
        % DengueModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';
    NbItsPMCMC = 6000;

    
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.Value = 0.05;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters = DefinePriors(Parameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParametersKalman = Parameters;

disp('EKF Maxim. alg. all pars for Hess')
Names = ParametersKalman.Names.All;
for i = 1:length(Names);
    ParametersKalman.(Names{i}).Estimated = 0;
end

ParametersKalman.SInit.Estimated = 1;
ParametersKalman.R1Init.Estimated = 1;
ParametersKalman.R2Init.Estimated = 1;
ParametersKalman.I1Init.Estimated = 1;
ParametersKalman.I2Init.Estimated = 1;
ParametersKalman.I12Init.Estimated = 1;
ParametersKalman.I21Init.Estimated = 1;
ParametersKalman.Q1Init.Estimated = 1;
ParametersKalman.Q2Init.Estimated = 1;
ParametersKalman.r0Init.Estimated = 1;
ParametersKalman.vm1.Estimated = 1;
ParametersKalman.qm1.Estimated = 1;
ParametersKalman.eta.Estimated = 1;
ParametersKalman.e.Estimated = 1;
ParametersKalman.d.Estimated = 1;
ParametersKalman.ade.Estimated = 1;
ParametersKalman.rep2.Estimated = 1;
ParametersKalman.phi.Estimated = 0;
ParametersKalman.SigmaRW.Estimated = 1;
ParametersKalman = DefineEstimatedParametersIndexes(ParametersKalman);
ParametersKalman = DefineTransfFunctions(ParametersKalman);
ParametersKalman = UpdateParsNoTransfToTransf(ParametersKalman);
ParametersKalman = DefinePriors(ParametersKalman);
Names = ParametersKalman.Names.Estimated;

disp('EKF Maxim. alg. for Hess')
ParametersKalman.Correction = 1;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersKalman.(Names{i}).TransfValue ;
end
% ParametersKalman.RWinEKF = 1;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,DengueModel,ParametersKalman),Initialization,optimset('MaxIter',5000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',2000));
Names = ParametersKalman.Names.Estimated;
for i = 1:length(Names)
    ParametersKalman.(Names{i}).TransfValue = (x(i));
end
ParametersKalman = UpdateParsTransfToNoTransf(ParametersKalman);
TellParsValues(ParametersKalman)

ParametersKalman.Finalx = x;
ResKal = KalmanNumericDerivativesWithPrior(Data,DengueModel,ParametersKalman);
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
    [x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,DengueModel,ParametersKalman),Initialization,optimset('MaxIter',5000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',2000));
    Names = ParametersKalman.Names.Estimated;
    for i = 1:length(Names)
        ParametersKalman.(Names{i}).TransfValue = (x(i));
    end
    ParametersKalman = UpdateParsTransfToNoTransf(ParametersKalman);
    TellParsValues(ParametersKalman)

    ParametersKalman.Finalx = x;
    ResKal = KalmanNumericDerivativesWithPrior(Data,DengueModel,ParametersKalman);
    Test = mean(eig(-ResKal.Hess)>0)==1;
    Cov = (-ResKal.Hess)^-1;
    Test
end

Parameters = ParametersKalman;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('SMC Maxim. alg. only sigma')
Parameters.NbParticules = 1000;
Parameters.StableCUseConstraint = 0;
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
    Parameters.SigmaRW.Value = 0.001;
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
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizewithPrior(x,Data,DengueModel,Parameters),Initialization,optimset('MaxIter',50,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)

% Result = EstimationEKFGen(Data, DengueModel, Parameters)
% ResultSMC = EstimationSMCsmoothGen(Data, DengueModel, Parameters)
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
% Parameters.eDengue.Estimated = 1;
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
% [x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizewithPrior(x,Data,DengueModel,Parameters),Initialization,optimset('MaxIter',80,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',5000));
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
Names = Parameters.Names.Estimated;
for i = 1:length(Names);
    Parameters.(Names{i}).Estimated = 0;
end
if strcmp(Parameters.DiffusionType,'Bertallanfy')
    Parameters.InitialIPropF.Estimated = 1;
    Parameters.InitialIPropM.Estimated = 1;
    Parameters.TotMFactor.Estimated = 1;
    Parameters.Alpham1.Estimated = 1;
    Parameters.MuFm1.Estimated = 1;
    Parameters.MuMm1.Estimated = 1;
    Parameters.BetaMFPerAct.Estimated = 1;
    Parameters.BetaFMPerAct.Estimated = 1;
    Parameters.NumberActsPerClient.Estimated = 1;
    Parameters.eDengue.Estimated = 1;
    Parameters.CF1.Estimated = 1;
    Parameters.CF2.Estimated = 1;
    Parameters.BRbase.Estimated = 1;
    Parameters.BRmu.Estimated = 1;
    Parameters.BRtinfl.Estimated = 1;
    Parameters.BRm.Estimated = 1;
    Parameters.SigmaRW.Estimated = 1;
else
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
    Parameters.eDengue.Estimated = 1;
    Parameters.CF1.Estimated = 1;
    Parameters.CF2.Estimated = 1;
    Parameters.SigmaRW.Estimated = 1;
end
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
% try
%     if Parameters.SwitchToIBM == 1
%         Parameters.DiffusionType = 'IBM';
%     end
% catch    
%     Parameters.DiffusionType = 'Add';
% end
Parameters.aim = 0.23;
Parameters.Epsil = 1;
Parameters.ModelType = 'SMC';
Parameters.Correction = 1.5;
TempPar = ProposeInitialParameter(Data, DengueModel, Parameters);
Parameters.KeepAll = 1;
Parameters.AdaptC = 0.999;
% [ParametersPMCMC, TempPar] = CalibrateMethod( Data, DengueModel, ParametersPMCMC, TempPar);
Res = RunEstimationMethod(Data, DengueModel,Parameters,TempPar,2000);
Res.Parameters = Parameters;

Cov =  cov(Res.TransfThetas');
Parameters.G = Cov^-1;
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [7 8 9];
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
% try
%     if Parameters.SwitchToIBM == 1
%         Parameters.DiffusionType = 'IBM';
%     end
% catch    
%     Parameters.DiffusionType = 'Add';
% end
Parameters.aim = 0.23;
Parameters.Epsil = 1;
Parameters.ModelType = 'SMC';
Parameters.Correction = 1.5;
TempPar = ProposeInitialParameter(Data, DengueModel, Parameters);
Parameters.KeepAll = 1;
Parameters.AdaptC = 0.99;
% [ParametersPMCMC, TempPar] = CalibrateMethod( Data, DengueModel, ParametersPMCMC, TempPar);
Res = RunEstimationMethod(Data, DengueModel,Parameters,TempPar,2000);
Res.Parameters = Parameters;


Cov =  cov(Res.TransfThetas');
Parameters.G = Cov^-1;
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [7 8 9];
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
% try
%     if Parameters.SwitchToIBM == 1
%         Parameters.DiffusionType = 'IBM';
%     end
% catch    
%     Parameters.DiffusionType = 'Add';
% end
Parameters.aim = 0.23;
Parameters.Epsil = 1;
Parameters.ModelType = 'SMC';
Parameters.Correction = 1;
TempPar = ProposeInitialParameter(Data, DengueModel, Parameters);
Parameters.KeepAll = 1;
Parameters.AdaptC = 0.999;
% [ParametersPMCMC, TempPar] = CalibrateMethod( Data, DengueModel, ParametersPMCMC, TempPar);
Res = RunEstimationMethod(Data, DengueModel,Parameters,TempPar,NbItsPMCMC);
Res.Parameters = Parameters;

try
    SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';
    save([SavePath Parameters.NameToSave '.mat'],'Res')
    disp('saved')
end


% 


% 
% disp('Step 3: Sigma Opt')
% % generate samples from posterior
% Parameters.NbParticules = 5000;
% Parameters.DiffusionType = 'Add';
% Res = EstimationSMCfiltGen(Data, DengueModel, Parameters);
% Res.LogLik
% 
% NbIts = 100;
% Paths = zeros(NbIts,length(Parameters.PathsToKeep),sum(Data.NbComputingSteps));
% LogLiks = [];
% for i = 1:NbIts
%     disp(i)
%     Temp = EstimationSMCsmoothGen(Data, DengueModel, Parameters);
%     ind = ceil(rand(1,1)*Parameters.NbParticules);
%     Paths(i,:,:) = Temp.CompletePaths(ind,:,:);
%     LogLiks(i) = Temp.LogLik;
% end
% Res.Paths = Paths;