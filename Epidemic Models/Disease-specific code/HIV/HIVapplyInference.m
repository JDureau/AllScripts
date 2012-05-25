function Res = HIVapplyInference(Data,Parameters)

HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIV_Initialize;
HIVModel.SMC_projection = @HIV_SMC_projection;
SavePath = '/users/ecologie/dureau/src/AllData/Avahan';

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

        if or(strcmp(Parameters.DiffusionType,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
%             HIVModel.LikFunction = 'not(Res.Crash)*Res.WentOutOrNot.*normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),(Parameters.ObsMax(IndTime-1)+Parameters.ObsMin(IndTime-1))*100/2,(Parameters.ObsMax(IndTime-1)-Parameters.ObsMin(IndTime-1))*100/4)';
            HIVModel.LikFunction = 'not(Res.Crash)*Res.WentOutOrNot.*normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),(Parameters.ObsMax(IndTime-1)+Parameters.ObsMin(IndTime-1))*100/2,sqrt(Parameters.Obs(IndTime-1)*100*(100-Parameters.Obs(IndTime-1)*100)/400))';
        else
%             HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),(Parameters.ObsMax(IndTime-1)+Parameters.ObsMin(IndTime-1))*100/2,(Parameters.ObsMax(IndTime-1)-Parameters.ObsMin(IndTime-1))*100/4)';
            HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),(Parameters.ObsMax(IndTime-1)+Parameters.ObsMin(IndTime-1))*100/2,sqrt(Parameters.Obs(IndTime-1)*100*(100-Parameters.Obs(IndTime-1)*100)/400))';
        end

        temp7 = zeros(1,9);
        temp8 = zeros(1,9);
        temp7(1,7) = 1;
        temp8(1,8) = 1;
        temps{7} = temp7;
        temps{8} = temp8;
        HIVModel.ObservationJacobian = {};
        for i = 1:length(ObsVars)
            HIVModel.ObservationJacobian{i+1} = temps{ObsVars(i)};
        end
        HIVModel.ObservationMeasurementNoise = {};
        Obs = Parameters.Obs*100;
        for i = 1:length(ObsVars)
%             HIVModel.ObservationMeasurementNoise{i+1} = ((Parameters.ObsMax(i)-Parameters.ObsMin(i))*100/4)^2;%(Data.Observations(ObsVars(i),i+1)*(100-Data.Observations(ObsVars(i),i+1))/400);
            HIVModel.ObservationMeasurementNoise{i+1} = (Obs(i)*(100-Obs(i))/400);
       end
        NbItsPMCMC = 150000;
        Parameters.TempName = ['Temp_' Parameters.NameToSave '_' Parameters.DiffusionType '.mat'];

    end
    
    
catch

    %%% CORRECT LIKE IN REAL DATA
    Parameters.NbVariables = 9;
    temp7 = zeros(1,9);
    temp8 = zeros(1,9);
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
        HIVModel.ObservationMeasurementNoise{i+1} = (Data.Observations(ObsVars(i+1),i+1)*(100-Data.Observations(ObsVars(i+1),i+1))/400); % error here, shouldn't be sqrt...
    end
    
    if or(strcmp(Parameters.DiffusionType,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
        HIVModel.LikFunction = 'not(Res.Crash)*normpdf(Res.WentOutOrNot.*Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)/400)))';
    else
        HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)/400)))';
    end
        % HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';
    NbItsPMCMC = 30000;
end

    
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.Problem = 'ImperialHIV2';
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.Value = 0.1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters = DefinePriors(Parameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  disp('Initialize maxim. alg.')
%  Names = Parameters.Names.All;
%  for i = 1:length(Names);
%      Parameters.(Names{i}).Estimated = 0;
%  end
%  if strcmp(Parameters.DiffusionType,'Logistic')
%     Parameters.InitialIPropF.Estimated = 1;
%     Parameters.InitialIPropM.Estimated = 1;
%     Parameters.TotMFactor.Estimated = 1;
%     Parameters.Alpham1.Estimated = 1;
%     Parameters.MuFm1.Estimated = 1;
%     Parameters.MuMm1.Estimated = 1;
%     Parameters.BetaMFPerAct.Estimated = 1;
%     Parameters.BetaFMPerAct.Estimated = 1;
%     Parameters.NumberActsPerClient.Estimated = 1;
%     Parameters.eHIV.Estimated = 1;
%     Parameters.CF1.Estimated = 1;
%     Parameters.CF2.Estimated = 1;
%     Parameters.CUinit.Estimated = 1;
%     Parameters.CUdelta.Estimated = 1;
%     Parameters.k.Estimated = 1;
%     Parameters.SigmaRW.Estimated = 1;
% else
%     Parameters.InitialFt.Estimated = 1;
%     Parameters.InitialIPropF.Estimated = 1;
%     Parameters.InitialIPropM.Estimated = 1;
%     Parameters.TotMFactor.Estimated = 1;
%     Parameters.Alpham1.Estimated = 1;
%     Parameters.MuFm1.Estimated = 1;
%     Parameters.MuMm1.Estimated = 1;
%     Parameters.BetaMFPerAct.Estimated = 1;
%     Parameters.BetaFMPerAct.Estimated = 1;
%     Parameters.NumberActsPerClient.Estimated = 1;
%     Parameters.eHIV.Estimated = 1;
%     Parameters.CF1.Estimated = 1;
%     Parameters.CF2.Estimated = 1;
%     Parameters.SigmaRW.Estimated = 1;
% end
% Parameters = DefineEstimatedParametersIndexes(Parameters);
% Parameters = DefineTransfFunctions(Parameters);
% Parameters = UpdateParsNoTransfToTransf(Parameters);
% Parameters = DefinePriors(Parameters);
% 
%  test = 0;
%  for k = 1:length(Names)
%      Parameters.(Names{k}).Sample = 1;
%  end
%  NbIts = 0;
%  Temp = EstimationEKFGen(Data, HIVModel, Parameters);
%  while not(test)
%      Parameters = SampleParameters(Parameters);
%      try
% %          TellParsValues(Parameters);
%          Temp = EstimationEKFGen(Data, HIVModel, Parameters);
%          if Temp.LogLik>-1000
%              test = 1;
%          end
%      end
%      NbIts = NbIts+1;
%      if NbIts == 10000
%          test = 1;
%      end
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParametersKalman = Parameters;

disp('EKF Maxim. alg. all pars for Hess')
Names = ParametersKalman.Names.All;
for i = 1:length(Names);
    ParametersKalman.(Names{i}).Estimated = 0;
end
if strcmp(ParametersKalman.DiffusionType,'Bertallanfy')
 
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
    ParametersKalman.BRbase.Estimated = 1;
    ParametersKalman.BRmu.Estimated = 1;
    ParametersKalman.BRtinfl.Estimated = 1;
    ParametersKalman.BRmm1.Estimated = 1;
    ParametersKalman.BRsigma.Estimated = 1;
elseif strcmp(ParametersKalman.DiffusionType,'BertallanfyConstr')
    ParametersKalman.BRmm1.Estimated = 1;
    ParametersKalman.BRbase.Estimated = 1;
    ParametersKalman.BRmu.Estimated = 1;
    ParametersKalman.BRtinfl.Estimated = 1;
    ParametersKalman.SigmaRW.Estimated = 1;
elseif strcmp(ParametersKalman.DiffusionType,'Add')
    ParametersKalman.InitialFt.Estimated = 1;
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
elseif strcmp(ParametersKalman.DiffusionType,'AddConstr')
    ParametersKalman.InitialFt.Estimated = 1;
    ParametersKalman.SigmaRW.Estimated = 1;
end
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
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,ParametersKalman),Initialization,optimset('MaxIter',2000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',4000));
Names = ParametersKalman.Names.Estimated;
for i = 1:length(Names)
    ParametersKalman.(Names{i}).TransfValue = (x(i));
end
ParametersKalman = UpdateParsTransfToNoTransf(ParametersKalman);
ParametersKalman = HIVModel.InitializeParameters(ParametersKalman);
TellParsValues(ParametersKalman)

ParametersKalman.Finalx = x;
ResKal = KalmanNumericDerivativesWithPrior(Data,HIVModel,ParametersKalman);
try
    Test = mean(eig(-ResKal.Hess)>0)==1;
    Test = Test*sum(sum(isreal(ResKal.Hess)));

    Cov = (-ResKal.Hess)^-1;
catch
    Test = 0;
end
Test

cpt = 1;
while not(Test)
    ParametersKalman.Correction = 1;
    Initialization = [];
    for i = 1:length(Names)
        Initialization(i) = ParametersKalman.(Names{i}).TransfValue ;
    end
    % ParametersKalman.RWinEKF = 1;
    [x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,ParametersKalman),Initialization,optimset('MaxIter',5000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
    Names = ParametersKalman.Names.Estimated;
    for i = 1:length(Names)
        ParametersKalman.(Names{i}).TransfValue = (x(i));
    end
    ParametersKalman = UpdateParsTransfToNoTransf(ParametersKalman);
    TellParsValues(ParametersKalman)

    ParametersKalman.Finalx = x;
    ResKal = KalmanNumericDerivativesWithPrior(Data,HIVModel,ParametersKalman);
    try
        Test = mean(eig(-ResKal.Hess)>0)==1;
        Test = Test*sum(sum(isreal(ResKal.Hess)));
        Cov = (-ResKal.Hess)^-1;
    end
    Test
    cpt = cpt+1;
    if cpt >10
        Test = 1;
        die
    end
end

ParametersKalman.KalCov = Cov;


Temp = struct();
Temp.ParametersKalman = ParametersKalman;
save([SavePath '/' Parameters.TempName],'Temp')

%  die


SavePath = '/users/ecologie/dureau/src/AllData/Avahan';
% SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan/Temp'

load([SavePath '/' Parameters.TempName])

Parameters = Temp.ParametersKalman;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('SMC Maxim. alg. only sigma')
% Parameters.NbParticules = 1000;
% Parameters.StableCUseConstraint = 0;
% Parameters.NoPaths = 0;
% Parameters.PathsToKeep = [7 8 9];
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Quick PMCMC')
Names = Parameters.Names.Estimated;
for i = 1:length(Names);
    Parameters.(Names{i}).Estimated = 0;
end
if strcmp(Parameters.DiffusionType,'Bertallanfy')
    Parameters.BRmm1.Estimated = 1;    
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
    Parameters.BRbase.Estimated = 1;
    Parameters.BRmu.Estimated = 1;
    Parameters.BRtinfl.Estimated = 1;
    Parameters.BRsigma.Estimated = 1;
    
  
    
elseif strcmp(Parameters.DiffusionType,'BertallanfyConstr')
    Parameters.BRmm1.Estimated = 1;
    Parameters.BRbase.Estimated = 1;
    Parameters.BRmu.Estimated = 1;
    Parameters.BRtinfl.Estimated = 1;
    Parameters.SigmaRW.Estimated = 1;
elseif strcmp(Parameters.DiffusionType,'Add')
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
elseif strcmp(Parameters.DiffusionType,'AddConstr')
    Parameters.InitialFt.Estimated = 1;
    Parameters.SigmaRW.Estimated = 1;
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters = DefinePriors(Parameters);
Names = Parameters.Names.Estimated;

% SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';
% if strcmp(Parameters.DiffusionType,'Bertallanfy')
%     load([SavePath '/BestBertCov.mat'])
% elseif  strcmp(Parameters.DiffusionType,'BertallanfyConstr')
%     load([SavePath '/BestBertConstrCov.mat'])
% else
%     Cov =  2.38^2/length(Names)*(-ResKal.Hess)^-1;%Parameters.CovInit;
% %     die
% end
Cov =  2.38^2/length(Names)*Parameters.KalCov;

disp('First MCMC')
Parameters.NbParticules  = 1000;
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
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
TempPar = ProposeInitialParameter(Data, HIVModel, Parameters);
% Parameters.KeepAll = 1;
Parameters.AdaptC = 0.99;
Parameters.AdMet = 0;
Parameters.AdMetBeta = 0.05;
% [ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
Res = RunEstimationMethod(Data, HIVModel,Parameters,TempPar,3000);
Res.Parameters = Parameters;

for i = 1:5
    disp(['MCMC ' num2str(i)])
    if Res.AccRate>0.05
        Cov =  2.38^2/length(Names)*cov(Res.TransfThetas');
        Parameters.G = Cov^-1;
    end
    Parameters.NoPaths = 1;
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
    TempPar = Res.TempPar;
%     Parameters.Kee��pAll = 1;
    Parameters.AdaptC = 0.98;
    % [ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
    Res = RunEstimationMethod(Data, HIVModel,Parameters,TempPar,3000);
    Res.Parameters = Parameters;
end

load([SavePath '/' Parameters.TempName])
Temp.ParametersPMCMC = Parameters;
Temp.ResPMCMCNoPaths = Res;
save([SavePath '/' Parameters.TempName],'Temp')


% die

load([SavePath '/' Parameters.TempName])
Res = Temp.ResPMCMCNoPaths;
Parameters = Res.Parameters;
Names = Parameters.Names.Estimated;

% if 1%try 
% 
%     load([Parameters.NameToSave])
%     Res = Ress{1};
%     Parameters = Res.Parameters;
%     Names = Parameters.Names.Estimated;
% end

disp(['Final MCMC'])

Cov =  2.38^2/length(Names)*cov(Res.TransfThetas');
% OldCov = Cov;
% Cov =  2.38^2/length(Names)*DetCov;
% Cov(Parameters.SigmaRW.Index,Parameters.SigmaRW.Index) = OldCov(13,13);
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
TempPar = Res.TempPar;
% Parameters.KeepAll = 1;
Parameters.AdaptC = 0.99;
% [ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
Parameters.SaveSpace = 1;
Res = RunEstimationMethod(Data, HIVModel,Parameters,TempPar,NbItsPMCMC);
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