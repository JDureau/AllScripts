function Res = HIVapplyInference(Data,Parameters)

HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIV_Initialize;
HIVModel.SMC_projection = @HIV_SMC_projection;
SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';

if 1%try 
    if Parameters.RealData

        Parameters.PathsToKeep = [7 8 9];

        Parameters.NbVariables = 9;
        Parameters.Problem = 'ImperialHIV';
        Parameters.ObservationLength = 25*12;
        Parameters.ComputationTStep = 0.5;
        Parameters.TypeWork = 'Normal';

        
        tmp = [];
        for i = 1:length(Parameters.ObsVars)
            tmp2 = Parameters.ObsVars{i};
            tmp(i) = tmp2(1);
        end
        indsclients = find(tmp == 8);
        indsFSWs    = find(tmp == 7);
        if Parameters.TakeClients                      
            inds = sort([indsclients(1:min(length(indsclients),Parameters.NbRounds)) indsFSWs(1:min(length(indsFSWs),Parameters.NbRounds))]);
        else
            inds = indsFSWs(1:min(length(indsFSWs),Parameters.NbRounds));
        end
        inds

        Obs = Parameters.Obs(inds);
        ObsVars = Parameters.ObsVars(inds);
        ObsYears = Parameters.ObsYears(inds);
        
        Data.Observations = zeros(10,length(ObsVars)+1);
        for i = 1:length(ObsVars)
            for j = 1:length(ObsVars{i})
                Data.Observations(ObsVars{i}(j),1+i) = (Obs{i}(j))*100;
            end
        end
        if and(Parameters.NbRounds == 3,max(length(indsclients),length(indsFSWs))<3)
            die
        end
        if and(Parameters.NbRounds == 2,max(length(indsclients),length(indsFSWs))<2)
            die
        end
        
        Instants = round((ObsYears-1985)*12);
        Data.Instants = round([0 Instants]/(Parameters.ComputationTStep));
        Data.ObservedVariables = [ 0 ObsVars];
        Data.NbComputingSteps = [0 diff(Data.Instants)];

        if or(strcmp(Parameters.DiffusionType,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
%             HIVModel.LikFunction = 'not(Res.Crash)*Res.WentOutOrNot.*normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),(Parameters.ObsMax(IndTime-1)+Parameters.ObsMin(IndTime-1))*100/2,(Parameters.ObsMax(IndTime-1)-Parameters.ObsMin(IndTime-1))*100/4)';
            HIVModel.LikFunction = 'not(Res.Crash)*Res.WentOutOrNot.*binopdf(round(Parameters.NbSamples(IndTime-1)*Parameters.Obs(IndTime-1)),Parameters.NbSamples(IndTime-1),Variables(:,Data.ObservedVariables(:,IndTime))/100)'; % On Oct 11, corrected here with Nb Samples (before there was 425 and 400... watch out)
        else
%             HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),(Parameters.ObsMax(IndTime-1)+Parameters.ObsMin(IndTime-1))*100/2,(Parameters.ObsMax(IndTime-1)-Parameters.ObsMin(IndTime-1))*100/4)';
%             HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),(Parameters.ObsMax(IndTime-1)+Parameters.ObsMin(IndTime-1))*100/2,sqrt(Parameters.Obs(IndTime-1)*100*(100-Parameters.Obs(IndTime-1)*100)/400))';
            HIVModel.LikFunction = 'binopdf(round(Parameters.NbSamples(IndTime-1)*Parameters.Obs{IndTime-1}),Parameters.NbSamples(IndTime-1),Variables(:,Data.ObservedVariables{IndTime})/100)';
        end

        temp7 = zeros(1,9);
%         temp7 = zeros(1,9);
        temp8 = zeros(1,9);
        temp7(1,7) = 1;
        temp8(1,8) = 1;
        temps{7} = temp7;
        temps{8} = temp8;
        HIVModel.ObservationJacobian = {};
        for i = 1:length(ObsVars)
            HIVModel.ObservationJacobian{i+1} = temps{ObsVars{i}(1)};
        end
        HIVModel.ObservationMeasurementNoise = {};
        for i = 1:length(ObsVars)
            for j = 1:length(ObsVars{i})
%             HIVModel.ObservationMeasurementNoise{i+1} = ((Parameters.ObsMax(i)-Parameters.ObsMin(i))*100/4)^2;%(Data.Observations(ObsVars(i),i+1)*(100-Data.Observations(ObsVars(i),i+1))/400);
                HIVModel.ObservationMeasurementNoise{i+1}(j) = (Parameters.Obs{i}(j)*100*(100-Parameters.Obs{i}(j)*100)/Parameters.NbSamples(i));
            end
       end
        NbItsPMCMC = 100000;
        Parameters.TempName = ['Temp_' Parameters.NameToSave '_' Parameters.DiffusionType '_' num2str(Parameters.NbRounds) '_' num2str(Parameters.TakeClients) '.mat'];

    else
        die
    end
    
    
else
    Parameters.NbVariables = 9;
    Parameters.Problem = 'ImperialHIV';
    Parameters.ObservationLength = 25*12;
    Parameters.ComputationTStep = 0.5;
    Parameters.TypeWork = 'Normal';

    Parameters.PathsToKeep = [9];

    
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
        HIVModel.LikFunction1 = 'not(Res.Crash)*binopdf(round(425*Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)/100),425,Variables(:,Data.ObservedVariables(:,IndTime))/100)';
        HIVModel.LikFunction2 = 'not(Res.Crash)*binopdf(round(425*Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)/100),425,Variables(:,Data.ObservedVariables(:,IndTime))/100)';
    else
%         HIVModel.LikFunction = 'not(Res.Crash)*normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)/400)))';
        HIVModel.LikFunction = 'binopdf(round(425*Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)/100),425,Variables(:,Data.ObservedVariables(:,IndTime))/100)';
    end
        % HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';
    NbItsPMCMC = 60000;
end

    
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.Problem = 'ImperialHIV2';
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.Value = 0.1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters = DefinePriors(Parameters);


% load([SavePath '/' Parameters.TempName])
try
    load([SavePath '/' Parameters.TempName])
    AlreadySomething = 1;
    Parameters = Temp.Parameters;
catch
    AlreadySomething = 0;
end
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

if not(AlreadySomething)
    ParametersKalman = Parameters;

    disp('EKF Maxim. alg. all pars for Hess')
    Names = ParametersKalman.Names.All;
    for i = 1:length(Names);
        ParametersKalman.(Names{i}).Estimated = 0;
    end
    if strcmp(ParametersKalman.DiffusionType,'Bertallanfy')

        ParametersKalman.InitialIPropF.Estimated = 1;
        ParametersKalman.InitialIPropM.Estimated = 1;
        ParametersKalman.Alpham1.Estimated = 1;
        ParametersKalman.MuFm1.Estimated = 1;
        ParametersKalman.MuMm1.Estimated = 1;
        ParametersKalman.BetaMFPerAct.Estimated = 1;
        ParametersKalman.BetaFMPerAct.Estimated = 1;
        ParametersKalman.NumberActsPerClient.Estimated = 1;
        ParametersKalman.eHIV.Estimated = 1;
        ParametersKalman.CF1.Estimated = 1;
        ParametersKalman.CF2.Estimated = 1;
        ParametersKalman.CM.Estimated = 1;
        ParametersKalman.BRbase.Estimated = 1;
        ParametersKalman.BRmu.Estimated = 1;
        ParametersKalman.BRtinfl.Estimated = 1;
        ParametersKalman.BRmm1.Estimated = 1;
        ParametersKalman.BRsigma.Value = 0.1;
        ParametersKalman.BRsigma.Min = -10^14;
        ParametersKalman.BRsigma.Max = 10^14;
        ParametersKalman.BRsigma.MinLim = 0;
        ParametersKalman.BRsigma.MaxLim = 2;
        ParametersKalman.BRsigma.Estimated = 1;
        ParametersKalman.BRsigma.TransfType = 'Logit';
        ParametersKalman.BRsigma.Init = 0;
    %     ParametersKalman.BRmm1.Value = 3;
    %     ParametersKalman.BRmu.Value = 0.8;
    elseif strcmp(ParametersKalman.DiffusionType,'Sigmoid')

        ParametersKalman.InitialIPropF.Estimated = 1;
        ParametersKalman.InitialIPropM.Estimated = 1;
        ParametersKalman.Alpham1.Estimated = 1;
        ParametersKalman.MuFm1.Estimated = 1;
        ParametersKalman.MuMm1.Estimated = 1;
        ParametersKalman.BetaMFPerAct.Estimated = 1;
        ParametersKalman.BetaFMPerAct.Estimated = 1;
        ParametersKalman.NumberActsPerClient.Estimated = 1;
        ParametersKalman.eHIV.Estimated = 1;
        ParametersKalman.CF1.Estimated = 1;
        ParametersKalman.CF2.Estimated = 1;
        ParametersKalman.CM.Estimated = 1;
        ParametersKalman.Sigmbase.Estimated = 1;
        ParametersKalman.Sigmmu.Estimated = 1;
        ParametersKalman.Sigmtinfl.Estimated = 1;
        ParametersKalman.Sigmrate.Estimated = 1;
        ParametersKalman.Sigmsigma.Value = 0.1;
        ParametersKalman.Sigmsigma.Min = -10^14;
        ParametersKalman.Sigmsigma.Max = 10^14;
        ParametersKalman.Sigmsigma.MinLim = 0;
        ParametersKalman.Sigmsigma.MaxLim = 2;
        ParametersKalman.Sigmsigma.Estimated = 1;
        ParametersKalman.Sigmsigma.TransfType = 'Logit';
        ParametersKalman.Sigmsigma.Init = 0;
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
        ParametersKalman.Alpham1.Estimated = 1;
        ParametersKalman.MuFm1.Estimated = 1;
        ParametersKalman.MuMm1.Estimated = 1;
        ParametersKalman.BetaMFPerAct.Estimated = 1;
        ParametersKalman.BetaFMPerAct.Estimated = 1;
        ParametersKalman.NumberActsPerClient.Estimated = 1;
        ParametersKalman.eHIV.Estimated = 1;
        ParametersKalman.CF1.Estimated = 1;
        ParametersKalman.CF2.Estimated = 1;
        ParametersKalman.CM.Estimated = 1;
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

    ParametersKalmanBeforeOpt = ParametersKalman;
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



    ParametersKalman = ParametersKalmanBeforeOpt;
    ParametersKalman.Correction = 0;
    Initialization = [];
    for i = 1:length(Names)
        Initialization(i) = ParametersKalman.(Names{i}).TransfValue ;
    end
    % ParametersKalman.RWinEKF = 1;
    [x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,ParametersKalman),Initialization,optimset('MaxIter',5000,'TolX',1e-2,'TolFun',1e-8,'MaxFunEvals',5000));
    Names = ParametersKalman.Names.Estimated;
    for i = 1:length(Names)
        ParametersKalman.(Names{i}).TransfValue = (x(i));
    end
    ParametersKalman = UpdateParsTransfToNoTransf(ParametersKalman);
    TellParsValues(ParametersKalman)
    ParametersKalman.Correction = 1;

    ParametersKalman.KalCov = Cov;


    Temp = struct();
    Temp.ParametersKalman = ParametersKalman;
    save([SavePath '/' Parameters.TempName],'Temp')

    %  die


    SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';
    % SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan/Temp'

    load([SavePath '/' Parameters.TempName])

    Parameters = Temp.ParametersKalman;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('SMC Maxim. alg. only sigma')
% Parameters.NbParticules = 1000;
% Parameters.StableCUseConstraint = 0;
% Parameters.NoPaths = 0;
% Parameters.PathsToKeep = [7 8 9];
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(AlreadySomething)
    disp('Quick PMCMC')
    Names = Parameters.Names.Estimated;
    for i = 1:length(Names);
        Parameters.(Names{i}).Estimated = 0;
    end
    if strcmp(Parameters.DiffusionType,'Bertallanfy')
        Parameters.BRmm1.Estimated = 1;    
        Parameters.InitialIPropF.Estimated = 1;
        Parameters.InitialIPropM.Estimated = 1;
        Parameters.Alpham1.Estimated = 1;
        Parameters.MuFm1.Estimated = 1;
        Parameters.MuMm1.Estimated = 1;
        Parameters.BetaMFPerAct.Estimated = 1;
        Parameters.BetaFMPerAct.Estimated = 1;
        Parameters.NumberActsPerClient.Estimated = 1;
        Parameters.eHIV.Estimated = 1;
        Parameters.CF1.Estimated = 1;
        Parameters.CF2.Estimated = 1;
        Parameters.CM.Estimated = 1;
        Parameters.BRbase.Estimated = 1;
        Parameters.BRmu.Estimated = 1;
        Parameters.BRtinfl.Estimated = 1;
        Parameters.BRsigma.Estimated = 1;
    elseif strcmp(Parameters.DiffusionType,'Sigmoid')
        Parameters.InitialIPropF.Estimated = 1;
        Parameters.InitialIPropM.Estimated = 1;
        Parameters.Alpham1.Estimated = 1;
        Parameters.MuFm1.Estimated = 1;
        Parameters.MuMm1.Estimated = 1;
        Parameters.BetaMFPerAct.Estimated = 1;
        Parameters.BetaFMPerAct.Estimated = 1;
        Parameters.NumberActsPerClient.Estimated = 1;
        Parameters.eHIV.Estimated = 1;
        Parameters.CF1.Estimated = 1;
        Parameters.CF2.Estimated = 1;
        Parameters.CM.Estimated = 1;
        Parameters.Sigmbase.Estimated = 1;
        Parameters.Sigmmu.Estimated = 1;
        Parameters.Sigmtinfl.Estimated = 1;
        Parameters.Sigmrate.Estimated = 1;
        Parameters.Sigmsigma.Estimated = 1;
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
        Parameters.Alpham1.Estimated = 1;
        Parameters.MuFm1.Estimated = 1;
        Parameters.MuMm1.Estimated = 1;
        Parameters.BetaMFPerAct.Estimated = 1;
        Parameters.BetaFMPerAct.Estimated = 1;
        Parameters.NumberActsPerClient.Estimated = 1;
        Parameters.eHIV.Estimated = 1;
        Parameters.CF1.Estimated = 1;
        Parameters.CF2.Estimated = 1;
        Parameters.CM.Estimated = 1;
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
    Cov =  Parameters.KalCov;

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
%     Res.Parameters = Parameters;

    for i = 1:3
        disp(['MCMC ' num2str(i)])
        if Res.AccRate>0.05
            Cov =  cov(Res.TransfThetas');
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
    %     Parameters.Kee??pAll = 1;
        Parameters.AdaptC = 0.999;
        % [ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
        Res = RunEstimationMethod(Data, HIVModel,Parameters,TempPar,5000);
%         Res.Parameters = Parameters;
    end

    load([SavePath '/' Parameters.TempName])
    Temp.ParametersPMCMC = Parameters;
    Temp.ResPMCMCNoPaths = Res;
    save([SavePath '/' Parameters.TempName],'Temp')
end

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

Cov =  cov(Res.TransfThetas');


% load([SavePath Parameters.NameToSave])
% Cov = cov(Res.TransfThetas');


% OldCov = Cov;
% Cov =  2.38^2/length(Names)*DetCov;
% Cov(Parameters.SigmaRW.Index,Parameters.SigmaRW.Index) = OldCov(13,13);
Parameters.G = Cov^-1;
Parameters.NoPaths = 0;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
% try
%     if Parameters.SwitchToIBM == 1
%         Parameters.DiffusionType = 'IBM';
%     end
% catch    
%     Parameters.DiffusionType = 'Add';
% end

Parameters.SaveForMarcMCMC = 0;

Parameters.aim = 0.23;
Parameters.Epsil = 1;
Parameters.ModelType = 'SMC';
Parameters.Correction = 1;
TempPar = Res.TempPar;
% Parameters.KeepAll = 1;
Parameters.AdaptC = 0.999;
% [ParametersPMCMC, TempPar] = CalibrateMethod( Data, HIVModel, ParametersPMCMC, TempPar);
Parameters.SaveSpace = 1;
Res = RunEstimationMethod(Data, HIVModel,Parameters,TempPar,NbItsPMCMC);
Res.Parameters = Parameters;

try
    save([SavePath Parameters.NameToSave '_' Parameters.DiffusionType '_' num2str(Parameters.NbRounds) '_' num2str(Parameters.TakeClients) '.mat'],'Res')
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