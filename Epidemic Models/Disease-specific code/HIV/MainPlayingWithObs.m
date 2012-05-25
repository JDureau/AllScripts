% Main Obs process HIV


Parameters = struct();

% Initializing all model parameters
Parameters.TotalFSW.Value = 2144;
Parameters.TotalFSW.Min = 804;
Parameters.TotalFSW.Max = 3752;
Parameters.TotalFSW.Estimates = 0;
Parameters.TotalFSW.TransfType = 'Log';
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialFt.Value = 0.08;
Parameters.InitialFt.Min = -10^14;
Parameters.InitialFt.Max = 10^14;
Parameters.InitialFt.MaxLim = 0.9;
Parameters.InitialFt.MinLim = 0.03;
Parameters.InitialFt.Estimated = 1;
Parameters.InitialFt.TransfType = 'Logit';
Parameters.InitialFt.Init = 1;
Parameters.SecondFt.Value = 0.08;
Parameters.SecondFt.Min = -10^14;
Parameters.SecondFt.Max = 10^14;
Parameters.SecondFt.MaxLim = 0.9;
Parameters.SecondFt.MinLim = 0.03;
Parameters.SecondFt.Estimated = 1;
Parameters.SecondFt.TransfType = 'Logit';
Parameters.ThirdFt.Value = 0.08;
Parameters.ThirdFt.Min = -10^14;
Parameters.ThirdFt.Max = 10^14;
Parameters.ThirdFt.MaxLim = 0.96;
Parameters.ThirdFt.MinLim = 0.03;
Parameters.ThirdFt.Estimated = 1;
Parameters.ThirdFt.TransfType = 'Logit';
Parameters.FourthFt.Value = 0.7;
Parameters.FourthFt.Min = -10^14;
Parameters.FourthFt.Max = 10^14;
Parameters.FourthFt.MaxLim = 0.96;
Parameters.FourthFt.MinLim = 0.03;
Parameters.FourthFt.Estimated = 1;
Parameters.FourthFt.TransfType = 'Logit';
Parameters.FifthFt.Value = 0.7;
Parameters.FifthFt.Min = -10^14;
Parameters.FifthFt.Max = 10^14;
Parameters.FifthFt.MaxLim = 0.96;
Parameters.FifthFt.MinLim = 0.03;
Parameters.FifthFt.Estimated = 1;
Parameters.FifthFt.TransfType = 'Logit';
Parameters.SixthFt.Value = 0.7;
Parameters.SixthFt.Min = -10^14;
Parameters.SixthFt.Max = 10^14;
Parameters.SixthFt.MaxLim = 0.96;
Parameters.SixthFt.MinLim = 0.03;
Parameters.SixthFt.Estimated = 1;
Parameters.SixthFt.TransfType = 'Logit';
Parameters.InitialIPropF.Value = 0.002;
Parameters.InitialIPropF.Min = 0;
Parameters.InitialIPropF.Max = 0.04;
Parameters.InitialIPropF.MaxLim = 0.10;
Parameters.InitialIPropF.MinLim = 0;
Parameters.InitialIPropF.Estimated = 1;
Parameters.InitialIPropF.TransfType = 'Logit';
Parameters.InitialIPropF.Init = 1;
Parameters.InitialIPropM.Value = 0.002;
Parameters.InitialIPropM.Min = 0;
Parameters.InitialIPropM.Max = 0.02;
Parameters.InitialIPropM.MaxLim = 0.1;
Parameters.InitialIPropM.MinLim = 0;
Parameters.InitialIPropM.Estimated = 1;
Parameters.InitialIPropM.TransfType = 'Logit';
Parameters.InitialIPropM.Init = 1;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotMFactor.Value = 12.3;
Parameters.TotMFactor.Min = -10^14;
Parameters.TotMFactor.Max = 10^14;
Parameters.TotMFactor.MinLim = 7;
Parameters.TotMFactor.MaxLim = 19;
Parameters.TotMFactor.Estimated = 1;
Parameters.TotMFactor.TransfType = 'Logit';
Parameters.TotMFactor.Init = 1;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.Alpham1.Value = 93;
Parameters.Alpham1.Min = -10^14;
Parameters.Alpham1.Max = 10^14;
Parameters.Alpham1.MinLim = 87; % was 92
Parameters.Alpham1.MaxLim = 138.5;
Parameters.Alpham1.Estimated = 1;
Parameters.Alpham1.TransfType = 'Logit';
Parameters.Alpham1.PlotInv = 1;
Parameters.MuFm1.Value = 54;
Parameters.MuFm1.Min = -10^14;
Parameters.MuFm1.Max =  10^14;
Parameters.MuFm1.MinLim = 45; % was 53 switched to 45 (street based might lower the average)
Parameters.MuFm1.MaxLim = 66;
Parameters.MuFm1.Estimated = 1;
Parameters.MuFm1.TransfType = 'Logit';
Parameters.MuFm1.PlotInv = 1;
Parameters.MuMm1.Value = 90;
Parameters.MuMm1.Min = -10^14;
Parameters.MuMm1.Max = 10^14;
Parameters.MuMm1.MinLim = 84;
Parameters.MuMm1.MaxLim = 240;
Parameters.MuMm1.Estimated = 1;
Parameters.MuMm1.PlotInv = 1;
Parameters.MuMm1.TransfType = 'Logit';
Parameters.BetaMFPerAct.Value = 0.0032;
Parameters.BetaMFPerAct.Min = -10^14;
Parameters.BetaMFPerAct.Max =  10^14;
Parameters.BetaMFPerAct.MinLim = 0.0006*2;
Parameters.BetaMFPerAct.MaxLim = 0.0011*5;
Parameters.BetaMFPerAct.Estimated = 1;
Parameters.BetaMFPerAct.TransfType = 'Logit';
Parameters.BetaFMPerAct.Value = 0.0048;
Parameters.BetaFMPerAct.Min = -10^14;
Parameters.BetaFMPerAct.Max =  10^14;
Parameters.BetaFMPerAct.MinLim = 0.0001*2;
Parameters.BetaFMPerAct.MaxLim = 0.0014*5;
Parameters.BetaFMPerAct.Estimated = 1;
Parameters.BetaFMPerAct.TransfType = 'Logit';
Parameters.NumberActsPerClient.Value = 1.7;
Parameters.NumberActsPerClient.Min = -10^14;
Parameters.NumberActsPerClient.Max =  10^14;
Parameters.NumberActsPerClient.MinLim = 1;
Parameters.NumberActsPerClient.MaxLim = 3;
Parameters.NumberActsPerClient.Estimated = 1;
Parameters.NumberActsPerClient.TransfType = 'Logit';
Parameters.BetaFM.Value = 1-(1-Parameters.BetaFMPerAct.Value)^Parameters.NumberActsPerClient.Value;
Parameters.BetaMF.Value = 1-(1-Parameters.BetaMFPerAct.Value)^Parameters.NumberActsPerClient.Value;
Parameters.eHIV.Value = 0.86;
Parameters.eHIV.Min = -10^14;
Parameters.eHIV.Max =  10^14;
Parameters.eHIV.MaxLim = 0.9;
Parameters.eHIV.MinLim = 0.8;
Parameters.eHIV.Estimated = 1;
Parameters.eHIV.TransfType = 'Logit';
Parameters.CF1.Value = 16.3;
Parameters.CF1.Min = -10^14;
Parameters.CF1.Max =  10^14;
Parameters.CF1.MinLim = 15.22;
Parameters.CF1.MaxLim = 17.3;
Parameters.CF1.Estimated = 1;
Parameters.CF1.TransfType = 'Logit';
Parameters.CF2.Value = 53.3;
Parameters.CF2.Min = -10^14;
Parameters.CF2.Max = 10^14;
Parameters.CF2.MinLim = 47.37;
Parameters.CF2.MaxLim =  56.6;
Parameters.CF2.Estimated = 1;
Parameters.CF2.TransfType = 'Logit';
Parameters.SigmaRW.Value = 0.1;
Parameters.SigmaRW.Min = -10^14;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.MinLim = 0;
Parameters.SigmaRW.MaxLim = 50;
Parameters.SigmaRW.Estimated = 0;
Parameters.SigmaRW.TransfType = 'Log';
Parameters.InitialDeriv = 0;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.StableCUseConstraint = 0;
TellParsValues(Parameters)

Parameters.NbVariables = 10;
Parameters.SigmaObs = 0.1;
Parameters.Problem = 'ImperialHIV';
Parameters.DiffusionType = 'Affine';
Parameters.ObservationLength = 25*12;
Parameters.ComputationTStep = 0.5;


% t0 = jan 85.
Data.Observations = zeros(10,5);
Data.Observations(7,2) = 26.11;
Data.Observations(7,3) = 24.24;
Data.Observations(8,4) = 5.4;
Data.Observations(7,5) = 11.10;
Data.Instants = round([0 236 264 286 292]/(Parameters.ComputationTStep));
Data.ObservedVariables = [ 0 7 7 8 7];
Data.NbComputingSteps = [0 diff(Data.Instants)];

HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIVInitialize;
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*Parameters.SigmaObs).*(Res.WentOutOrNot)';
HIVModel.SMC_projection = @HIV_SMC_projection;




SavePath = 'S:\Results';
Res = load([SavePath '\HIV_Mysore_PlayingPriors.mat']);
Res = Res.Ress;
Res = Res{17,3};

Parameters = Res.Parameters;
Trajs = {};

% cst low
Parameters.ComputationTStep = 0.5;
Parameters.NbTSteps = 1200;
Fts = 0.1*ones(1,Parameters.NbTSteps);
HIV_CreateData(Fts,Parameters,HIVModel,Data)

ind = 1;
Trajs{ind} = struct();
Trajs{ind}.Fts = Fts;
Trajs{ind}.Parameters = Parameters;

%cst high: extinction
Parameters.ComputationTStep = 0.5;
Parameters.NbTSteps = 3000;
Fts = 0.5*ones(1,Parameters.NbTSteps);
HIV_CreateData(Fts,Parameters,HIVModel,Data)

ind = 2;
Trajs{ind} = struct();
Trajs{ind}.Fts = Fts;
Trajs{ind}.Parameters = Parameters;

% sigmoid: realistic
Parameters.ComputationTStep = 0.5;
Parameters.NbTSteps = 100;
tmp = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ComputationTStep*Parameters.NbTSteps;
Fts = sigmoid((tmp-40)/2)*0.6 +0.01;
plot( Fts)
HIV_CreateData(Fts,Parameters,HIVModel,Data)

ind = 3;
Trajs{ind} = struct();
Trajs{ind}.Fts = Fts;
Trajs{ind}.Parameters = Parameters;


% sigmoid: higher and slower
Parameters.ComputationTStep = 0.5;
Parameters.NbTSteps = 120;
tmp = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ComputationTStep*Parameters.NbTSteps;
Fts = sigmoid((tmp-40)/5)*0.6 +0.3;
plot( Fts)
HIV_CreateData(Fts,Parameters,HIVModel,Data)

ind = 4;
Trajs{ind} = struct();
Trajs{ind}.Fts = Fts;
Trajs{ind}.Parameters = Parameters;

% random walk 
Parameters.ComputationTStep = 0.5;
Parameters.NbTSteps = 1200;
tmp = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ComputationTStep*Parameters.NbTSteps;
Fts = min(1,max(0,0.1+cumsum(0.01*randn(1,Parameters.NbTSteps)*sqrt(Parameters.ComputationTStep))));
plot( Fts)
HIV_CreateData(Fts,Parameters,HIVModel,Data)

ind = 5;
Trajs{ind} = struct();
Trajs{ind}.Fts = Fts;
Trajs{ind}.Parameters = Parameters;


ObsPatterns{1}.NbObs = 4;
ObsPatterns{1}.NbPatterns = 15;
ObsPatterns{2}.NbObs = 5;
ObsPatterns{2}.NbPatterns = 15;
for i = 1:2
    for j = 1:ObsPatterns{i}.NbPatterns 
        rdinds = sort(rand(1,ObsPatterns{i}.NbObs));
        rdclientorfsw = round(rand(1,ObsPatterns{i}.NbObs));
        for indtraj = 1:length(Trajs)
            Parameters = Trajs{indtraj}.Parameters;
            Fts = Trajs{indtraj}.Fts;
            ObsInds = ceil(rdinds*Parameters.NbTSteps);
            Data = HIV_CreateData(Fts,Parameters,HIVModel,Data);
            Data.Observations = zeros(10,ObsPatterns{i}.NbObs+2);
            Data.Instants = zeros(1,ObsPatterns{i}.NbObs+2);
            for indobs = 1:ObsPatterns{i}.NbObs
                Data.Observations(7+rdclientorfsw(indobs),indobs+1) = Data.BuiltTraj(ObsInds(indobs),7+rdclientorfsw(indobs));
                Data.Instants(indobs+1) = ObsInds(indobs);
                Data.ObservedVariables(indobs+1) = (7+rdclientorfsw(indobs));
            end
            Data.Observations(8,ObsPatterns{i}.NbObs+2)= 0;
            Data.Instants(ObsPatterns{i}.NbObs+2) = size(Data.BuiltTraj,1);
            Data.ObservedVariables(ObsPatterns{i}.NbObs+2) = 8; 
            Data.NbComputingSteps = [0 diff(Data.Instants)];
            Trajs{indtraj}.Datas{i,j} = Data;
            Model = HIVModel;
            tmp1 = zeros(1,10);
            tmp2 = zeros(1,10);
            tmp1(1,7) = 1;
            tmp2(1,8) = 1;
            temp = {};
            temp{1} = tmp1;
            temp{2} = tmp2;
            
            Model.ObservationJacobian = {};
            for indobs = 1:ObsPatterns{i}.NbObs
                Model.ObservationJacobian{1+indobs} = temp{1+rdclientorfsw(indobs)};
                Model.ObservationMeasurementNoise{1+indobs} = (Parameters.SigmaObs*Data.Observations(7+rdclientorfsw(indobs),1+indobs))^2;
            end
            Trajs{indtraj}.Models{i,j} = Model;
        end
    end
end
            


% go            
for indObsPattern = 1:2
    for indtraj = 1:length(Trajs)
        for j = 1:ObsPatterns{indObsPattern}.NbPatterns      
            Data = Trajs{indtraj}.Datas{indObsPattern,j};
            Model = Trajs{indtraj}.Models{indObsPattern,j};
            Parameters = Trajs{indtraj}.Parameters;
            Names = Parameters.Names.All;
            for i = 1:length(Names);
                Parameters.(Names{i}).Estimated = 0;
            end
            Parameters.InitialFt.Value = Trajs{indtraj}.Fts(1);
            Parameters.SigmaRW.Value = 0.1;
            Parameters.SigmaRW.Min = 0;
            Parameters.SigmaRW.Max = 10^14;
            Parameters.SigmaRW.Estimated = 1;
            Parameters.SigmaRW.TransfType = 'Log';
            Parameters = DefineEstimatedParametersIndexes(Parameters);
            Parameters = DefineTransfFunctions(Parameters);
            Parameters = UpdateParsNoTransfToTransf(Parameters);
            Parameters = DefinePriors(Parameters);
            % optimize RW parameters with SMC
            Parameters.NoPaths = 0;
            Names = Parameters.Names.Estimated;
            Initialization = [];
            for i = 1:length(Names)
                Initialization(i) = Parameters.(Names{i}).TransfValue ;
            end
            [x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizewithPrior(x,Data,Model,Parameters),Initialization,optimset('MaxIter',50,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));
            Names = Parameters.Names.Estimated;
            for i = 1:length(Names)
                Parameters.(Names{i}).TransfValue = (x(i));
            end
            Parameters = UpdateParsTransfToNoTransf(Parameters);
            TellParsValues(Parameters)

            Parameters.TypeWork = 'PlayObs';
            Parameters.NbParticules = 5000;
            Parameters.DiffusionType = 'Add';
            Res = EstimationSMCfiltGen(Data, Model, Parameters);

            NbIts = 40;
            Paths = zeros(NbIts,length(Parameters.PathsToKeep),sum(Data.NbComputingSteps));
            LogLiks = [];
            for i = 1:NbIts
                disp(i)
                Temp = EstimationSMCsmoothGen(Data, Model, Parameters);
                ind = ceil(rand(1,1)*Parameters.NbParticules);
                Paths(i,:,:) = Temp.CompletePaths(ind,:,:);
                LogLiks(i) = Temp.LogLik;
            end
            Res.Paths = Paths;
            PlotResHIV(Res,Parameters)

            Trajs{indtraj}.Ress{indObsPattern,j} = Res;      
        end
    end
end


SavePath = 'S:\Results';
% SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
% save([SavePath '\HIV_Mysore_ForProj_02I0_Calib.mat'],'ResRW1_RW_Calib_Rds_1_2_3')
save([SavePath '\HIV_PlayingObsPatterns.mat'],'Trajs')

load([SavePath '\HIV_PlayingObsPatterns.mat'])

clf
indtraj = 5;
indObsPattern = 2;
NbPats = 4;%ObsPatterns{indObsPattern}.NbPatterns ; 
tmp = ceil(NbPats/3);
Data = Trajs{indtraj}.Datas{indObsPattern,1};
TStep = Trajs{indtraj}.Parameters.ComputationTStep;
NbTSteps = Trajs{indtraj}.Parameters.NbTSteps;
xis = TStep:Parameters.ComputationTStep:NbTSteps*TStep;
subplot(tmp*6,1,1:tmp)
plot(xis,Data.BuiltTraj(:,7),'m')
title('FSWs')
subplot(tmp*6,1,tmp+1:2*tmp)
plot(xis,Data.BuiltTraj(:,8),'k')
title('Clients')
subplot(tmp*6,1,2*tmp+1:3*tmp)
plot(xis,Data.BuiltTraj(:,9))
ylim([0 1])
title('Condom Use')
scores = [];
for i = 1:ObsPatterns{indObsPattern}.NbPatterns
    Res = Trajs{indtraj}.Ress{indObsPattern,i};
    x = 0;
    for j = 1:size(Res.Paths,1)
        x = x + (squeeze(Res.Paths(j,3,:)) - Data.BuiltTraj(1:size(Res.Paths,3),9)).^2;
    end        
    scores(i) = mean(x);
end
[bof,inds] = sort(scores);
inds = inds([1 2 ObsPatterns{indObsPattern}.NbPatterns-1 ObsPatterns{indObsPattern}.NbPatterns])
for i = 1:4
    ind = inds(i);
    Data = Trajs{indtraj}.Datas{indObsPattern,ind};
    subplot(tmp*6,1,3*tmp+i)
    tmpinds = find(Data.ObservedVariables == 7);
    if not(isempty(tmpinds))
        plot(Data.Instants(tmpinds),0,'o','MarkerEdgeColor','m','MarkerFaceColor','m')
    end
    hold on
    tmpinds = find(Data.ObservedVariables == 8);
    if indObsPattern == 2
        tmpinds = tmpinds(1:end-1);
    else
        tmpinds = tmpinds(1:end-2);
    end
    if not(isempty(tmpinds))
        plot(Data.Instants(tmpinds),0,'o','MarkerEdgeColor','k','MarkerFaceColor','k')
    end
    hold off
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',0)
    set(gca,'yticklabel',scores(ind))
%     clf
%     plot(Data.BuiltTraj(1:size(Res.Paths,3),9),'g')
%     hold on
%     Res = Trajs{indtraj}.Ress{indObsPattern,ind};
%     plot(mean(squeeze(Res.Paths(:,3,:))))
%     hold off
%     pause()
end    


inds = [1 3 4 5];
for i = 1:4
    indtraj = inds(i);
    indObsPattern = 2;
    NbPats = 4;%ObsPatterns{indObsPattern}.NbPatterns ; 
    tmp = ceil(NbPats/3);
    Data = Trajs{indtraj}.Datas{indObsPattern,1};
    TStep = Trajs{indtraj}.Parameters.ComputationTStep;
    NbTSteps = Trajs{indtraj}.Parameters.NbTSteps;
    xis = TStep:Parameters.ComputationTStep:NbTSteps*TStep;
    subplot(4,1,i)
    plot(xis,Data.BuiltTraj(:,9))
    ylim([0 1])
    title('Condom Use')
    xlabel('time')
end
    
