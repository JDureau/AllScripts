% Main Sigm process HIV



%% load paths

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
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\HIV')

%%


SavePath = 'S:\Results';
temp = load([SavePath '\ParametersMysore.mat']);
ParametersMysore = temp.Parameters;


SavePath = 'S:\Results';
temp = load([SavePath '\ParametersBelgaum.mat']);
ParametersBelgaum = temp.Parameters;



ParametersBelgaum.NbVariables = 10;
ParametersBelgaum.SigmaObs = 0.1;
ParametersBelgaum.Problem = 'ImperialHIV';
ParametersBelgaum.DiffusionType = 'Affine';
ParametersBelgaum.ObservationLength = 25*12;
ParametersBelgaum.ComputationTStep = 0.5;

ParametersMysore.NbVariables = 10;
ParametersMysore.SigmaObs = 0.1;
ParametersMysore.Problem = 'ImperialHIV';
ParametersMysore.DiffusionType = 'Affine';
ParametersMysore.ObservationLength = 25*12;
ParametersMysore.ComputationTStep = 0.5;
Parameters = ParametersMysore;

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
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';
HIVModel.SMC_projection = @HIV_SMC_projection;


% SavePath = 'S:\Results';
% Res = load([SavePath '\HIV_Mysore_PlayingPriors.mat']);
% Res = Res.Ress;
% Res = Res{17,3};


% Parameters = Res.Parameters;
Names = Parameters.Names.All;
for i = 1:length(Names)
    Parameters.(Names{i}).Estimated = 1;
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.TypeWork = 'Normal';

Parameters.ObsNoise = 0;
Ress = {};
Model = HIVModel;
% 
% trajs = [];
% ampls = [];
% for indTest = 1:NbTests
%     disp(indTest)
%     ResGen = GenerateDataForGivenSig(Data,Parameters,HIVModel);
%     trajs(indTest,:) = ResGen.Data.BuiltTraj(:,9);    
%     ampls(indTest) = ResGen.Data.BuiltTraj(end,9)-ResGen.Data.BuiltTraj(1,9);
% end
% subplot(2,1,1)
% plot(trajs')
% subplot(2,1,2)
% hist(ampls)

NbTests = 400;

ResGens = {};

for indTest = 15:NbTests
    disp(indTest)
    ResGen = GenerateDataForGivenSig(Data,ParametersMysore,HIVModel);
    ResGens{indTest} = ResGen;
% end
%     
% 
% ampls = [];
% infl = [];
% for indTest = 1:100%NbTests
%     ampls(indTest) = ResGens{indTest}.ampl;
%     infl(indTest) = ResGens{indTest}.infl;
% %     plot(ResGens{indTest}.Data.BuiltTraj(:,9))
% %     pause()
% end
% 
% plot(infl,ampls,'.')

    
    Data = ResGen.Data;
%     Parameters = ResGen.Parameters;
    % first, model parameters are fixed to MLE, we just estimate the diffusion
    % parameters
    
    Parameters = ResGen.Parameters;
%     
%     temp7 = zeros(1,10);
%     temp8 = zeros(1,10);
%     temp7(1,7) = 1;
%     temp8(1,8) = 1;
%     temps{7} = temp7;
%     temps{8} = temp8;
%     HIVModel.ObservationJacobian = {};
%     ObsVars = Data.ObservedVariables;
%     for i = 1:length(ObsVars)-1
%         HIVModel.ObservationJacobian{i+1} = temps{ObsVars(i+1)};
%     end
%     HIVModel.ObservationMeasurementNoise = {};
%     for i = 1:length(ObsVars)-1
%         HIVModel.ObservationMeasurementNoise{i+1} = (Data.Observations(ObsVars(i+1),i+1)*(100-Data.Observations(ObsVars(i+1),i+1))/400);
%     end
% 
%     % optimize all parameters.
%     Names = Parameters.Names.All;
%     for i = 1:length(Names);
%         Parameters.(Names{i}).Estimated = 0;
%     end
%     Parameters.InitialFt.Estimated = 1;
%     Parameters.SecondFt.Estimated = 1;
%     Parameters.ThirdFt.Estimated = 1;
%     Parameters.FourthFt.Estimated = 1;
%     Parameters.FifthFt.Estimated = 1;
%     Parameters.SixthFt.Estimated = 1;
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
%     Parameters = DefineEstimatedParametersIndexes(Parameters);
%     Parameters = DefineTransfFunctions(Parameters);
%     Parameters = UpdateParsNoTransfToTransf(Parameters);
%     Parameters = DefinePriors(Parameters);
%     Parameters.Problem = 'ImperialHIV';
%     Parameters.DiffusionType = 'AffineAdd';
%     Names = Parameters.Names.Estimated;
%     
%     test = 0;
%     for k = 1:length(Names)
%         Parameters.(Names{k}).Sample = 1;
%     end
%     while not(test)
%         Parameters = SampleParameters(Parameters);
%         try
%             Temp = EstimationEKFGen(Data, HIVModel, Parameters);
%             if Temp.LogLik>-1000
%                 test = 1;
%             end
%         end
%         NbIts = NbIts+1;
%         if NbIts == 10000
%             test = 1;
%         end
%     end
%         
%     
%     Initialization = [];
%     for i = 1:length(Names)
%         Initialization(i) = Parameters.(Names{i}).TransfValue ;
%     end
%     [x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,Parameters),Initialization,optimset('MaxIter',5000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',2000));
%     Names = Parameters.Names.Estimated;
%     for i = 1:length(Names)
%         Parameters.(Names{i}).TransfValue = (x(i));
%     end
%     Parameters = UpdateParsTransfToNoTransf(Parameters);
%     TellParsValues(Parameters)
% 
% 
%     
    
    Parameters.NbParticules = 1000;
    Parameters.StableCUseConstraint = 0;
    Parameters.NoPaths = 0;
    Parameters.PathsToKeep = [7 8 9];
    Parameters.DiffusionType = 'Add';
    Parameters.MCMCType = 'jiji';
    Parameters.GMeth = 'cst given';
    Names = Parameters.Names.All;
    for i = 1:length(Names);
        Parameters.(Names{i}).Estimated = 0;
    end

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


    

    %% Filter

    disp('Step 3: Sigma Opt')

    Parameters.NbParticules = 5000;
    Parameters.DiffusionType = 'Add';
    Res = EstimationSMCfiltGen(Data, Model, Parameters);
    Res.LogLik

    NbIts = 100;
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

    Ress{indTest} = Res;
%     PlotResHIV(Res,Parameters)
%     temp7 = zeros(1,10);
%     temp8 = zeros(1,10);
%     temp7(1,7) = 1;
%     temp8(1,8) = 1;
%     HIVModel.ObservationJacobian = {};
%     HIVModel.ObservationJacobian{2} = temp7;
%     HIVModel.ObservationJacobian{3} = temp7;
%     HIVModel.ObservationJacobian{4} = temp8;
%     HIVModel.ObservationJacobian{5} = temp7;
%     HIVModel.ObservationMeasurementNoise = {};
%     HIVModel.ObservationMeasurementNoise{2} = (DataGen.Observations(7,2)*(100-DataGen.Observations(7,2))/(400));
%     HIVModel.ObservationMeasurementNoise{3} = (DataGen.Observations(7,3)*(100-DataGen.Observations(7,3))/(400));
%     HIVModel.ObservationMeasurementNoise{4} = (DataGen.Observations(8,4)*(100-DataGen.Observations(8,4))/(400));
%     HIVModel.ObservationMeasurementNoise{5} = (DataGen.Observations(7,5)*(100-DataGen.Observations(7,5))/(400));
% 
%     Ress{i} = QuickHIVfilteredPost(DataGen,Parameters,HIVModel);
    
    SavePath = 'S:\Results\';
    save([SavePath 'VIH_PlayingWithSigm_NoParEst3.mat'],'Ress')

end

Resss = {};
inds = {'','2','3'};
for i = 1:3
    SavePath = 'S:\Results\';
    load([SavePath 'VIH_PlayingWithSigm' inds{i} '.mat'])
    Resss{i} = Ress;
end


MeanErrors = [];
Errors = [];
ampl = [];
slope = [];
inflpt = [];
for i = 1:3
    for j = 1:length(Resss{i})
        Res = Resss{i}{j};
        errors = zeros(1,584);
        for k = 1:584
            errors(k) = mean((Res.Paths(:,3,k)-Res.Data.BuiltTraj(k,9)).^2);
        end
        Res.Errors = errors;
        Res.MeanError = mean(errors);
        MeanErrors(i,j) = mean(errors);
        Errors(i,j,:) = errors;
        ampl(i,j) = Res.Data.BuiltTraj(584,9)-Res.Data.BuiltTraj(1,9);
        tmp = diff(Res.Data.BuiltTraj(:,9));
        [b,ind] = max(tmp);
        inflpt(i,j) = ind;
        slope(i,j) = (Res.Data.BuiltTraj(ind+1,9)-Res.Data.BuiltTraj(ind-1,9))/(2*Parameters.ComputationTStep);
    end
end

clf
[fi,xi] = ksdensity(MeanErrors(1,:));
plot(xi,fi)
hold on
[fi,xi] = ksdensity(MeanErrors(2,:));
plot(xi,fi,'k')
[fi,xi] = ksdensity(MeanErrors(3,:));
plot(xi,fi,'g')
hold off

[b,ind] = min(MeanErrors(3,:));
PlotResHIV(Resss{3}{ind},Resss{3}{ind}.Parameters,Resss{3}{ind}.Data.BuiltTraj(:,9))
[b,ind] = max(MeanErrors(3,:));
PlotResHIV(Resss{3}{ind},Resss{3}{ind}.Parameters,Resss{3}{ind}.Data.BuiltTraj(:,9))

hist([MeanErrors(3,:) MeanErrors(2,:)])
xlabel('Mean Squared error')
title('Distribution of Mean Squared Errors')

xs = [];
ys = [];
for j = 1:length(Resss{3})
    Res = Resss{i}{j};
    for k = 1:584
        xs(end+1) = k;
        ys(end+1) = Errors(3,j,k);
    end
end
scattercloud(xs*Parameters.ComputationTStep,ys)
xlabel('time')
ylabel('squared errors at time t')
title('Distribution of squared errors along time')

% role of ampl
plot(ampl(3,:),MeanErrors(3,:),'.')
hold on
plot(ampl(2,:),MeanErrors(2,:),'.')
plot(ampl(1,:),MeanErrors(1,:),'.')
hold off
xlabel('Amplitude of Ft trajectory')
ylabel('Mean Squarred Error')
title('Role of amplitude')

% role of t
plot(inflpt(3,:),MeanErrors(3,:),'.')
hold on
plot(inflpt(2,:),MeanErrors(2,:),'.')
plot(inflpt(1,:),MeanErrors(1,:),'.')
hold off
xlabel('Inflection point trajectory')
ylabel('Mean Squarred Error')
title('Role of inflection point')


PlotScatteredQuantity([([inflpt(3,:) inflpt(2,:)]*Parameters.ComputationTStep)' [ampl(3,:) ampl(2,:)]'],-[MeanErrors(3,:) MeanErrors(2,:)],6)
xlabel('Inflection Point')
ylabel('Amplitude')
title('Combined role of amplitude and time of inflection', 'FontWeight','bold')
colorbar
colorbar('YTickLabel', {'Worse','','','','','Average','','','','Best'})



ind = 59;
traj = Resss{3}{ind}.Data.BuiltTraj(:,9);
RealDists = [];
EstDists  = [];
for i = 2:3
    for j = 1:100
        RealDists(end+1) = mean((traj-Resss{i}{j}.Data.BuiltTraj(:,9)).^2);
        errors = zeros(1,584);
        Res = Resss{i}{j};
        for k = 1:584
            errors(k) = mean((Res.Paths(:,3,k)-traj(k)).^2);
        end
        EstDists(end+1) = mean(errors);
    end
end

clf
plot(RealDists,EstDists,'.')
xlabel('Real Mean Squarred Distance to Avahan-like Trajectories')
ylabel('Estimated Mean Squarred Distance to Avahan-like Trajectories')
title('How well does this procedure discriminate trajectories like the ones estimated for Avahan?', 'FontWeight','bold')


RealInflPts = [];
RealEnd = [];
RealBaseLine = [];
EstInflPts = [];
EstEnd = [];
EstBaseLine = [];
InteresInds = [];
for i = 2:3
    for j = 1:100
        Res = Resss{i}{j};
        MeanTraj = mean(Res.Paths(:,3,:));
        RealInflPts(end+1) = inflpt(i,j);
        tmp = diff(MeanTraj);
        [b,ind] = max(tmp);
        EstInflPts(end+1) = ind;
        RealBaseLine(end+1) = Res.Data.BuiltTraj(1,9);
        EstBaseLine(end+1) = MeanTraj(1);
        RealEnd(end+1) = Res.Data.BuiltTraj(end,9);
        EstEnd(end+1) = MeanTraj(end);
        if and(ampl(i,j)>0.5,inflpt>420)
            InteresInds = length(EstEnd);
        end
    end
end


plot(RealDists,((-RealInflPts+EstInflPts)*Parameters.ComputationTStep),'.')
xlabel('Mean Squarred Distance to Avahan Estimated Trajectories')
ylabel('Error in inflection point estimation (in months) (Estimated-Real)')
title('How well is the inflection point estimated for Avahan-like trajectories?', 'FontWeight','bold')

plot(RealDists,abs((RealBaseLine-EstBaseLine)),'.')
xlabel('Mean Squarred Distance to Avahan Estimated Trajectories')
ylabel('Error in condom use baseline estimation ')
title('How well is the condom use baseline estimated for Avahan-like trajectories?', 'FontWeight','bold')

plot(RealDists,abs((RealEnd-EstEnd)),'.')
xlabel('Mean Squarred Distance to Avahan Estimated Trajectories')
ylabel('Error in condom use 2010 estimation ')
title('How well is the condom use in 2010 estimated for Avahan-like trajectories?', 'FontWeight','bold')





SavePath = 'S:\Results\';
load([SavePath 'VIH_PlayingWithSigm2.mat'])

MeanVals = [];
AmplCI = [];
Ress  = Resss{3};
for i = 1:3
    for j = 1:length(Resss{i})
        if ampl(i,j)<0.2
            PlotResHIV(Resss{i}{j},Resss{i}{j}.Parameters,Resss{i}{j}.Data.BuiltTraj(:,9))
            title([i j])
            MeanVals(end+1) = mean(Resss{i}{j}.Data.BuiltTraj(:,9));
            AmplCI(end+1) = mean(quantile(Resss{i}{j}.Paths(:,3,:),0.75)-quantile(Resss{i}{j}.Paths(:,3,:),0.25));
            pause()
        end
    end        
end

clf
plot(abs(MeanVals-0.5),AmplCI,'.')
xlabel('Distance between condom use contant value and 0.5')
ylabel('Width of the mean posterior 50% confidence interval')
title('Estimation of constant condom use improve when this value is far from 0.5', 'FontWeight','bold')


hist(max(squeeze(Resss{1}{7}.Paths(:,3,:))')-min(squeeze(Resss{1}{7}.Paths(:,3,:))'))

PlotResHIV(Resss{1}{7},Resss{1}{7}.Parameters,Resss{1}{7}.Data.BuiltTraj(:,9))

hist(max(squeeze(Resss{2}{17}.Paths(:,3,:))')-min(squeeze(Resss{2}{17}.Paths(:,3,:))'))

PlotResHIV(Resss{2}{17},Resss{2}{17}.Parameters,Resss{2}{17}.Data.BuiltTraj(:,9))

clf
inds = ceil(rand(1,2)*100);
plot(squeeze(Resss{2}{17}.Paths(inds,3,:))','.')
title('fig 12: 2 random sample paths')

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
    
