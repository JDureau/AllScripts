%% Tilting the priors


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

%% Mysore

% Setting

% Load Mysore data and priors


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
Parameters.SigmaRW.Estimated = 1;
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
Parameters.TypeWork = 'Normal';


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

temp7 = zeros(1,10);
temp8 = zeros(1,10);
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





Names = Parameters.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters.DiffusionType = 'AffineAdd';

[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,Parameters),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
Initialization = x;
Names = Parameters.Names.Estimated;
ParametersKalman = Parameters;
for i = 1:length(Names)
    ParametersKalman.(Names{i}).TransfValue = (x(i));
    Parameters.(Names{i}).TransfValue = (x(i));
end
ParametersKalman = UpdateParsTransfToNoTransf(ParametersKalman);
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(ParametersKalman)
Names = ParametersKalman.Names.Estimated;
for i = 1:length(Names)
    ParametersKalman.(Names{i}).SamplStd = 0.001*ParametersKalman.(Names{i}).Value;
end
ResKal = KalmanNumericDerivativesWithPrior(Data,HIVModel,ParametersKalman);
Test = mean(eig(-ResKal.Hess)>0)==1;
Cov = (-ResKal.Hess)^-1;
Test
% 
% Names = Parameters.Names.Estimated;
% Cov = zeros(length(Names),length(Names));
% for i = 1:length(Names)
%     if and(strcmp(Names{i}(end-1:end),'Ft'),not(strcmp(Names{i}(1:2),'In')))
%         Parameters.(Names{i}).Estimated = 0;
%     else
%         Parameters.(Names{i}).Estimated = 1;
%         Parameters.(Names{i}).Value = ParametersKalman.(Names{i}).Value;
%     end
% end
% Parameters.SigmaRW.Estimated = 1;
% Parameters.SigmaRW.Value = 0.1;
% Parameters = DefineEstimatedParametersIndexes(Parameters);
% Parameters = DefineTransfFunctions(Parameters);
% Parameters = DefinePriors(Parameters);
% Parameters = UpdateParsNoTransfToTransf(Parameters);
% Parameters.StableCUseConstraint = 0;
% Names = Parameters.Names.Estimated;
% Cov = zeros(length(Names),length(Names));
% 
% inds = [];
% indsKal = [];
% for i = 1:length(Names)
%     if ParametersKalman.(Names{i}).Estimated
%         indsKal(i) =  ParametersKalman.(Names{i}).Index;
%         inds(i) = Parameters.(Names{i}).Index;
%     end
% end
% tmp = (-ResKal.Hess)^-1;
% Cov(inds,inds) = tmp(indsKal,indsKal);
% indsig = Parameters.SigmaRW.Index;
% Cov(indsig,indsig) = 0.1;
%     
% 
% Parameters.MIFNbIterations = 100;
% Parameters.MIFNbParticules = 10000;
% Parameters.MIFCoolingParameters = 0.000005^(1/Parameters.MIFNbIterations);
% Parameters.MIFb = 1;
% Parameters.MIFSigmas = [];
% Parameters.NoPaths = 0;
% Parameters.PathsToKeep = [1:9]';
% Parameters.MCMCType = 'gtgt';
% Names = Parameters.Names.Estimated;
% Parameters.MIFCov = 0.0005*(Cov);
% 
% ResMIF = MIFHIV(Data,Parameters,HIVModel);
% 
% 
% Names = Parameters.Names.Estimated;
% for i = 1:length(Names)
%     ind = Parameters.(Names{i}).Index;
%     Parameters.(Names{i}).TransfValue = ResMIF.ThetasRecord(ind,end);
% end
% Parameters= UpdateParsTransfToNoTransf(Parameters);

Names = Parameters.Names.All;
for i = 1:length(Names);
    Parameters.(Names{i}).Estimated = 0;
end
Parameters.SigmaRW.Value = 0.1;
Parameters.SigmaRW.Min = 0;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Log';
Parameters.NbParticules = 1000;
Parameters.StableCUseConstraint = 0;
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [7 8 9];
Parameters.DiffusionType = 'Add';
Parameters.MCMCType = 'jiji';
Parameters.GMeth = 'cst given';
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
Parameters.NbParticules = 1000;

[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizewithPrior(x,Data,Model,Parameters),Initialization,optimset('MaxIter',50,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)




Res = EstimationSMCfiltGen(Data, Model, Parameters);
Res.LogLik
PlotResHIV(Res)

NbIts = 20;
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

PlotResHIV(Res)



%% Playing with priors, one at a time.


% Namestmp = Parameters.Names.Estimated;
% Names = {};
% for i = 2:length(Namestmp)
%     Names{i-1} = Namestmp{i};
% end
Namestmp = Parameters.Names.Estimated;
Names = {};
for i = 7:length(Namestmp)-1  % be careful
    Names{i-6} = Namestmp{i};
end

ParametersTilted = Parameters;

Ress = {};
lambda = 1;
for i = 1:length(Names)
    ind = Parameters.(Names{i}).Index;
    
    if Parameters.(Names{i}).Max == 10^14
        ParametersTilted = Parameters;
        ampl = Parameters.(Names{i}).MaxLim - Parameters.(Names{i}).MinLim;
        mean = (Parameters.(Names{i}).MaxLim + Parameters.(Names{i}).MinLim)/2;


        % plus
        if strcmp(Names{i},'eHIV')
            ParametersTilted.(Names{i}).MaxLim = min(1,Parameters.(Names{i}).MaxLim + lambda*ampl);
        else
            ParametersTilted.(Names{i}).MaxLim = Parameters.(Names{i}).MaxLim + lambda*ampl;
        end
        ParametersTilted.(Names{i}).MinLim = max(0.1*mean,Parameters.(Names{i}).MinLim + lambda*ampl);
        ParametersTilted.(Names{i}).Estimated = 1;
        ParametersTilted.(Names{i}).Value = (ParametersTilted.(Names{i}).MaxLim + ParametersTilted.(Names{i}).MinLim)/2;
        ParametersTilted = DefineEstimatedParametersIndexes(ParametersTilted);
        ParametersTilted = DefineTransfFunctions(ParametersTilted);
        ParametersTilted = DefinePriors(ParametersTilted);
        ParametersTilted = UpdateParsNoTransfToTransf(ParametersTilted);
        Ress{ind,1} = QuickHIVfilteredPost(Data,ParametersTilted,HIVModel);

        % minus
        ParametersTilted.(Names{i}).MaxLim = max(0.2*mean,Parameters.(Names{i}).MaxLim - lambda*ampl);
        ParametersTilted.(Names{i}).MinLim = max(0.1*mean,Parameters.(Names{i}).MinLim -lambda*ampl);
        ParametersTilted.(Names{i}).Value = (ParametersTilted.(Names{i}).MaxLim + ParametersTilted.(Names{i}).MinLim)/2;
        ParametersTilted = DefineEstimatedParametersIndexes(ParametersTilted);
        ParametersTilted = DefineTransfFunctions(ParametersTilted);
        ParametersTilted = DefinePriors(ParametersTilted);
        ParametersTilted = UpdateParsNoTransfToTransf(ParametersTilted);
        Ress{ind,2} = QuickHIVfilteredPost(Data,ParametersTilted,HIVModel);
    else
        ParametersTilted = Parameters;
        ampl = Parameters.(Names{i}).Max - Parameters.(Names{i}).Min;
        mean = (Parameters.(Names{i}).Max + Parameters.(Names{i}).Min)/2;


        % plus
        if strcmp(Names{i},'eHIV')
            ParametersTilted.(Names{i}).Max = min(1,Parameters.(Names{i}).Max + lambda*ampl);
        else
            ParametersTilted.(Names{i}).Max = Parameters.(Names{i}).Max + lambda*ampl;
        end
        ParametersTilted.(Names{i}).Min = max(0.1*mean,Parameters.(Names{i}).Min + lambda*ampl);

        ParametersTilted.(Names{i}).MaxLim = max(ParametersTilted.(Names{i}).MaxLim,ParametersTilted.(Names{i}).Max);       
        ParametersTilted.(Names{i}).MinLim = max(ParametersTilted.(Names{i}).MinLim,ParametersTilted.(Names{i}).Min);
        ParametersTilted.(Names{i}).Value =(ParametersTilted.(Names{i}).MaxLim + ParametersTilted.(Names{i}).MinLim)/2;
        ParametersTilted = DefineEstimatedParametersIndexes(ParametersTilted);
        ParametersTilted = DefineTransfFunctions(ParametersTilted);
        ParametersTilted = DefinePriors(ParametersTilted);
        ParametersTilted = UpdateParsNoTransfToTransf(ParametersTilted);
        Ress{ind,1} = QuickHIVfilteredPost(Data,ParametersTilted,HIVModel);

        % minus
        ParametersTilted.(Names{i}).Max = max(0.2*mean,Parameters.(Names{i}).Max - lambda*ampl);
        ParametersTilted.(Names{i}).Min = max(0.1*mean,Parameters.(Names{i}).Min -lambda*ampl);
        ParametersTilted.(Names{i}).MaxLim = max(ParametersTilted.(Names{i}).MaxLim,ParametersTilted.(Names{i}).Max);       
        ParametersTilted.(Names{i}).MinLim = min(ParametersTilted.(Names{i}).MinLim,ParametersTilted.(Names{i}).Min);
        ParametersTilted.(Names{i}).Value = ParametersTilted.(Names{i}).Max ;
        ParametersTilted = DefineEstimatedParametersIndexes(ParametersTilted);
        ParametersTilted = DefineTransfFunctions(ParametersTilted);
        ParametersTilted = DefinePriors(ParametersTilted);
        ParametersTilted = UpdateParsNoTransfToTransf(ParametersTilted);
        Ress{ind,2} = QuickHIVfilteredPost(Data,ParametersTilted,HIVModel);
    end
end

PlotResHIV(Ress{7,2})

for i = 1:length(Names)
    ind = Parameters.(Names{i}).Index;
    Ress{ind,3} = Res;
end


SavePath = 'S:\Results';
% SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
% save([SavePath '\HIV_Mysore_ForProj_02I0_Calib.mat'],'ResRW1_RW_Calib_Rds_1_2_3')
save([SavePath '\HIV_Mysore_PlayingPriors.mat'],'Ress')

load([SavePath '\HIV_Mysore_PlayingPriors.mat'])

% add the no-transformation in blue

clear mean
clf

indsmooth = 15;
for i = 1:length(Names)
    subplot(4,3,i)
    ind = Parameters.(Names{i}).Index;
    ParametersTilted1 = Ress{ind,1}.Parameters;
    ParametersTilted2 = Ress{ind,2}.Parameters;
    Parameters = Ress{ind,3}.Parameters;
    vals = {};
    plot(smooth(mean(squeeze(Ress{ind,1}.Paths(:,3,:))),indsmooth),'g')
    hold on
    plot(smooth(mean(squeeze(Ress{ind,3}.Paths(:,3,:))),indsmooth),'LineWidth',2)
    plot(smooth(mean(squeeze(Ress{ind,2}.Paths(:,3,:))),indsmooth),'r')
    plot(smooth(mean(squeeze(Ress{ind,1}.Paths(:,3,:))),indsmooth),'w','LineWidth',2)
    plot(smooth(mean(squeeze(Ress{ind,2}.Paths(:,3,:))),indsmooth),'w','LineWidth',2)
    plot(smooth(quantile(squeeze(Ress{ind,3}.Paths(:,3,:)),0.75),indsmooth),':','LineWidth',2)
    plot(smooth(quantile(squeeze(Ress{ind,3}.Paths(:,3,:)),0.25),indsmooth),':','LineWidth',2)
    for j = 1:3
        if Ress{ind,j}.LogLik<-20
            vals{j} = 'N.C.';
            disp(Ress{ind,j}.LogLik)
        else
            if j == 1
                vals{j} = num2str(ParametersTilted1.(Names{i}).Value,2);
                plot(smooth(mean(squeeze(Ress{ind,j}.Paths(:,3,:))),indsmooth),'g','LineWidth',2)
                plot(smooth(quantile(squeeze(Ress{ind,j}.Paths(:,3,:)),0.75),indsmooth),':g','LineWidth',2)
                plot(smooth(quantile(squeeze(Ress{ind,j}.Paths(:,3,:)),0.25),indsmooth),':g','LineWidth',2)
            elseif j == 2
                vals{j} = num2str(ParametersTilted2.(Names{i}).Value,2);
                plot(smooth(mean(squeeze(Ress{ind,j}.Paths(:,3,:))),indsmooth),'r','LineWidth',2)
                plot(smooth(quantile(squeeze(Ress{ind,j}.Paths(:,3,:)),0.75),indsmooth),':r','LineWidth',2)
                plot(smooth(quantile(squeeze(Ress{ind,j}.Paths(:,3,:)),0.25),indsmooth),':r','LineWidth',2)
            else
                vals{j} = num2str(Parameters.(Names{i}).Value,2);
            end
        end
    end
    hold off
    xlabel('time')
    ylabel('Condom Use Frequency')
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    title([Names{i} '  (' vals{2} ' - ' vals{3} ' - ' vals{1} ')'],'FontWeight','bold')
    legend('+','prior','-')
    xlim([0 600])
end




%% Belgaum

% Setting


Parameters = struct();

% Initializing all model parameters
Parameters.TotalFSW.Value = 1742;
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
Parameters.InitialIPropF.Value = 0.02;
Parameters.InitialIPropF.Min = 0;
Parameters.InitialIPropF.Max = 0.04;
Parameters.InitialIPropF.MaxLim = 0.10;
Parameters.InitialIPropF.MinLim = 0;
Parameters.InitialIPropF.Estimated = 1;
Parameters.InitialIPropF.TransfType = 'Logit';
Parameters.InitialIPropM.Value = 0.001;
Parameters.InitialIPropM.Min = 0;
Parameters.InitialIPropM.Max = 0.02;
Parameters.InitialIPropM.MaxLim = 0.1;
Parameters.InitialIPropM.MinLim = 0;
Parameters.InitialIPropM.Estimated = 1;
Parameters.InitialIPropM.TransfType = 'Logit';
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotMFactor.Value = 12;
Parameters.TotMFactor.Min = -10^14;
Parameters.TotMFactor.Max = 10^14;
Parameters.TotMFactor.MinLim = 7;
Parameters.TotMFactor.MaxLim = 35;
Parameters.TotMFactor.Estimated = 1;
Parameters.TotMFactor.TransfType = 'Logit';
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.Alpham1.Value = 110;
Parameters.Alpham1.Min = -10^14;
Parameters.Alpham1.Max = 10^14;
Parameters.Alpham1.MinLim = 92;
Parameters.Alpham1.MaxLim = 138.5;
Parameters.Alpham1.Estimated = 1;
Parameters.Alpham1.TransfType = 'Logit';
Parameters.Alpham1.PlotInv = 1;
Parameters.MuFm1.Value = 120;
Parameters.MuFm1.Min = -10^14;
Parameters.MuFm1.Max =  10^14;
Parameters.MuFm1.MinLim = 94; 
Parameters.MuFm1.MaxLim = 250;
Parameters.MuFm1.Estimated = 1;
Parameters.MuFm1.TransfType = 'Logit';
Parameters.MuFm1.PlotInv = 1;
Parameters.MuMm1.Value = 120;
Parameters.MuMm1.Min = -10^14;
Parameters.MuMm1.Max = 10^14;
Parameters.MuMm1.MinLim = 83;
Parameters.MuMm1.MaxLim = 143;
Parameters.MuMm1.Estimated = 1;
Parameters.MuMm1.PlotInv = 1;
Parameters.MuMm1.TransfType = 'Logit';
Parameters.BetaMFPerAct.Value = 0.0029;
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
Parameters.CF1.Value = 23;
Parameters.CF1.Min = -10^14;
Parameters.CF1.Max =  10^14;
Parameters.CF1.MinLim = 22.1;
Parameters.CF1.MaxLim = 25.1;
Parameters.CF1.Estimated = 1;
Parameters.CF1.TransfType = 'Logit';
Parameters.CF2.Value = 95;
Parameters.CF2.Min = -10^14;
Parameters.CF2.Max = 10^14;
Parameters.CF2.MinLim = 84.8;
Parameters.CF2.MaxLim =  103.8;
Parameters.CF2.Estimated = 1;
Parameters.CF2.TransfType = 'Logit';
Parameters.SigmaRW.Value = 10;
Parameters.SigmaRW.Min = -10^14;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.MinLim = 0;
Parameters.SigmaRW.MaxLim = 50;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Logit';
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
Data.Observations(7,2) = 33.9;
Data.Observations(7,3) = 27.3;
Data.Observations(8,4) = 6.2;
Data.Observations(7,5) = 22.3;
Data.Instants = round([0 250 283 286 309]/(Parameters.ComputationTStep));
Data.ObservedVariables = [ 0 7 7 8 7];
Data.NbComputingSteps = [0 diff(Data.Instants)];


HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIVInitialize;
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*Parameters.SigmaObs).*(Res.WentOutOrNot)';
HIVModel.SMC_projection = @HIV_SMC_projection;

temp7 = zeros(1,10);
temp8 = zeros(1,10);
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



Names = Parameters.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters.DiffusionType = 'AffineAdd';

[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,HIVModel,Parameters),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));
Initialization = x;
Names = Parameters.Names.Estimated;
ParametersKalman = Parameters;
for i = 1:length(Names)
    ParametersKalman.(Names{i}).TransfValue = (x(i));
    Parameters.(Names{i}).TransfValue = (x(i));
end
ParametersKalman = UpdateParsTransfToNoTransf(ParametersKalman);
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(ParametersKalman)
Names = ParametersKalman.Names.Estimated;
for i = 1:length(Names)
    ParametersKalman.(Names{i}).SamplStd = 0.001*ParametersKalman.(Names{i}).Value;
end
ResKal = KalmanNumericDerivativesWithPrior(Data,HIVModel,ParametersKalman);
Test = mean(eig(-ResKal.Hess)>0)==1;
Cov = (-ResKal.Hess)^-1;
Test
% 
% Names = Parameters.Names.Estimated;
% Cov = zeros(length(Names),length(Names));
% for i = 1:length(Names)
%     if and(strcmp(Names{i}(end-1:end),'Ft'),not(strcmp(Names{i}(1:2),'In')))
%         Parameters.(Names{i}).Estimated = 0;
%     else
%         Parameters.(Names{i}).Estimated = 1;
%         Parameters.(Names{i}).Value = ParametersKalman.(Names{i}).Value;
%     end
% end
% Parameters.SigmaRW.Estimated = 1;
% Parameters.SigmaRW.Value = 0.1;
% Parameters = DefineEstimatedParametersIndexes(Parameters);
% Parameters = DefineTransfFunctions(Parameters);
% Parameters = DefinePriors(Parameters);
% Parameters = UpdateParsNoTransfToTransf(Parameters);
% Parameters.StableCUseConstraint = 0;
% Names = Parameters.Names.Estimated;
% Cov = zeros(length(Names),length(Names));
% 
% inds = [];
% indsKal = [];
% for i = 1:length(Names)
%     if ParametersKalman.(Names{i}).Estimated
%         indsKal(i) =  ParametersKalman.(Names{i}).Index;
%         inds(i) = Parameters.(Names{i}).Index;
%     end
% end
% tmp = (-ResKal.Hess)^-1;
% Cov(inds,inds) = tmp(indsKal,indsKal);
% indsig = Parameters.SigmaRW.Index;
% Cov(indsig,indsig) = 0.1;
%     
% 
% Parameters.MIFNbIterations = 100;
% Parameters.MIFNbParticules = 10000;
% Parameters.MIFCoolingParameters = 0.000005^(1/Parameters.MIFNbIterations);
% Parameters.MIFb = 1;
% Parameters.MIFSigmas = [];
% Parameters.NoPaths = 0;
% Parameters.PathsToKeep = [1:9]';
% Parameters.MCMCType = 'gtgt';
% Names = Parameters.Names.Estimated;
% Parameters.MIFCov = 0.0005*(Cov);
% 
% ResMIF = MIFHIV(Data,Parameters,HIVModel);
% 
% 
% Names = Parameters.Names.Estimated;
% for i = 1:length(Names)
%     ind = Parameters.(Names{i}).Index;
%     Parameters.(Names{i}).TransfValue = ResMIF.ThetasRecord(ind,end);
% end
% Parameters= UpdateParsTransfToNoTransf(Parameters);

Names = Parameters.Names.All;
for i = 1:length(Names);
    Parameters.(Names{i}).Estimated = 0;
end
Parameters.SigmaRW.Value = 0.1;
Parameters.SigmaRW.Min = 0;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Log';
Parameters.NbParticules = 1000;
Parameters.StableCUseConstraint = 0;
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [7 8 9];
Parameters.DiffusionType = 'Add';
Parameters.MCMCType = 'jiji';
Parameters.GMeth = 'cst given';
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
Parameters.NbParticules = 1000;

[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizewithPrior(x,Data,Model,Parameters),Initialization,optimset('MaxIter',50,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)




Res = EstimationSMCfiltGen(Data, Model, Parameters);
Res.LogLik
PlotResHIV(Res)

NbIts = 20;
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

PlotResHIV(Res)


%% Playing with priors, one at a time.

Parameters.Problem = 'ImperialHIV';
Namestmp = Parameters.Names.Estimated;
Names = {};
for i = 7:length(Namestmp)-1
    Names{i-6} = Namestmp{i};
end

ParametersTilted = Parameters;

Ress = {};
lambda = 1;
for i = 1:length(Names)
    ind = Parameters.(Names{i}).Index;
    
    if Parameters.(Names{i}).Max == 10^14
        ParametersTilted = Parameters;
        ampl = Parameters.(Names{i}).MaxLim - Parameters.(Names{i}).MinLim;
        mean = (Parameters.(Names{i}).MaxLim + Parameters.(Names{i}).MinLim)/2;


        % plus
        if strcmp(Names{i},'eHIV')
            ParametersTilted.(Names{i}).MaxLim = min(1,Parameters.(Names{i}).MaxLim + lambda*ampl);
        else
            ParametersTilted.(Names{i}).MaxLim = Parameters.(Names{i}).MaxLim + lambda*ampl;
        end
        ParametersTilted.(Names{i}).MinLim = max(0.1*mean,Parameters.(Names{i}).MinLim + lambda*ampl);

        ParametersTilted.(Names{i}).Value = (ParametersTilted.(Names{i}).MaxLim + ParametersTilted.(Names{i}).MinLim)/2;
        ParametersTilted = DefineEstimatedParametersIndexes(ParametersTilted);
        ParametersTilted = DefineTransfFunctions(ParametersTilted);
        ParametersTilted = DefinePriors(ParametersTilted);
        ParametersTilted = UpdateParsNoTransfToTransf(ParametersTilted);
        Ress{ind,1} = QuickHIVfilteredPost(Data,ParametersTilted,HIVModel);

        % minus
        ParametersTilted.(Names{i}).MaxLim = max(0.2*mean,Parameters.(Names{i}).MaxLim - lambda*ampl);
        ParametersTilted.(Names{i}).MinLim = max(0.1*mean,Parameters.(Names{i}).MinLim -lambda*ampl);
        ParametersTilted.(Names{i}).Value = (ParametersTilted.(Names{i}).MaxLim + ParametersTilted.(Names{i}).MinLim)/2;
        ParametersTilted = DefineEstimatedParametersIndexes(ParametersTilted);
        ParametersTilted = DefineTransfFunctions(ParametersTilted);
        ParametersTilted = DefinePriors(ParametersTilted);
        ParametersTilted = UpdateParsNoTransfToTransf(ParametersTilted);
        Ress{ind,2} = QuickHIVfilteredPost(Data,ParametersTilted,HIVModel);
    else
        ParametersTilted = Parameters;
        ampl = Parameters.(Names{i}).Max - Parameters.(Names{i}).Min;
        mean = (Parameters.(Names{i}).Max + Parameters.(Names{i}).Min)/2;


        % plus
        if strcmp(Names{i},'eHIV')
            ParametersTilted.(Names{i}).Max = min(1,Parameters.(Names{i}).Max + lambda*ampl);
        else
            ParametersTilted.(Names{i}).Max = Parameters.(Names{i}).Max + lambda*ampl;
        end
        ParametersTilted.(Names{i}).Min = max(0.1*mean,Parameters.(Names{i}).Min + lambda*ampl);

        ParametersTilted.(Names{i}).MaxLim = max(ParametersTilted.(Names{i}).MaxLim,ParametersTilted.(Names{i}).Max);       
        ParametersTilted.(Names{i}).MinLim = max(ParametersTilted.(Names{i}).MinLim,ParametersTilted.(Names{i}).Min);
        ParametersTilted.(Names{i}).Value = max((ParametersTilted.(Names{i}).Max+ParametersTilted.(Names{i}).Min)/2,ParametersTilted.(Names{i}).Min + 0.1*ampl);
        ParametersTilted = DefineEstimatedParametersIndexes(ParametersTilted);
        ParametersTilted = DefineTransfFunctions(ParametersTilted);
        ParametersTilted = DefinePriors(ParametersTilted);
        ParametersTilted = UpdateParsNoTransfToTransf(ParametersTilted);
        Ress{ind,1} = QuickHIVfilteredPost(Data,ParametersTilted,HIVModel);

        % minus
        ParametersTilted.(Names{i}).Max = max(0.2*mean,Parameters.(Names{i}).Max - lambda*ampl);
        ParametersTilted.(Names{i}).Min = max(0.1*mean,Parameters.(Names{i}).Min -lambda*ampl);
        ParametersTilted.(Names{i}).MaxLim = max(ParametersTilted.(Names{i}).MaxLim,ParametersTilted.(Names{i}).Max);       
        ParametersTilted.(Names{i}).MinLim = min(ParametersTilted.(Names{i}).MinLim,ParametersTilted.(Names{i}).Min);
        ParametersTilted.(Names{i}).Value = max((ParametersTilted.(Names{i}).Max+ParametersTilted.(Names{i}).Min)/2,ParametersTilted.(Names{i}).Max - 0.1*ampl);
        ParametersTilted = DefineEstimatedParametersIndexes(ParametersTilted);
        ParametersTilted = DefineTransfFunctions(ParametersTilted);
        ParametersTilted = DefinePriors(ParametersTilted);
        ParametersTilted = UpdateParsNoTransfToTransf(ParametersTilted);
        Ress{ind,2} = QuickHIVfilteredPost(Data,ParametersTilted,HIVModel);
    end
end

PlotResHIV(Ress{7,2})

for i = 1:length(Names)
    ind = Parameters.(Names{i}).Index;
    Ress{ind,3} = Res;
end

PlotResHIV(Ress{7,2})




SavePath = 'S:\Results';
% SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
% save([SavePath '\HIV_Mysore_ForProj_02I0_Calib.mat'],'ResRW1_RW_Calib_Rds_1_2_3')
save([SavePath '\HIV_Belgaum_PlayingPriors2.mat'],'Ress')

load([SavePath '\HIV_Belgaum_PlayingPriors2.mat'])

clear mean
clf

indsmooth = 15;
for i = 1:length(Names)
    subplot(4,3,i)
    ind = Parameters.(Names{i}).Index;
    ParametersTilted1 = Ress{ind,1}.Parameters;
    ParametersTilted2 = Ress{ind,2}.Parameters;
    Parameters = Ress{ind,3}.Parameters;
    vals = {};
    plot(smooth(mean(squeeze(Ress{ind,1}.Paths(:,3,:))),indsmooth),'g')
    hold on
    plot(smooth(mean(squeeze(Ress{ind,3}.Paths(:,3,:))),indsmooth),'LineWidth',2)
    plot(smooth(mean(squeeze(Ress{ind,2}.Paths(:,3,:))),indsmooth),'r')
    plot(smooth(mean(squeeze(Ress{ind,1}.Paths(:,3,:))),indsmooth),'w','LineWidth',2)
    plot(smooth(mean(squeeze(Ress{ind,2}.Paths(:,3,:))),indsmooth),'w','LineWidth',2)
    plot(smooth(quantile(squeeze(Ress{ind,3}.Paths(:,3,:)),0.75),indsmooth),':','LineWidth',2)
    plot(smooth(quantile(squeeze(Ress{ind,3}.Paths(:,3,:)),0.25),indsmooth),':','LineWidth',2)
    for j = 1:3
        if Ress{ind,j}.LogLik<-20
            vals{j} = 'N.C.';
            disp(Ress{ind,j}.LogLik)
        else
            if j == 1
                vals{j} = num2str(ParametersTilted1.(Names{i}).Value,2);
                plot(smooth(mean(squeeze(Ress{ind,j}.Paths(:,3,:))),indsmooth),'g','LineWidth',2)
                plot(smooth(quantile(squeeze(Ress{ind,j}.Paths(:,3,:)),0.75),indsmooth),':g','LineWidth',2)
                plot(smooth(quantile(squeeze(Ress{ind,j}.Paths(:,3,:)),0.25),indsmooth),':g','LineWidth',2)
            elseif j == 2
                vals{j} = num2str(ParametersTilted2.(Names{i}).Value,2);
                plot(smooth(mean(squeeze(Ress{ind,j}.Paths(:,3,:))),indsmooth),'r','LineWidth',2)
                plot(smooth(quantile(squeeze(Ress{ind,j}.Paths(:,3,:)),0.75),indsmooth),':r','LineWidth',2)
                plot(smooth(quantile(squeeze(Ress{ind,j}.Paths(:,3,:)),0.25),indsmooth),':r','LineWidth',2)
            else
                vals{j} = num2str(Parameters.(Names{i}).Value,2);
            end
        end
    end
    hold off
    xlabel('time')
    ylabel('Condom Use Frequency')
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    title([Names{i} '  (' vals{2} ' - ' vals{3} ' - ' vals{1} ')'],'FontWeight','bold')
    legend('+','prior','-')
    xlim([0 600])
end

