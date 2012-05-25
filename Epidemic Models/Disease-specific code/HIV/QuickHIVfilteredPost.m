function Res = QuickHIVfilteredPost(Data,Parameters,Model)
tic



Parameters.MIFNbIterations = 100;
Parameters.MIFNbParticules = 10000;
Parameters.MIFCoolingParameters = 0.000005^(1/Parameters.MIFNbIterations);
Parameters.MIFb = 1;
Parameters.MIFSigmas = [];
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:9]';
Parameters.MCMCType = 'gtgt';
Names = Parameters.Names.Estimated;
Parameters.MIFCov = 0.0005*(Cov);

ResMIF = MIFHIV(Data,Parameters,HIVModel);



% first, model parameters are fixed to MLE, we just estimate the diffusion
% parameters
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
Res = EstimationSMCsmoothGen(Data, Model, Parameters);

Parameters.InitialFt.Value = min(Parameters.InitialFt.MaxLim,max(Parameters.InitialFt.MinLim,Res.PosteriorMeansRecord(9,2)));
Parameters.InitialFt.Estimated = 1;
Parameters.SecondFt.Value = min(Parameters.SecondFt.MaxLim,max(Parameters.InitialFt.MinLim,Res.PosteriorMeansRecord(9,2)));
Parameters.SecondFt.Estimated = 1;
Parameters.ThirdFt.Value = min(Parameters.ThirdFt.MaxLim,max(Parameters.InitialFt.MinLim,Res.PosteriorMeansRecord(9,3)));
Parameters.ThirdFt.Estimated = 1;
Parameters.FourthFt.Value = min(Parameters.FourthFt.MaxLim,max(Parameters.InitialFt.MinLim,Res.PosteriorMeansRecord(9,3)));
Parameters.FourthFt.Estimated = 1;
Parameters.FifthFt.Value = min(Parameters.FifthFt.MaxLim,max(Parameters.InitialFt.MinLim,Res.PosteriorMeansRecord(9,4)));
Parameters.FifthFt.Estimated = 1;
Parameters.SixthFt.Value = min(Parameters.SixthFt.MaxLim,max(Parameters.InitialFt.MinLim,Res.PosteriorMeansRecord(9,5)));
Parameters.SixthFt.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

%% Kalman Opt
disp('Step 1: Kalman Opt')

Names = Parameters.Names.All;
Initialization = [];
for i = 1:length(Names)
    Parameters.(Names{i}).Estimated = 1;
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
Parameters.SigmaRW.Estimated = 0;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters.DiffusionType = 'AffineAdd';

[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,Model,Parameters),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-5,'MaxFunEvals',10000));
Initialization = x;
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)


%% MIF

disp('Step 2: SigmaRW')


% first, model parameters are fixed to MLE, we just estimate the diffusion
% parameters
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


% 
% 
% Names = Parameters.Names.Estimated;
% Cov = zeros(length(Names),length(Names));
% for i = 1:length(Names)
%     if and(strcmp(Names{i}(end-1:end),'Ft'),not(strcmp(Names{i}(1:2),'In')))
%         Parameters.(Names{i}).Estimated = 0;
%     else
%         Cov(i,i) = 0.1*Parameters.(Names{i}).TransfValue;
%     end
% end
% Parameters = DefineEstimatedParametersIndexes(Parameters);
% Parameters = DefineTransfFunctions(Parameters);
% Parameters = DefinePriors(Parameters);
% Parameters = UpdateParsNoTransfToTransf(Parameters);
% 
% Parameters.MIFNbIterations = 50;
% Parameters.MIFNbParticules = 5000;
% Parameters.MIFCoolingParameters = 0.000005^(1/Parameters.MIFNbIterations);
% Parameters.MIFb = 1;
% Parameters.MIFSigmas = [];
% Parameters.NoPaths = 0;
% Parameters.PathsToKeep = [1:9]';
% Parameters.MCMCType = 'gtgt';
% Names = Parameters.Names.Estimated;
% 
% ResMIF2 = MIFHIV(Data,Parameters,Model);
% 
% 
% Names = Parameters.Names.Estimated;
% for i = 1:length(Names)
%     ind = Parameters.(Names{i}).Index;
%     Parameters.(Names{i}).TransfValue = ResMIF.ThetasRecord(ind,end);
% end
% Parameters= UpdateParsTransfToNoTransf(Parameters);
% 
% Res = EstimationSMCfiltGen(Data, Model, Parameters);
% Res.LogLik
% Res.Paths = [];
% Res.Paths = ResMIF.RecordStates;
% PlotResHIV(Res)


%% Filter

disp('Step 3: Sigma Opt')

Parameters.NbParticules = 5000;
Parameters.DiffusionType = 'Add';
Res = EstimationSMCfiltGen(Data, Model, Parameters);
Res.LogLik

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


Res.ComputTime = toc;
Res.Data = Data;








