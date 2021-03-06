function Res = SEIRinferenceSim(Data)

% provides estimates through StoPF, StoPS, DetPF, DetPS, DetEKF
Res = struct();


Parameters = Data.Parameters;

Parameters.ComputationTStep = 0.1;
Data.Instants = [0:size(Data.Observations,2)-1]*7/Parameters.ComputationTStep;



% Initialize Model

SEIRModel = struct();
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),coeff*Data.Observations(5,IndTime),coeff*Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
SEIRModel.SMC_projection = @SEIR_SMCsto_projection;

Parameters = SEIRModel.InitializeParameters(Parameters);



Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).Estimated = 0;
end
%     Parameters.betainit.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;

Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
% Parameters.InitialCov.Estimated = 1;


Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;

Initialization = [];
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
Parameters.Correction = 0;
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',50,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters) 

ParametersKal = Parameters;

ResKal = EstimationEKFGen(Data, SEIRModel, Parameters);

Res.EKF = ResKal;


Parameters.ComputationTStep = 0.1;
Data.Instants = [0:size(Data.Observations,2)-1]*7/Parameters.ComputationTStep;

Parameters.NoPaths = 1;
Parameters.PathsToKeep = [1:7]';
Parameters.NbParticules = 2000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;

Initialization = [];
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
Parameters.Correction = 0;
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',50,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  

Res.DisPF = EstimationSMCFiltGen(Data, SEIRModel, Parameters);

SMCPaths = [];
NbIts = 40; % I check this was enough for convergence on an example, both for MSE and estimate CI.
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:7]';
Parameters.NbParticules = 1000;
for i = 1:NbIts
    i
    Temp = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);
    ind = ceil(rand(1,1)*Parameters.NbParticules);
    SMCPaths(end+1,:) = Temp.CompletePaths(ind,6,:);
end


Res.DisPS = SMCPaths;




