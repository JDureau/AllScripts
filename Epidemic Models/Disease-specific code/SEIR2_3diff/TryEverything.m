function Res = TryEverything(DataTemp,SEIRModel,Parameters)

% Adapt priors
Parameters.DiffusionType = 'Add';

Parameters.TotalPopulation = 1;
Parameters.EInitProp.Value = max(eps,1.5*DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.MinLim = 0.00001*max(eps,DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.MaxLim = 10*max(eps,DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Value = max(eps,1.5*DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.MinLim = 0.00001*max(eps,DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.MaxLim = 10*max(eps,DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.betainit.Value = 0.6;
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.DiffusionType = 'Add';
Parameters = SEIRModel.InitializeParameters(Parameters);

Parameters.InitialCovFact.Value = 0.002;
Parameters.InitialCovFact.Min = -10^14;
Parameters.InitialCovFact.Max =  10^14;
Parameters.InitialCovFact.Estimated =  0;
Parameters.InitialCovFact.TransfType = 'Log';
Parameters.InitialCovFact.Init = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


% Run EKF Add 1 (0.05) & 2 (0.00)
temp = zeros(1,7);
temp(1,5) = 1;
SEIRModel.ObservationJacobian = {};
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(DataTemp.Instants)
    SEIRModel.ObservationJacobian{i} = temp;
    SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs*DataTemp.Observations(5,i))^2;
end
SEIRModel.InitializeParameters = @SEIRInitialize;

nbtests = 2000;
cpt = 0;
test = 0;
while not(test)
    Parameters = SampleParameters(Parameters);
    cpt = cpt+1;
    try
        Tmp = EstimationEKFGen(DataTemp, SEIRModel, Parameters);
        if (Tmp.LogLik)>200
            test = 1;
        end
    end
    if cpt>nbtests
        disp('doesn''t work with this trajectory')
        Res = struct();
        return
    end
end
Parameters = SEIRModel.InitializeParameters(Parameters);
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
fval = Inf;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,DataTemp,SEIRModel,Parameters),Initialization,optimset('MaxIter',3000,'TolX',1e-8,'TolFun',1e-3));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)
% Parameters = SEIRModel.InitializeParameters(Parameters);
% for i = 1:length(Names)
%     Initialization(i) = Parameters.(Names{i}).TransfValue ;
% end
% fval = Inf;
% Parameters.InitialCovFact = 0.00;
% [x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,DataTemp,SEIRModel,Parameters),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-3));
% Names = Parameters.Names.Estimated;
% for i = 1:length(Names)
%     Parameters.(Names{i}).TransfValue = (x(i));
% end
% Parameters = UpdateParsTransfToNoTransf(Parameters);
% TellParsValues(Parameters)

% Generate data
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:7]';
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr'
tempSMCsmooth = EstimationSMCsmoothGen(DataTemp, SEIRModel, Parameters);
yis3 = squeeze(exp(tempSMCsmooth.CompletePaths(1,6,:)));

Parameters.ObsNoise = 0;
DataGen = SEIR_CreateData(yis3,Parameters,DataTemp, SEIRModel);
ParametersGen = Parameters;
DataGen.RealBetaTraj = yis3;

Res.DataTemp = DataTemp;
Res.DataGen = DataGen;
Res.ParametersGen = Parameters;
if not(max(DataTemp.Observations(5,:))<0.01)
   return 
end

%% Apply Everything With Add
Parameters.DiffusionType = 'Add';

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).Estimated = 0;
end
Parameters.SigmaRW.Estimated = 1;
Parameters.InitialCovFact.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

% 1 : filt EKF
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
fval = Inf;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,DataGen,SEIRModel,Parameters),Initialization,optimset('MaxIter',3000,'TolX',1e-8,'TolFun',1e-5));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)

Res.ResEKF_Add =  EstimationEKFGen(DataGen, SEIRModel, Parameters);

% 2 : filt SMC
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
fval = Inf;
Parameters.NoPaths = 1;
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,DataTemp,SEIRModel,Parameters),Initialization,optimset('MaxIter',100,'TolX',1e-8,'TolFun',1e-5));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)

Parameters.NoPaths = 0;
Res.ResSMCfilt_Add = EstimationSMCfiltGen(DataGen, SEIRModel, Parameters);
Res.ResSMCfilt_Add.Paths = Res.ResSMCfilt_Add.Paths(ceil(rand(1,50)*Parameters.NbParticules),:,:);
Res.ResSMCfilt_Add.ResampledParticules = [];

% 3 : smooth SMC
Parameters.NbParticules = 1000;
TempSMCsmooth = struct();
Nbits = 60;
TempSMCsmooth.Paths = []; 
for i = 1:Nbits
    disp(['smooth ' num2str(i)]) 
    tempSMCsmooth = EstimationSMCsmoothGen(DataGen, SEIRModel, Parameters);
    ind = ceil(rand(1,1)*Parameters.NbParticules);
    TempSMCsmooth.Paths(i,:,:) = tempSMCsmooth.CompletePaths(ind,:,:);
end
Res.ResSMCsmooth_Add = TempSMCsmooth;

% TempSMCsmooth.Data = DataGen;
% TempSMCsmooth.Parameters = Parameters;
% PlotResGoogle(TempSMCsmooth,8)

%% Apply Everything With IBM
Parameters.DiffusionType = 'IBM';

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).Estimated = 0;
end
Parameters.SigmaRW.Estimated = 1;
Parameters.InitialCovFact.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Parameters.InitialCovFact.Value = 0.20;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


try 
    
    
    tmpLogLik = -Inf;
    nbtests = 10000;
    cpt = 0;
    test = 0;
    while not(test)
        Parameters.InitialCovFact.Value = 0.05 + 0.15*rand(1,1);
        Parameters.SigmaRW.Value = 2*rand(1,1);
        Parameters = UpdateParsNoTransfToTransf(Parameters);
        cpt = cpt+1;
        try
            Tmp = EstimationEKFGen(DataGen, SEIRModel, Parameters);
            if (Tmp.LogLik)>200
                test = 1;
            end
        end
        if cpt>nbtests
            disp('doesn''t work with this trajectory')
            Res = struct();
            return
        end
        if (Tmp.LogLik)>tmpLogLik
            tmpLogLik = Tmp.LogLik;
            BestPars = Parameters;
        end
        disp(cpt)
        disp(Tmp.LogLik)
    end
    
    Parameters = BestPars;
    
    % 1 : filt EKF
    Names = Parameters.Names.Estimated;
    for i = 1:length(Names)
        Initialization(i) = Parameters.(Names{i}).TransfValue ;
    end
    fval = Inf;
    [x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,DataGen,SEIRModel,Parameters),Initialization,optimset('MaxIter',10000,'TolX',1e-8,'TolFun',1e-4));
    Names = Parameters.Names.Estimated;
    for i = 1:length(Names)
        Parameters.(Names{i}).TransfValue = (x(i));
    end
    Parameters = UpdateParsTransfToNoTransf(Parameters);
    TellParsValues(Parameters)

    Res.ResEKF_IBM =  EstimationEKFGen(DataGen, SEIRModel, Parameters);



    % 2 : filt SMC
    Parameters.InitialCovFact.Estimated = 0;
    Parameters = DefineEstimatedParametersIndexes(Parameters);
    Parameters = DefineTransfFunctions(Parameters);
    Parameters = DefinePriors(Parameters);
    Parameters = UpdateParsNoTransfToTransf(Parameters);
    Names = Parameters.Names.Estimated;
    for i = 1:length(Names)
        Initialization(i) = Parameters.(Names{i}).TransfValue ;
    end
    fval = Inf;
    Parameters.NoPaths = 1;
    [x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,DataGen,SEIRModel,Parameters),Initialization,optimset('MaxIter',100,'TolX',1e-8,'TolFun',1e-3));
    Names = Parameters.Names.Estimated;
    for i = 1:length(Names)
        Parameters.(Names{i}).TransfValue = (x(i));
    end
    Parameters = UpdateParsTransfToNoTransf(Parameters);
    TellParsValues(Parameters)

    Parameters.NoPaths = 0;

    Res.ResSMCfilt_IBM = EstimationSMCfiltGen(DataGen, SEIRModel, Parameters);
    Res.ResSMCfilt_IBM.Paths = Res.ResSMCfilt_IBM.Paths(ceil(rand(1,50)*Parameters.NbParticules),:,:);
    Res.ResSMCfilt_IBM.ResampledParticules = [];

    % 3 : smooth SMC
    Parameters.NbParticules = 1000;
    TempSMCsmooth = struct();
    Nbits = 60;
    TempSMCsmooth.Paths = []; 
    for i = 1:Nbits
        disp(['smooth ' num2str(i)]) 
        tempSMCsmooth = EstimationSMCsmoothGen(DataGen, SEIRModel, Parameters);
        ind = ceil(rand(1,1)*Parameters.NbParticules);
        TempSMCsmooth.Paths(i,:,:) = tempSMCsmooth.CompletePaths(ind,:,:);
    end
    Res.ResSMCsmooth_IBM = TempSMCsmooth;
    
%     TempSMCsmooth.Data = DataGen;
%     TempSMCsmooth.Parameters = Parameters;
%     PlotResGoogle(TempSMCsmooth,8)
end
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

