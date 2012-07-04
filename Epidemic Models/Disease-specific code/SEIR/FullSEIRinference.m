function [] = FullSEIRinference(Data,DiffType,ObsType,Name,IndModel)


% IndModel = 1 -> normal
% IndModel = 2 -> 2diff
% IndModel = 3 -> 2diffb
% IndModel = 4 -> 3diff
        
NbIters = 15000;
NbItersPrep = 5000;



SavePath = '/users/ecologie/dureau/src/AllData/ResultsMarc/';

switch IndModel
    case 1
        tmp = load([SavePath '/ParametersSEIR.mat']);
    case 2
        tmp = load([SavePath '/ParametersSEIR2_1diff.mat']);
    case 3
        tmp = load([SavePath '/ParametersSEIR2_1diffb.mat']);    
    case 4
        tmp = load([SavePath '/ParametersSEIR2_2diff.mat']);
    case 5
        tmp = load([SavePath '/ParametersSEIR2_2diffb.mat']);
    case 6
        tmp = load([SavePath '/ParametersSEIR2_3diff.mat']);
    case 7
        tmp = load([SavePath '/ParametersSEIR2_4diff.mat']);
end



switch IndModel
    case 1
        NameToSave = ['MarcData_StructModel_normal.mat'];
    case 2
        NameToSave = ['MarcData_StructModel_1diff.mat'];
    case 3
        NameToSave = ['MarcData_StructModel_1diffb.mat'];
    case 4
        NameToSave = ['MarcData_StructModel_2diff.mat'];
    case 5
        NameToSave = ['MarcData_StructModel_2diffb.mat'];
    case 6
        NameToSave = ['MarcData_StructModel_3diff.mat'];
    case 7
        NameToSave = ['MarcData_StructModel_4diff.mat'];
end

Parameters = tmp.Parameters;

Data.Instants = [1:size(Data.Observations,2)]*size(Data.Observations,1)/Parameters.ComputationTStep;
if IndModel>=2
    Data.ObservedVariables = diag([9 10])*ones(2,length(Data.Instants));
else
    Data.ObservedVariables = 5*ones(1,length(Data.Instants));
end
Data.NbComputingSteps = [0 diff(Data.Instants)];

if strcmp(ObsType,'Fixed')
    Parameters.SigmaObs.Value = 0.1;
elseif strcmp(ObsType,'Estimated')
    Parameters.SigmaObs.Value = sqrt(log(1.01));
    Parameters.SigmaObs.Min = -10^14;
    Parameters.SigmaObs.Max =  10^14;
    Parameters.SigmaObs.MinLim = sqrt(log(1.0001));
    Parameters.SigmaObs.MaxLim = sqrt(log(1.16));
    Parameters.SigmaObs.Estimated = 1;
    Parameters.SigmaObs.TransfType = 'Logit';
    Parameters.SigmaObs.Sample = 1;
elseif strcmp(ObsType,'CoeffStudy')
    Parameters.SigmaObs.Value = 0.1;
    Parameters.MultCoeff.Value = Data.Coeff;
    Parameters.MultCoeff.Min =  -10^14;
    Parameters.MultCoeff.Max =   10^14;
    Parameters.MultCoeff.MinLim = 5;
    Parameters.MultCoeff.MaxLim = 50;
    Parameters.MultCoeff.Estimated = 0;
    Parameters.MultCoeff.TransfType = 'Logit';
    Parameters.MultCoeff.Sample = 1;
end
Parameters.Correction = 1;

% load([SavePath '/Temp_' NameToSave])

try
    load([SavePath '/Temp_' NameToSave])
    AlreadySomething = 1;
catch
    AlreadySomething = 0;
end
 AlreadySomething 

SEIRModel = struct();

switch IndModel
    case 1
        SEIRModel.EKF_projection = @SEIR_EKF_projection;
        SEIRModel.InitializeParameters = @SEIRInitialize;
        SEIRModel.SMC_projection = @SEIR_SMC_projection;
        SEIRModel.LikFunction = 'normpdf(log(Variables(:,5)),transpose(log(coeff*Data.Observations(5,IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),sqrt(log(Parameters.SigmaObs.Value^2+1)))';%Parameters.SigmaObs.Value)';
    case 2
        SEIRModel.EKF_projection = @SEIR2_1diff_EKF_projection;
        SEIRModel.InitializeParameters = @SEIR2_1diff_Initialize;
        SEIRModel.SMC_projection = @SEIR2_1diff_SMC_projection;
        SEIRModel.LikFunction = 'mvnpdf(log(Variables(:,[9 10])),transpose(log(coeff*Data.Observations([9 10],IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),diag([log(Parameters.SigmaObs.Value^2+1) log(Parameters.SigmaObs.Value^2+1)]))';
    case 3
        SEIRModel.EKF_projection = @SEIR2_1diffb_EKF_projection;
        SEIRModel.InitializeParameters = @SEIR2_1diffb_Initialize;
        SEIRModel.SMC_projection = @SEIR2_1diffb_SMC_projection;
        SEIRModel.LikFunction = 'mvnpdf(log(Variables(:,[9 10])),transpose(log(coeff*Data.Observations([9 10],IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),diag([log(Parameters.SigmaObs.Value^2+1) log(Parameters.SigmaObs.Value^2+1)]))';
    case 4
        SEIRModel.EKF_projection = @SEIR2_2diff_EKF_projection;
        SEIRModel.InitializeParameters = @SEIR2_2diff_Initialize;
        SEIRModel.SMC_projection = @SEIR2_2diff_SMC_projection;
        SEIRModel.LikFunction = 'mvnpdf(log(Variables(:,[9 10])),transpose(log(coeff*Data.Observations([9 10],IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),diag([log(Parameters.SigmaObs.Value^2+1) log(Parameters.SigmaObs.Value^2+1)]))';
    case 5
        SEIRModel.EKF_projection = @SEIR2_2diffb_EKF_projection;
        SEIRModel.InitializeParameters = @SEIR2_2diffb_Initialize;
        SEIRModel.SMC_projection = @SEIR2_2diffb_SMC_projection;
        SEIRModel.LikFunction = 'mvnpdf(log(Variables(:,[9 10])),transpose(log(coeff*Data.Observations([9 10],IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),diag([log(Parameters.SigmaObs.Value^2+1) log(Parameters.SigmaObs.Value^2+1)]))';
    case 6
        SEIRModel.EKF_projection = @SEIR2_3diff_EKF_projection;
        SEIRModel.InitializeParameters = @SEIR2_3diff_Initialize;
        SEIRModel.SMC_projection = @SEIR2_3diff_SMC_projection;
        SEIRModel.LikFunction = 'mvnpdf(log(Variables(:,[9 10])),transpose(log(coeff*Data.Observations([9 10],IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),diag([log(Parameters.SigmaObs.Value^2+1) log(Parameters.SigmaObs.Value^2+1)]))';
    case 7
        SEIRModel.EKF_projection = @SEIR2_4diff_EKF_projection;
        SEIRModel.InitializeParameters = @SEIR2_4diff_Initialize;
        SEIRModel.SMC_projection = @SEIR2_4diff_SMC_projection;
        SEIRModel.LikFunction = 'mvnpdf(log(Variables(:,[9 10])),transpose(log(coeff*Data.Observations([9 10],IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),diag([log(Parameters.SigmaObs.Value^2+1) log(Parameters.SigmaObs.Value^2+1)]))';
end


% SEIRModel.LikFunction = 'normpdf(log(Variables(:,5)),log(coeff*Data.Observations(5,IndTime))-0.5*Parameters.SigmaObs.Value^2,Parameters.SigmaObs.Value)';
% SEIRModel.LikFunction = 'normpdf(Variables(:,5),coeff*Data.Observations(5,IndTime),coeff*Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
Parameters = SEIRModel.InitializeParameters(Parameters);

% temp = zeros(1,7);
% temp(1,5) = 1;
% SEIRModel.ObservationJacobian = {};
% SEIRModel.ObservationMeasurementNoise = {};
% for i = 1:length(Data.Instants)
%     SEIRModel.ObservationJacobian{i} = temp;
%     SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs.Value*Data.Observations(5,i))^2;
% end

if not(AlreadySomething)
    Names = Parameters.Names.Estimated;
    for j = 1:length(Names)
        Parameters.(Names{j}).Estimated = 0;
    end
    % Parameters.betainit.Estimated = 1;
    switch IndModel
        case 1
            Parameters.SigmaRW.Estimated = 1;
            Parameters.EInitProp.Estimated = 1;
            Parameters.IInitProp.Estimated = 1;
            Parameters.RInitProp.Estimated = 1;
        case 2
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
    %         Parameters.beta11init.Estimated = 1;
    %         Parameters.beta22init.Estimated = 1;
            Parameters.adultsmult.Estimated = 1;
            Parameters.kidsmult.Estimated = 1;
            Parameters.adultsadd.Estimated = 1;
            Parameters.kidsadd.Estimated = 1;
        case 3
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
    %         Parameters.beta11init.Estimated = 1;
    %         Parameters.beta22init.Estimated = 1;
            Parameters.adultsmult.Estimated = 1;
            Parameters.adultsadultsmult.Estimated = 1;
            Parameters.kidsmult.Estimated = 1;
            Parameters.adultsadd.Estimated = 1;
            Parameters.kidsadd.Estimated = 1;
            Parameters.adultsadultsadd.Estimated = 1;
        case 4
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.SigmaRW22.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
    %         Parameters.beta11init.Estimated = 1;
    %         Parameters.beta22init.Estimated = 1;
            Parameters.adultsmult.Estimated = 1;
        case 5
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.SigmaRW22.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
    %         Parameters.beta11init.Estimated = 1;
    %         Parameters.beta22init.Estimated = 1;
            Parameters.adultsmult.Estimated = 1;
            Parameters.kidsmult.Estimated = 1;
        case 6
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.SigmaRW22.Estimated = 1;
            Parameters.SigmaRW12.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
            Parameters.adultsmult.Estimated = 1;
        case 7
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.SigmaRW22.Estimated = 1;
            Parameters.SigmaRW12.Estimated = 1;
            Parameters.SigmaRW21.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
    end


    Parameters = DefineEstimatedParametersIndexes(Parameters);
    Parameters = DefineTransfFunctions(Parameters);
    Parameters = DefinePriors(Parameters);
    Parameters = UpdateParsNoTransfToTransf(Parameters);

    Test = 0;
    try
        Temp = EstimationEKFGen(Data, SEIRModel, Parameters);    
        if (Temp.LogLik>-600)
            Test = 1;
        end
    end
    NbIts = 0;
    while not(Test)
        Parameters = SampleParameters(Parameters);
        switch IndModel
            case 1
                Parameters.SigmaRW.Value = rand(1,1)*2;
            case 2
                Parameters.SigmaRW11.Value = 0.1+0.01*rand(1,1)*2;
                Parameters.adultsmult.Value = 1.7+0.1*rand(1,1)*2;
                Parameters.kidsmult.Value = 0.96+0.1*rand(1,1)*2;
                Parameters.adultsadd.Value = 0.19+0.1*rand(1,1)*2;
                Parameters.kidsadd.Value = 0.98+0.1*rand(1,1)*2;
                Parameters.R2InitProp.Value = 0.5;
                Parameters.R1InitProp.Value = 0.1;
            case 3
                Parameters.SigmaRW11.Value = rand(1,1)*2;
                Parameters.adultsmult.Value = rand(1,1)*2;
                Parameters.adultsadultsmult.Value = rand(1,1)*2;
                Parameters.kidsmult.Value = rand(1,1)*2;
                Parameters.adultsadd.Value = rand(1,1);
                Parameters.adultsadultsadd.Value = rand(1,1);
                Parameters.kidsadd.Value = rand(1,1);
                Parameters.R2InitProp.Value = 0.5;
                Parameters.R2InitProp.Value = 0.1;
            case 4
                Parameters.SigmaRW11.Value = rand(1,1)*2;
                Parameters.SigmaRW22.Value = rand(1,1)*2;
                Parameters.adultsmult.Value = rand(1,1);
            case 5
                Parameters.SigmaRW11.Value = rand(1,1)*2;
                Parameters.SigmaRW22.Value = rand(1,1)*2;
                Parameters.adultsmult.Value = rand(1,1);
                Parameters.kidsmult.Value = rand(1,1);
            case 6
                Parameters.SigmaRW11.Value = rand(1,1)*2;
                Parameters.SigmaRW22.Value = rand(1,1)*2;
                Parameters.SigmaRW12.Value = rand(1,1)*2;
                Parameters.adultsmult.Value = rand(1,1);
            case 7
                Parameters.SigmaRW11.Value = rand(1,1)*2;
                Parameters.SigmaRW22.Value = rand(1,1)*2;
                Parameters.SigmaRW12.Value = rand(1,1)*2;
                Parameters.SigmaRW21.Value = rand(1,1)*2;
        end
        Parameters = UpdateParsNoTransfToTransf(Parameters);
    %     Parameters.InitialCov = rand(1,1)*0.4;

        try
            Temp = EstimationEKFGen(Data, SEIRModel, Parameters);
            Temp.LogLik
            if (Temp.LogLik>-600)
                Test = 1;
            end
        end
        NbIts = NbIts + 1;
        if NbIts>40000
            disp('Can''t initialize IBM')
            die
        end
    end

    Names = Parameters.Names.Estimated;
    for j = 1:length(Names)
        Parameters.(Names{j}).Estimated = 0;
    end
    switch IndModel
        case 1
            Parameters.SigmaRW.Estimated = 1;
            Parameters.betainit.Estimated = 1;
            Parameters.Problem = 'Marc1diff';
        case 2
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.beta11init.Estimated = 1;
            Parameters.beta22init.Estimated = 1;
            Parameters.adultsmult.Estimated = 1;
            Parameters.kidsmult.Estimated = 1;
            Parameters.adultsadd.Estimated = 1;
            Parameters.kidsadd.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
            Parameters.Problem = 'Marc2diff';
            
        case 3
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.beta11init.Estimated = 1;
            Parameters.beta22init.Estimated = 1;
            Parameters.adultsmult.Estimated = 1;
            Parameters.kidsmult.Estimated = 1;
            Parameters.adultsadultsmult.Estimated = 1;
            Parameters.adultsadd.Estimated = 1;
            Parameters.kidsadd.Estimated = 1;
            Parameters.adultsadultsadd.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
            Parameters.Problem = 'Marc2diffb';

        case 4
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.SigmaRW22.Estimated = 1;
            Parameters.beta11init.Estimated = 1;
            Parameters.beta22init.Estimated = 1;
            Parameters.beta12init.Estimated = 1;
            Parameters.adultsmult.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
            Parameters.Problem = 'Marc2diff';

        case 5
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.SigmaRW22.Estimated = 1;
            Parameters.beta11init.Estimated = 1;
            Parameters.beta22init.Estimated = 1;
            Parameters.adultsmult.Estimated = 1;
            Parameters.kidsmult.Estimated = 1;
            Parameters.adultsadd.Estimated = 1;
            Parameters.kidsadd.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
            Parameters.Problem = 'Marc2diffb';

        case 6
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.SigmaRW22.Estimated = 1;
            Parameters.SigmaRW12.Estimated = 1;
            Parameters.beta11init.Estimated = 1;
            Parameters.beta22init.Estimated = 1;
            Parameters.beta12init.Estimated = 1;
            Parameters.adultsmult.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
            Parameters.Problem = 'Marc3diff';
         case 7
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.SigmaRW22.Estimated = 1;
            Parameters.SigmaRW12.Estimated = 1;
            Parameters.SigmaRW21.Estimated = 1;
            Parameters.beta11init.Estimated = 1;
            Parameters.beta22init.Estimated = 1;
            Parameters.beta12init.Estimated = 1;
            Parameters.beta21init.Estimated = 1;
            Parameters.E1InitProp.Estimated = 1;
            Parameters.I1InitProp.Estimated = 1;
            Parameters.R1InitProp.Estimated = 1;
            Parameters.E2InitProp.Estimated = 1;
            Parameters.I2InitProp.Estimated = 1;
            Parameters.R2InitProp.Estimated = 1;
            Parameters.Problem = 'Marc4diff';
    end
    % Parameters.EInitPropNoise.Estimated = 1;
    % Parameters.IInitPropNoise.Estimated = 1;
    % Parameters.RInitPropNoise.Estimated = 1;
    Parameters.km1.Estimated = 1;
    Parameters.gammam1.Estimated = 1;
    if strcmp(Parameters.DiffusionType,'IBM')
        Parameters.betaderinit.Estimated = 1;
    end
    % if strcmp(ObsType,'Estimated')
    %     Parameters.SigmaObs.Estimated = 1;
    % end
    Parameters = DefineEstimatedParametersIndexes(Parameters);
    Parameters = DefineTransfFunctions(Parameters);
    Parameters = DefinePriors(Parameters);
    Parameters = UpdateParsNoTransfToTransf(Parameters);

    Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);



    Parameters.km1.Estimated = 1;
    Parameters.gammam1.Estimated = 1;
    if strcmp(ObsType,'Estimated')
        Parameters.SigmaObs.Estimated = 1 ;
    elseif strcmp(ObsType,'CoeffEst')
        Parameters.MultCoeff.Estimated = 1 ;
    end
    if strcmp(Parameters.DiffusionType,'IBM')
        Parameters.betaderinit.Estimated = 0;
        %     Parameters.betaderinitNoise.Estimated = 1;
    end
    Parameters = DefineEstimatedParametersIndexes(Parameters);
    Parameters = DefineTransfFunctions(Parameters);
    Parameters = DefinePriors(Parameters);
    Parameters = UpdateParsNoTransfToTransf(Parameters);

    Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);
    KalHess = Parameters.KalHess;


    % Parameters.InitialCov = 0;
    % Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);


    Test = 0;

    NbIts = 0;

    while not(Test)

        Parameters = KalmOpt(Parameters,Data,SEIRModel,1500);

        try
            KalHess = Parameters.KalHess;
            Test = Parameters.KalmMaxReached;
        end

        NbIts = NbIts + 1;
        disp(NbIts)
        if NbIts>50
            return
        end
    end
    
    
    Temp = struct();
    Temp.Parameters = Parameters;
    save([SavePath '/Temp_' NameToSave],'Temp')
end

% SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan/Temp'

load([SavePath '/Temp_' NameToSave])
Parameters = Temp.Parameters;
KalHess = Parameters.KalHess;



Cov = (-KalHess)^-1;

Parameters.Correction = 0;
Parameters = KalmOpt(Parameters,Data,SEIRModel,1500);
Parameters.Correction = 1;

%% SMC Optimization

disp('Phase 2')

% Parameters.AdMet = 1;
% Parameters.AdMetBeta = 0.05;

Parameters.NoPaths = 1;
switch IndModel
    case 1
        Parameters.PathsToKeep = [1:7]';
        Parameters.NbParticules = 1000;
    case {2, 3, 4, 5}
        Parameters.PathsToKeep = [1:12]';
        Parameters.NbParticules = 3000;
    case 6 
        Parameters.PathsToKeep = [1:13]';
        Parameters.NbParticules = 6000;
    case 7
        Parameters.PathsToKeep = [1:14]';
        Parameters.NbParticules = 6000;
end


Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;


if strcmp(Parameters.DiffusionType,'IBM')
    Parameters.ComputationTStep = 1/3;
    Parameters.NbParticules = 10000;
else
%     Parameters.NbParticules = 1000;
    Parameters.ComputationTStep = 1/3;
end
    


Data.Instants = [0:size(Data.Observations,2)-1]*7/Parameters.ComputationTStep;
if IndModel>=2
    Data.ObservedVariables = diag([9 10])*ones(2,length(Data.Instants));
else
    Data.ObservedVariables = 5*ones(1,length(Data.Instants));
end
Data.NbComputingSteps = [0 diff(Data.Instants)];




% Running SMC from there
Cov = (-KalHess)^-1;
Parameters.G = Cov^-1;
Parameters.NoPaths = 0;
Parameters.ModelType='SMC';
Parameters.AdaptC = 0.999;
Parameters.NbVariables = 7;
Parameters.aim = 0.23;
Parameters.Epsil = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.PathsToKeep 
paths = [];
thetas = [];
names = Parameters.Names.Estimated;
for i = 1:1000
    i
    ResTemp = EstimationSMCsmoothGen(Data,SEIRModel,Parameters);
    RandInd = ceil(rand(1,1)*Parameters.NbParticules);
    paths(i,:,:) = squeeze(ResTemp.CompletePaths(RandInd,:,:));
    for j = 1:length(names)
        thetas(Parameters.(names{j}).Index,i) = Parameters.(names{j}).Value;
    end
end
Res.Paths = paths;
Res.Thetas = thetas;
save([SavePath '/TempSmoothed_' NameToSave],'Res')
clear Res




Cov = (-KalHess)^-1;
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.ModelType='SMC';
Parameters.AdaptC = 0.999;
Parameters.NbVariables = 7;
Parameters.aim = 0.23;
Parameters.Epsil = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.MCMCType
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
Parameters.ModelType='Kalman';
Parameters.AdaptC = 0.999;
Parameters.AdMet = 0;
Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,20000);
Cov = cov(Res.TransfThetas');
Parameters.G = Cov^-1;
Parameters.ModelType='Kalman';
Parameters.AdaptC = 0.999;
Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,20000);



save([SavePath '/Temp0_' NameToSave],'Res')

load([SavePath '/Temp0_' NameToSave])
Parameters = Res.Parameters;
% Parameters.NbParticules = 3000;

Parameters.ModelType = 'SMC';
dim = length(Parameters.Names.Estimated);
Cov = cov(Res.TransfThetas');
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
% if strcmp(Parameters.PMCMC,'Gibbs')
%     Parameters.NoPaths = 0;
% else
%     Parameters.NoPaths = 1;
% end
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
Parameters.AdaptC = 0.99;
Parameters.AdMet = 0;
Parameters.AdMetBeta = 0.05;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Res2 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbItersPrep);

save([SavePath '/Temp1_' NameToSave],'Res2')

load([SavePath '/Temp1_' NameToSave]) 

TempRes = Res2;

Parameters = TempRes.Parameters;
try
    Cov = cov(TempRes.TransfThetas');
end
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
Parameters.AdMet = 0;
Parameters.AdMetBeta = 0.05;
TempPar = TempRes.TempPar;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Res2 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbItersPrep);

% SavePath = 'S:\Results\';
% save([Name '_NoPaths.mat'],'Res2')

save([SavePath '/Temp2_' NameToSave],'Res2')

load([SavePath '/Temp2_' NameToSave]) 



TempRes = Res2;




Cov = cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
TempPar = TempRes.TempPar;
Parameters.NoPaths = 0;
Parameters.AdMet = 0;
Parameters.AdMetBeta = 0.05;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);

% if strcmp(Parameters.DiffusionType,'IBM')
%     Parameters.NoPaths = 0;
% else
%     Parameters.NoPaths = 0;
% end
if IndModel>1
    Parameters.PathsToKeep = [1:12]';
else
    Parameters.PathsToKeep = [1:7]';
end
Parameters.SaveSpace = 1;
Res3 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);


save([SavePath '/' NameToSave],'Res3')


