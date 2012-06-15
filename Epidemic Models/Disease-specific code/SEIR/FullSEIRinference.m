function [] = FullSEIRinference(Data,DiffType,ObsType,Name,IndModel)


% IndModel = 1 -> normal
% IndModel = 2 -> 2diff
% IndModel = 3 -> 2diffb
% IndModel = 4 -> 3diff
        
NbIters = 30000;
NbItersPrep = 300;



SavePath = '/users/ecologie/dureau/src/AllData/ResultsMarc/';

switch IndModel
    case 1
        tmp = load([SavePath '/ParametersSEIR.mat']);
    case 2
        tmp = load([SavePath '/ParametersSEIR2_1diff.mat']);
    case 3
        tmp = load([SavePath '/ParametersSEIR2_2diff.mat']);
    case 4
        tmp = load([SavePath '/ParametersSEIR2_2diffb.mat']);
    case 5
        tmp = load([SavePath '/ParametersSEIR2_3diff.mat']);
end



switch IndModel
    case 1
        NameToSave = ['MarcData_StructModel_normal.mat'];
    case 2
        NameToSave = ['MarcData_StructModel_1diff.mat'];
    case 3
        NameToSave = ['MarcData_StructModel_2diff.mat'];
    case 4
        NameToSave = ['MarcData_StructModel_2diffb.mat'];
    case 5
        NameToSave = ['MarcData_StructModel_3diff.mat'];
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

load([SavePath '/Temp_' NameToSave])

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
        SEIRModel.EKF_projection = @SEIR2_2diff_EKF_projection;
        SEIRModel.InitializeParameters = @SEIR2_2diff_Initialize;
        SEIRModel.SMC_projection = @SEIR2_2diff_SMC_projection;
        SEIRModel.LikFunction = 'mvnpdf(log(Variables(:,[9 10])),transpose(log(coeff*Data.Observations([9 10],IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),diag([log(Parameters.SigmaObs.Value^2+1) log(Parameters.SigmaObs.Value^2+1)]))';
    case 4
        SEIRModel.EKF_projection = @SEIR2_2diffb_EKF_projection;
        SEIRModel.InitializeParameters = @SEIR2_2diffb_Initialize;
        SEIRModel.SMC_projection = @SEIR2_2diffb_SMC_projection;
        SEIRModel.LikFunction = 'mvnpdf(log(Variables(:,[9 10])),transpose(log(coeff*Data.Observations([9 10],IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),diag([log(Parameters.SigmaObs.Value^2+1) log(Parameters.SigmaObs.Value^2+1)]))';
    case 5
        SEIRModel.EKF_projection = @SEIR2_3diff_EKF_projection;
        SEIRModel.InitializeParameters = @SEIR2_3diff_Initialize;
        SEIRModel.SMC_projection = @SEIR2_3diff_SMC_projection;
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
    %         Parameters.adultsadd.Estimated = 1;
    %         Parameters.kidsadd.Estimated = 1;
        case 3
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.SigmaRW22.Estimated = 1;
        case 4
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.SigmaRW22.Estimated = 1;

        case 5
            Parameters.SigmaRW11.Estimated = 1;
            Parameters.SigmaRW22.Estimated = 1;
            Parameters.SigmaRW12.Estimated = 1;
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
        if (Temp.LogLik>-300)
            Test = 1;
        end
    end
    NbIts = 0;
    while not(Test)
        Parameters = SampleParameters(Parameters);
        if IndModel>2
            Parameters.SigmaRW11.Value = rand(1,1)*2;
            Parameters.SigmaRW22.Value = rand(1,1)*2;
        elseif IndModel==2
            Parameters.SigmaRW11.Value = rand(1,1)*2;
            Parameters.adultsmult.Value = rand(1,1);
            Parameters.kidsmult.Value = rand(1,1);
            Parameters.R2InitProp.Value = 0.5;
            Parameters.R2InitProp.Value = 0.1;
        else
            Parameters.SigmaRW.Value = rand(1,1)*2;
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
        if NbIts>20000
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

        case 4
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
            Parameters.Problem = 'Marc2diff';

        case 5
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



    Names = Parameters.Names.Estimated;
    for j = 1:length(Names)
        Parameters.(Names{j}).Estimated = 0;
    end
    switch IndModel
        case 1
            Parameters.SigmaRW.Estimated = 1;
            Parameters.betainit.Estimated = 1;
            Parameters.EInitProp.Estimated = 1;
            Parameters.IInitProp.Estimated = 1;
            Parameters.RInitProp.Estimated = 1;
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

        case 4
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

        case 5
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
    end
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

    if IndModel == 1
        Temp = struct();
        Temp.Parameters = Parameters;
        save([SavePath '/Temp_' NameToSave],'Temp')
    end
end

% SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan/Temp'

load([SavePath '/Temp_' NameToSave])
Parameters = Temp.Parameters;
KalHess = Parameters.KalHess;



Cov = (-KalHess)^-1;

Parameters.Correction = 0;
% Parameters = KalmOpt(Parameters,Data,SEIRModel,1500);
Parameters.Correction = 1;


%% SMC Optimization

disp('Phase 2')

% Parameters.AdMet = 1;
% Parameters.AdMetBeta = 0.05;

Parameters.NoPaths = 1;
switch IndModel
    case 1
        Parameters.PathsToKeep = [1:7]';
    case {2, 3, 4}
        Parameters.PathsToKeep = [1:12]';
    case 5
        Parameters.PathsToKeep = [1:13]';
end

if IndModel >2
    Parameters.NbParticules = 4000;
else
    Parameters.NbParticules = 1000;
end
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;

% Initialization = [];
% Names = Parameters.Names.Estimated;
% for i = 1:length(Names)
%     Initialization(i) = Parameters.(Names{i}).TransfValue ;
% %         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
% end
% Parameters.Correction = 1;
% [x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizewithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',20,'TolX',1e-8,'TolFun',1e-8));
% Names = Parameters.Names.Estimated;
% for i = 1:length(Names)
%     Parameters.(Names{i}).TransfValue = (x(i));
% end
% Parameters = UpdateParsTransfToNoTransf(Parameters);
% TellParsValues(Parameters)  

% Temp = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);
% Temp.LogLik
% 
% % defining NbParts 
% nbparts = [  100 1000 2000 3000 4000  5000 ];
% % nbparts = [  5000 6000 7000 8000 10000];
% NbTests = 50;
% LogLikrecords = [];
% for i = 1:length(nbparts)
%     Parameters.NbParticules = nbparts(i);
%     for j = 1:NbTests
%         disp([i j])
%         Temp = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);
%         LogLikrecords(i,j) = Temp.LogLik;
%     end
% end
% clf
% hold on
% for i = 1:length(nbparts)
%     plot(nbparts(i),std(LogLikrecords(i,:)),'o')
% end
% hold off
% clf
% hold on
% for i = 1:length(nbparts)
%     plot(nbparts(i),mean(LogLikrecords(i,:)),'.')
%     plot(nbparts(i),mean(LogLikrecords(i,:))+std(LogLikrecords(i,:)),'.r')
%     plot(nbparts(i),mean(LogLikrecords(i,:))-std(LogLikrecords(i,:)),'.r')
% end
% hold off
% 
% 
% Parameters.NbParticules = 4000;
% 
% % defining TStep
% DataTests = Data;
% 
% res = [2 3 10].^-1;
% logliks = [];
% paths = [];
% for i = 1:length(res)
%     disp(i)
%     for j = 1:20
%         Parameters.ComputationTStep = res(i);
%         DataTests.Instants = [1:size(Data.Observations,2)]*7/Parameters.ComputationTStep;
%         DataTests.ObservedVariables = 5*ones(1,length(DataTests.Instants));
%         DataTests.NbComputingSteps = [0 diff(DataTests.Instants)];
%         Temp = EstimationSMCsmoothGen(DataTests, SEIRModel, Parameters);
%         logliks(i,j) = Temp.LogLik;
%         paths(i,j,:) = exp(Temp.PosteriorMeansRecord(6,:));
%     end
% end
% logliks'
% clf
% hold on
% for i = 1:length(res)
%     plot(res(i),mean(logliks(i,:)),'.')
%     plot(res(i),mean(logliks(i,:))+std(logliks(i,:)),'.r')
%     plot(res(i),mean(logliks(i,:))-std(logliks(i,:)),'.r')
% end
% hold off
% 
% cols = rand(length(res),3);
% clf
% hold on
% for i = 1:length(res)
%     plot(mean(squeeze(paths(i,:,:))),'col',cols(i,:),'LineWidth',2)
%     plot(quantile(squeeze(paths(i,:)),0.95),'col',cols(i,:))
%     plot(quantile(squeeze(paths(i,:)),0.05),'col',cols(i,:))
% end
% legend
% hold off


% if strcmp(Parameters.PMCMC,'Gibbs')
%     tmpCov = (-KalHess)^-1;
%     for i = 1:length(Names)
%         ind = Parameters.(Names{i}).Index;
%         Parameters.(Names{i}).GibbsSigma = sqrt(2.38^2*tmpCov(ind,ind));
%         Parameters.(Names{i}).GibbsSampler = @GibbsRWsampler;
%     end
% end


if strcmp(Parameters.DiffusionType,'IBM')
    Parameters.ComputationTStep = 1/3;
    Parameters.NbParticules = 10000;
else
    Parameters.NbParticules = 1000;
    Parameters.ComputationTStep = 1/3;
end
    
if IndModel >=2
    Parameters.NbParticules = 1000;
end

Data.Instants = [0:size(Data.Observations,2)-1]*7/Parameters.ComputationTStep;
if IndModel>=2
    Data.ObservedVariables = diag([9 10])*ones(2,length(Data.Instants));
else
    Data.ObservedVariables = 5*ones(1,length(Data.Instants));
end
Data.NbComputingSteps = [0 diff(Data.Instants)];




% Cov = (-KalHess)^-1;
% Parameters.G = Cov^-1;
% Parameters.NoPaths = 1;
% Parameters.ModelType='SMC';
% Parameters.AdaptC = 0.999;
% Parameters.NbVariables = 7;
% Parameters.aim = 0.23;
% Parameters.Epsil = 1;
% Parameters.MCMCType = 'Rand';
% Parameters.GMeth = 'cst given';
% Parameters.MCMCType
% TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
% Parameters.ModelType='Kalman';
% Parameters.AdaptC = 0.999;
% Parameters.AdMet = 0;
% Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,20000);
% Cov = cov(Res.TransfThetas');
% Parameters.G = Cov^-1;
% Parameters.ModelType='Kalman';
% Parameters.AdaptC = 0.999;
% Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,20000);
% 
% 
% 
% save([SavePath '/Temp0_' NameToSave],'Res')

load([SavePath '/Temp0_' NameToSave])
Parameters = Res.Parameters;
Parameters.NbParticules = 3000;

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
% Res2 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbItersPrep);

% save([SavePath '/Temp1_' NameToSave],'Res2')

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
% Res2 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbItersPrep);

% SavePath = 'S:\Results\';
% save([Name '_NoPaths.mat'],'Res2')

% save([SavePath '/Temp2_' NameToSave],'Res2')

load([SavePath '/Temp2_' NameToSave]) 



TempRes = Res2;


% Running SMC from there
Parameters = TempRes.Parameters;
Parameters.NoPaths = 0;
IndModel
switch IndModel
    case 1
        Parameters.PathsToKeep = [1:7]';
    case {2, 3, 4}
        'you'
        Parameters.PathsToKeep = [1:12]';
    case 5
        Parameters.PathsToKeep = [1:13]';
end
Parameters.PathsToKeep 
paths = [];
thetas = [];
names = Parameters.Names.Estimated;
for i = 1:10
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
Parameters.SaveSpace = 0;
Res3 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);


save([SavePath '/' NameToSave],'Res3')


