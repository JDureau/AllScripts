%% Main Spline


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
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\SEIR')


SavePath = 'S:\Results\';
Name = [SavePath 'SEIR_simData.mat'];
load(Name)

Betas = DataGen.RealBetaTraj;
xbetas = 1:length(Betas);

% Spline = SplineFitBeta(Betas,25);
% Fit = spline(Spline.tis,Spline.yis,xbetas);
% plot(Betas,'g')
% hold on
% plot(Fit,'k')
% hold off

% Spline = SplineFitPrev(DataGen,DataGen.Parameters,100);
% Fit = spline(Spline.tis,Spline.yis,xbetas);
% plot(Betas,'g')
% hold on
% plot(Fit,'k')
% hold off

Parameters = Data.Parameters;

SEIRModel = struct();
SEIRModel.EKF_projection = @SEIR_EKFSpline_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);


NoiseLevels = [0.05 0.1 0.15 0.2 0.25 0.3];
NbSplines = [10 100];
Splines = {};
SMCPathss = {};
for i = 1:length(NoiseLevels)
    for j = 1:length(NbSplines)
       Betas = DataGen.RealBetaTraj;
       Parameters = DataGen.Parameters;
       Parameters.ObsNoise = NoiseLevels(i);
       NoisyData = SEIR_CreateData(Betas,Parameters,DataGen,SEIRModel);
       
       
        SplineInit = SplineFitBeta(Data.RealBetaTraj,NbSplines(j));

        tsinit = SplineInit.tis;
        yisinit = SplineInit.yis;

        Parameters.NbSplines = NbSplines(j);        
        n = Parameters.NbSplines;
        Initialization(1:n) = tsinit;
        Initialization(n+1:2*n) = yisinit;
        Names = Parameters.Names.Estimated;
        for k = 1:length(Names)
            Initialization(2*n+k) = Parameters.(Names{k}).TransfValue ;
        %         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
        end
        SEIRModel.InitializeParameters = @SEIRInitialize;
        [x,fval,exitflag,output] = fminsearch(@(x) KalmanSplineToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxFunEvals',2000,'MaxIter',2000,'TolX',1e-7,'TolFun',1e-7));
        Names = Parameters.Names.Estimated;
        tis = x(1:n) ;
        yis = x(n+1:2*n) ;
        for k = 1:length(Names)
            Parameters.(Names{k}).TransfValue = (x(2*n+k));
        end
        Parameters = UpdateParsTransfToNoTransf(Parameters);
        Parameters.Finalx = x;
        TellParsValues(Parameters)

        Spline.tis = tis;
        Spline.yis = yis;
        Spline.BetaFit = spline(Spline.tis,Spline.yis,1:sum(Data.NbComputingSteps));
        Splines{i}{j} = Spline;
        
        clf
        plot(Data.RealBetaTraj,'g')
        hold on
        plot(Spline.BetaFit)
        hold off
    end
    
   

    Names = Parameters.Names.Estimated;
    for i = 1:length(Names)
        Parameters.(Names{i}).Estimated = 0;
    end
    Parameters.betainit.Estimated = 1;
    Parameters.SigmaRW.Estimated = 1;
    Parameters = DefineEstimatedParametersIndexes(Parameters);
    Parameters = DefineTransfFunctions(Parameters);
    Parameters = DefinePriors(Parameters);
    Parameters = UpdateParsNoTransfToTransf(Parameters);
    Parameters.MCMCType = 'frfr';
    Parameters.NoPaths = 1;
    Parameters.NbParticules = 2000;
    Names = Parameters.Names.Estimated;
    Initialization = [];
    Names = Parameters.Names.Estimated;
    for i = 1:length(Names)
        Initialization(i) = Parameters.(Names{i}).TransfValue ;
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
    
    SMCPaths = [];
    NbIts = 100;
    Parameters.NoPaths = 0;
    Parameters.PathsToKeep = [1:7]';
    Parameters.NbParticules = 1000;
    for i = 1:NbIts
        i
        Temp = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);
        ind = ceil(rand(1,1)*Parameters.NbParticules);
        SMCPaths(end+1,:) = Temp.CompletePaths(ind,6,:);
    end

    SMCPathss{IndDataSet} = SMCPaths;
    
    
    SavePath = 'S:\Results\';
    Name = [SavePath 'SEIR_NoiseSplines.mat'];
    save(Name,'Splines')
    Name = [SavePath 'SEIR_NoiseSMC.mat'];
    save(Name,'SMCPathss')
    
end

Betas = DataGen.RealBetaTraj;
for i = 1:length(NoiseLevels)
    for j = 2:length(NbSplines)
        Spline = Splines{i}{j};
        Fit = spline(Spline.tis,Spline.yis,xbetas);
        plot(Betas,'g')
        hold on
        plot(Fit,'k')
        hold off
        pause()
    end
end    





