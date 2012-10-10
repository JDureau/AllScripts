% Google Flu Trends SEIRS


%% Load Data
DataPath = 'H:\PhD Work\Data\GoogleFluTrends';

fid = fopen([DataPath '\France.txt'],'r');
A = textscan(fid,'%s');
fclose(fid);
A = A{:};

regions = {'Alsace','Aquitaine','Auvergne','Burgundy','Brittany','Centre','Champagne-Ardenne','Franche-Comte','Ile-de-France','Languedoc-Roussillon','Lorraine','Midi-Pyrenees','Nord-Pas-de-Calais','Normandy - Lower','Normandy - Upper','Pays de la Loire','Picardie','Poitou-Charentes','Provence-Alpes-Cote d Azur','Rhône-Alpes'};
Data.Dates = {};
Data.NewCases = {};
for i = 2:length(A)
    Data.Dates{i-1} = A{i}(1:10);
    inds = [regexp(A{i}, ',') length(A{i})+1];
    for j = 1:length(inds)-1
        if inds(j)+1<inds(j+1)
            Data.NewCases{i-1}(j) = str2num(A{i}(inds(j)+1:inds(j+1)-1));
        else
            Data.NewCases{i-1}(j) = 0;
        end
    end
end

ToPlot = {'Alsace','Lorraine','Ile-de-France'};
temp = [];
for j = 1:length(ToPlot)
    ind = find(strcmp(regions,ToPlot{j}));
    for i = 1:length(Data.NewCases)
        temp(j,i) = Data.NewCases{i}(ind);
    end
end
plot(temp')
legend(ToPlot)

%% Inference
cd('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Filtering')
addpath('H:\PhD Work\Matlab Scripts\General Tools')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Parameter Estimation')
addpath('H:\PhD Work\Matlab Scripts\Toolboxes\Resampling\pf_resampling')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Joint Sampling')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\MIF')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Model Selection')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Optimization Approach')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Generic Parameter Estimation')
addpath('H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Generic Parameter Estimation\Models')
addpath('H:\PhD Work\Matlab Scripts\Toolboxes')

Parameters = struct();

Parameters.TotalPopulation = 12*10^6;
Parameters.k.Value = 1/3.5;
Parameters.k.Min = 1/4;
Parameters.k.Max = 1/3;
Parameters.k.Estimated = 1;
Parameters.k.TransfType = 'Log';
Parameters.gamma.Value = 1/7.5;
Parameters.gamma.Min = 1/8;
Parameters.gamma.Max = 1/7;
Parameters.gamma.Estimated = 1;
Parameters.gamma.TransfType = 'Log';
Parameters.betainit.Value = 0.25;
Parameters.betainit.Min = 0;
Parameters.betainit.Max = 10^14;
Parameters.betainit.Estimated = 1;
Parameters.betainit.TransfType = 'Log';
Parameters.EInitProp.Value = 0.000001;
Parameters.EInitProp.Min = 0;
Parameters.EInitProp.Max = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.EInitProp.TransfType = 'Log';
Parameters.IInitProp.Value = 0.000001;
Parameters.IInitProp.Min = 0;
Parameters.IInitProp.Max = 1;
Parameters.IInitProp.Estimated = 0;
Parameters.IInitProp.TransfType = 'Log';
Parameters.RInitProp.Value = 0.01;
Parameters.RInitProp.Min = 0;
Parameters.RInitProp.Max = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.RInitProp.TransfType = 'Log';
Parameters.SigmaRW.Value = exp(-0.6);
Parameters.SigmaRW.Min = 0;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Log';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)

Parameters.NbVariables = 6;
Parameters.SigmaObs = 0.05;
Parameters.DiffusionType = 'Add';
Parameters.ObservationLength = 7*50;
Parameters.ComputationTStep = 1;


% t0 = 28 jun 2009.
Data.Observations = zeros(5,50);
ind = find(strcmp(regions,'Alsace'));
for i = 1:50
    Data.Observations(5,i) = Data.NewCases{i+300}(ind);
end
Data.Instants = [1:50]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];

SEIRModel = struct();
SEIRModel.EKF_projection = @SEIR_EKF_projection2;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);







temp = zeros(1,6);
temp(1,5) = 1;
SEIRModel.ObservationJacobian = {};
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.Instants)
    SEIRModel.ObservationJacobian{i} = temp;
    SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs*Data.Observations(5,i))^2;
end

% SMC Optimization
Parameters.NoPaths = 1;
Parameters.NbParticules = 200;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    ParametersC.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)



Names = Parameters.Names.Estimated;
Cov = zeros(length(Names),length(Names));
for i = 1:length(Names)
    Cov(i,i) = (Parameters.(Names{i}).TransfValue*0.1)^2;
end


Parameters.Epsil = 1;
Parameters.ComputeRWsamplCov = 0;
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.aim = 0.23;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
[Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
ResGoog_RW_Calib = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,10000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\ResGoog_RW_Calib.mat'],'ResGoog_RW_Calib')


TempPar = ResGoog_RW.TempPar;
Parameters.NoPaths = 0;
Parameters.PathsToKeep = 1:6;
TransfThetasSamples = ResGoog_RW.TransfThetas;
Parameters = DefineScalingPars(TransfThetasSamples,Parameters);
ScaledTransfThetasSamples = Scale(TransfThetasSamples,Parameters);
Parameters.GMeth = 'GMM';
Parameters.DensityModel = FindBestDensityModel(ScaledTransfThetasSamples',Parameters,1);
[Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
ResGoog_RW = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,2000);

SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
save([SavePath '\ResGoog_RW.mat'],'ResGoog_RW')


Cov = zeros(3,3)
Cov(1,1) = 0.2^2;
Cov(2,2) = 10^2;
Cov(3,3) = 0.25;


Parameters.NbParticules = 500;
Parameters.PathsToKeep = [1:6]';
Parameters.MCMCType = 'Add';
ResultSMC = EstimationSMCsmoothGen(Data, SEIRModel, Parameters)

subplot(2,1,1)
plot(Data.Instants,ResultSMC.PosteriorMeansRecord(3,:))
hold on
plot(Data.Instants,Data.Observations(1,:),'g')
hold off
set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
subplot(2,1,2)
plot(squeeze(mean(exp(ResultSMC.CompletePaths(:,4,:)))))
hold on
plot(cumsum(Data.NbComputingSteps),0*ones(size(Data.NbComputingSteps)),'.k')
set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
title(ResultSMC.LogLik)
pause(0.01)


Res = ResGoog_RW;
% plot ResRun
figure(1)
for i = 1:5
    subplot(6,1,i)
    plot(mean(squeeze(Res.Paths(:,i,:))))
    hold on
%     plot(quantile(squeeze(Res.Paths(:,i,:)),0.025),'r')
%     plot(quantile(squeeze(Res.Paths(:,i,:)),0.975),'r')
    hold off
end
subplot(6,1,6)
plot(mean(squeeze(exp(Res.Paths(:,6,:)))))
hold on
plot(quantile(squeeze(exp(Res.Paths(:,6,:))),0.025),'r')
plot(quantile(squeeze(exp(Res.Paths(:,6,:))),0.975),'r')
hold off    
        
figure(2)
Names = Parameters.Names.Estimated;
k = length(Names);
for i = 1:k
    subplot(ceil(sqrt(k)),ceil(sqrt(k)),i)
    plot(Res.TransfThetas(i,:))   
%     plot(Res.Thetas(i,:))
    title(Names{i})
end
    
        
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).Estimated = 0;
end
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
        
        
Parameters.NbParticules = 10000;
LogLiks = [];
for i = 1:50
    i
    ResultSMC = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);
    LogLiks(i) = ResultSMC.LogLik;
end
hist(LogLiks)

Temp = ResultSMC
subplot(2,1,1)
plot(Data.Instants,Temp.PosteriorMeansRecord(3,:))
hold on
plot(Data.Instants,Data.Observations(5,:),'g')
hold off
set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
subplot(2,1,2)
plot(Data.Instants,exp(Temp.PosteriorMeansRecord(6,:)))
set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
title(Temp.LogLik)
pause(0.01)
        