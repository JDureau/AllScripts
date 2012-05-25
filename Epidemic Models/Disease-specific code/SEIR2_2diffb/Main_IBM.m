% tests IBM

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

%%
DataPath = 'H:\My Documents\PhD Work\Data\HPA';

A = load([DataPath '\andre_estimates_31_01.txt']);

Data.Dates = {};
Data.NewCases = {};
for i = 1:size(A,1)
    Data.Dates{i} = i;
    for j = 1:7
    Data.NewCases{i}{j} = A(i,j)*10;
    end
end

InitialDate = struct();
InitialDate.Month = 6;
InitialDate.Day = 1;
InitialDate.Year = 2009;
Data.Dates = ApproxWeeklyDates(InitialDate,35);

plot(A)
legend('0-4','5-14','15-24','25-44','45-64','65+')

PopWeigths = [667600,2461800,5904100,6862500,14417400,12847800,7929300];
PopProps = PopWeigths/sum(PopWeigths)*100;
Weigthed = sum((A(:,1:7)*diag(PopWeigths.^-1)*diag(PopProps/100)*100000)');

plot(A*diag(PopWeigths.^-1)*100000)
hold on
plot(Weigthed,'k','LineWidth',2)
hold off
legend('<1','1-4','5-14','15-24','25-44','45-64','65+')


Parameters = struct();

Parameters.Problem = 'MarcFlu';
Parameters.NbVariables = 7;
Parameters.SigmaObs = 0.1;
Parameters.DiffusionType = 'Add';
Parameters.ObservationLength = 7*35;
Parameters.ComputationTStep = 0.1;
Parameters.TotalPopulation = 100000;

Data.Observations = zeros(6,35);
Data.Observations(5,:) = Weigthed*10;
Data.Instants = [1:35]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];

disp('Be careful with total population')

Parameters.km1.Value = 2;
Parameters.km1.Min = -10^14;
Parameters.km1.Max = 10^14;
Parameters.km1.MinLim = 0.1;
Parameters.km1.MaxLim = 2.1;
Parameters.km1.Estimated = 1;
Parameters.km1.TransfType = 'Logit';
Parameters.km1.PlotInv = 1;
Parameters.gammam1.Value = 2;
Parameters.gammam1.Min = -10^14;
Parameters.gammam1.Max = 10^14;
Parameters.gammam1.MinLim = 0.4;
Parameters.gammam1.MaxLim = 2.1;
Parameters.gammam1.Estimated = 1;
Parameters.gammam1.TransfType = 'Logit';
Parameters.gammam1.PlotInv = 1;
Parameters.betainit.Value = 1;
Parameters.betainit.Min = -10^14;
Parameters.betainit.Max = 10^14;
Parameters.betainit.Estimated = 1;
Parameters.betainit.TransfType = 'Log';
Parameters.betainit.Init = 1;
Parameters.betaderinit.Value = 0;
Parameters.betaderinit.Min = -10^14;
Parameters.betaderinit.Max = 10^14;
Parameters.betaderinit.MinLim = -10;
Parameters.betaderinit.MaxLim = 10;
Parameters.betaderinit.Estimated = 1;
Parameters.betaderinit.TransfType = 'Logit';
Parameters.betaderinit.Init = 1;
Parameters.EInitProp.Value = max(eps,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Min = -10^14;%0.2*max(eps,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Max = 10^14;%5*max(eps,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.MinLim = 0;
Parameters.EInitProp.MaxLim = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.EInitProp.TransfType = 'Logit';
Parameters.EInitProp.Init = 1;
Parameters.IInitProp.Value = max(eps,0.2*Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Min = -10^14;%0.1*max(eps,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Max = 10^14;%5*max(eps,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.MinLim = 0;
Parameters.IInitProp.MaxLim = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.IInitProp.TransfType = 'Logit';
Parameters.IInitProp.Init = 1;
Parameters.RInitProp.Value = 0.2;
Parameters.RInitProp.Min = 0;
Parameters.RInitProp.Max = 0.30;
Parameters.RInitProp.MinLim = 0;
Parameters.RInitProp.MaxLim = 0.6;
Parameters.RInitProp.Estimated = 1;
Parameters.RInitProp.TransfType = 'Logit';
Parameters.RInitProp.Init = 1;
Parameters.SigmaRW.Value = exp(-0.6);
Parameters.SigmaRW.Min = -10^140;
Parameters.SigmaRW.Max = 10^140;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Log';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)


    

SEIRModel = struct();
SEIRModel.UKF_projection = @SEIR_UKF_projection;
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);

ParametersInit = Parameters;

% EKF Optimization
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];



temp = zeros(1,6);
temp(1,5) = 1;
SEIRModel.ObservationJacobian = {};
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.Instants)
    SEIRModel.ObservationJacobian{i} = temp;
    SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs*Data.Observations(5,i))^2;
end


Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


DataPath = 'H:\My Documents\PhD Work\Data\GoogleFluTrends';
fid = fopen([DataPath '\France.txt'],'r');
A = textscan(fid,'%s');
fclose(fid);
A = A{:};

years = 2003:2010;
regions = {'France','Alsace','Aquitaine','Auvergne','Burgundy','Brittany','Centre','Champagne-Ardenne','Franche-Comte','Ile-de-France','Languedoc-Roussillon','Lorraine','Midi-Pyrenees','Nord-Pas-de-Calais','Normandy - Lower','Normandy - Upper','Pays de la Loire','Picardie','Poitou-Charentes','Provence-Alpes-Cote d Azur','Rhône-Alpes'};
DataFr.Dates = {};
DataFr.NewCases = {};
for i = 2:length(A)
    DataFr.Dates{i-1} = A{i}(1:10);
    inds = [regexp(A{i}, ',') length(A{i})+1];
    for j = 1:length(inds)-1
        if inds(j)+1<inds(j+1)
            DataFr.NewCases{i-1}(j) = str2num(A{i}(inds(j)+1:inds(j+1)-1));
        else
            DataFr.NewCases{i-1}(j) = 0;
        end
    end
end
DataFr.regions = regions;
ParametersInit = Parameters;

%%


region = regions{ceil(rand(1)*length(regions))};
year = years(ceil(rand(1)*length(years)));
disp([ region ' - ' num2str(year)])

DataTemp = ExtractGoogleTSeries(DataFr,Parameters,region,year);
plot(DataTemp.Observations(5,:))


Parameters.DiffusionType = 'Add';
Parameters.DiffusionType = 'IBM';


Parameters.TotalPopulation = 1;
Parameters.EInitProp.Value = max(eps,1.5*DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.MinLim = 0.001*max(eps,DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.MaxLim = 2*max(eps,DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Value = max(eps,1.5*DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.MinLim = 0.001*max(eps,DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.MaxLim = 2*max(eps,DataTemp.Observations(5,1))/Parameters.TotalPopulation;
Parameters.betainit.Value = 0.6;
Parameters = UpdateParsNoTransfToTransf(Parameters);
temp = zeros(1,7);
temp(1,5) = 1;
SEIRModel.ObservationJacobian = {};
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(DataTemp.Instants)
    SEIRModel.ObservationJacobian{i} = temp;
    SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs*DataTemp.Observations(5,i))^2;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
Names = Parameters.Names.Estimated;
Parameters = SEIRModel.InitializeParameters(Parameters);
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
fval = Inf;
Parameters.InitialCovFact = 0.00002;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,DataTemp,SEIRModel,Parameters),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-5));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)


Parameters.InitialCovFact = 0.00000000;
Temp = EstimationEKFGen(DataTemp, SEIRModel, Parameters);


Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:6]';
Parameters.NbParticules = 4000;
Temp = EstimationSMCsmoothGen(DataTemp, SEIRModel, Parameters);

PlotResGoogle(Temp,8)


traj = mean(squeeze(Temp.CompletePaths(:,6,:)));
plot(diff(traj)/Parameters.ComputationTStep)

betas = 0.7*ones(1,3570);%+0.0001*(1:3570);
t = 1:3570;
betas = 0.7*ones(1,3570) + 100*normpdf(t,1500,500);
plot(betas)
DataGen = SEIR_CreateData(betas,Parameters,DataTemp, SEIRModel);

% 
% Parameters.NbParticules = 1000;
Temp = EstimationEKFGen(DataGen, SEIRModel, Parameters);
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
end
fval = Inf;
Parameters.InitialCovFact = 0.3;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,DataGen,SEIRModel,Parameters),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-5));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)



Temp = EstimationSMCsmoothGen(DataGen, SEIRModel, Parameters);

PlotResGoogle(Temp,8)

ampls = 0:10;
DiffType = {'Add','IBM'};
    
for i = 1:length(ampls)
    for j = 1:2
        Parameters.DiffusionType = DiffType{j};
        
        betas = 0.7*ones(1,3570) + 100*normpdf(t,1500,500);
        plot(betas)
        DataGen = SEIR_CreateData(betas,Parameters,DataTemp, SEIRModel);

        Temp = EstimationSMCsmoothGen(DataGen, SEIRModel, Parameters);

        Ress{i,j} = Temp;
    end
end  
    
    
    
    


















