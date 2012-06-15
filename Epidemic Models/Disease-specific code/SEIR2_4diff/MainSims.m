% HPA data

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




%% Inference
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
tmp = load([SavePath 'ParametersSEIR_MarcAdd.mat']);
Parameters = tmp.Parameters;

TStep = Parameters.ComputationTStep;


years = 2003:2010;
regions = {'France','Alsace','Aquitaine','Auvergne','Burgundy','Brittany','Centre','Champagne-Ardenne','Franche-Comte','Ile-de-France','Languedoc-Roussillon','Lorraine','Midi-Pyrenees','Nord-Pas-de-Calais','Normandy - Lower','Normandy - Upper','Pays de la Loire','Picardie','Poitou-Charentes','Provence-Alpes-Cote d Azur','Rhône-Alpes'};

  

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






Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:7]';
Parameters.TotalPopulation = 100000;



test = 0;
while not(test)
    region = regions{ceil(rand(1)*length(regions))};
    year = years(ceil(rand(1)*length(years)));
    disp([ region ' - ' num2str(year)])
    try 
        Data = ExtractGoogleTSeries(DataFr,Parameters,region,year);
        tmp = find(not(Data.Observations(5,:)==0));
        if tmp(1)<5
            Data.Observations(5,1:tmp(1)-1) = Data.Observations(5,tmp(1));
            test = 1;
        end
    end
end
plot(Data.Observations(5,:))
temp = zeros(1,6);
temp(1,5) = 1;
Parameters.betainit.Value = 0.7;
Parameters.EInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Min = 0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Min = 0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);  

SEIRModel = struct();
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);
SEIRModel.ObservationJacobian = {};
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.Instants)
    SEIRModel.ObservationJacobian{i} = temp;
    SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs.Value*Data.Observations(5,i))^2;
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
% Parameters.betainitMean.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;

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
    Parameters.SigmaRW.Value = rand(1,1)*2;
%     Parameters.EInitPropNoise.Value = rand(1,1)*0.3;
%     Parameters.IInitPropNoise.Value = rand(1,1)*0.3;
%     Parameters.RInitPropNoise.Value = rand(1,1)*0.3;
    Parameters.InitialCov = rand(1,1)*0.4;
    Parameters = UpdateParsNoTransfToTransf(Parameters);
    try
        Temp = EstimationEKFGen(Data, SEIRModel, Parameters);
        if (Temp.LogLik>-300)
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
Parameters.SigmaRW.Estimated = 1;
Parameters.betainit.Estimated = 1;
% Parameters.EInitPropNoise.Estimated = 1;
% Parameters.IInitPropNoise.Estimated = 1;
% Parameters.RInitPropNoise.Estimated = 1;
Parameters.InitialCov = Parameters.InitialCov/2;
if strcmp(Parameters.DiffusionType,'IBM')
    Parameters.betaderinit.Estimated = 1;
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);



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
% Parameters.EInitPropNoise.Estimated = 1;
% Parameters.IInitPropNoise.Estimated = 1;
% Parameters.RInitPropNoise.Estimated = 1;
% Parameters.betainitNoise.Estimated = 1;
Parameters.InitialCov = 0.15;
if strcmp(Parameters.DiffusionType,'IBM')
    Parameters.betaderinit.Estimated = 1;
    Parameters.InitialCov = Parameters.InitialCov/2;
    %     Parameters.betaderinitNoise.Estimated = 1;
end
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);


Parameters.InitialCov = 0.1;
Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);

Parameters.InitialCov = 0.05;
Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);

Parameters.InitialCov = 0;
Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);



   

Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:7]';
ResultSMC = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);

betas = mean(squeeze(exp(ResultSMC.CompletePaths(:,6,:))));

Parameters.ObsNoise = 0.1;
DataGen = SEIR_CreateData(betas,Parameters,Data,SEIRModel);

plot(DataGen.Observations(5,:))

SavePath = 'S:\Results\';
Name = [SavePath 'SEIR_simData.mat'];
DataGen.Parameters = Parameters;
save(Name,'DataGen')

Name = [SavePath 'SimData_Add.mat'];
Difftype = 'Add';
ObsType = 'Fixed';
FullSEIRinference(Data,Difftype,ObsType,Name)


%% Many Sims

NbIts = 100;
DatasGens = {};
for IndIt = 21:NbIts
    

    Parameters.NoPaths = 0;
    Parameters.PathsToKeep = [1:7]';
    Parameters.TotalPopulation = 100000;



    test = 0;
    while not(test)
        region = regions{ceil(rand(1)*length(regions))};
        year = years(ceil(rand(1)*length(years)));
        disp([ region ' - ' num2str(year)])
        try 
            Data = ExtractGoogleTSeries(DataFr,Parameters,region,year);
            tmp = find(not(Data.Observations(5,:)==0));
            if tmp(1)<5
                Data.Observations(5,1:tmp(1)-1) = Data.Observations(5,tmp(1));
                test = 1;
            end
        end
    end
    plot(Data.Observations(5,:))
    temp = zeros(1,6);
    temp(1,5) = 1;
    Parameters.betainit.Value = 0.7;
    Parameters.EInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
    Parameters.EInitProp.Min = 0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
    Parameters.EInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
    Parameters.IInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
    Parameters.IInitProp.Min = 0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
    Parameters.IInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
    Parameters = DefineEstimatedParametersIndexes(Parameters);
    Parameters = DefineTransfFunctions(Parameters);
    Parameters = DefinePriors(Parameters);
    Parameters = UpdateParsNoTransfToTransf(Parameters);  

    try
        SEIRModel = struct();
        SEIRModel.EKF_projection = @SEIR_EKF_projection;
        SEIRModel.InitializeParameters = @SEIRInitialize;
        SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
        SEIRModel.SMC_projection = @SEIR_SMC_projection;
        Parameters = SEIRModel.InitializeParameters(Parameters);
        SEIRModel.ObservationJacobian = {};
        SEIRModel.ObservationMeasurementNoise = {};
        for i = 1:length(Data.Instants)
            SEIRModel.ObservationJacobian{i} = temp;
            SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs.Value*Data.Observations(5,i))^2;
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
        % Parameters.betainitMean.Estimated = 1;
        Parameters.SigmaRW.Estimated = 1;

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
            Parameters.SigmaRW.Value = rand(1,1)*2;
        %     Parameters.EInitPropNoise.Value = rand(1,1)*0.3;
        %     Parameters.IInitPropNoise.Value = rand(1,1)*0.3;
        %     Parameters.RInitPropNoise.Value = rand(1,1)*0.3;
            Parameters.InitialCov = rand(1,1)*0.4;
            Parameters = UpdateParsNoTransfToTransf(Parameters);
            try
                Temp = EstimationEKFGen(Data, SEIRModel, Parameters);
                if (Temp.LogLik>-300)
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
        Parameters.SigmaRW.Estimated = 1;
        Parameters.betainit.Estimated = 1;
        % Parameters.EInitPropNoise.Estimated = 1;
        % Parameters.IInitPropNoise.Estimated = 1;
        % Parameters.RInitPropNoise.Estimated = 1;
        Parameters.InitialCov = Parameters.InitialCov/2;
        if strcmp(Parameters.DiffusionType,'IBM')
            Parameters.betaderinit.Estimated = 1;
        end
        Parameters = DefineEstimatedParametersIndexes(Parameters);
        Parameters = DefineTransfFunctions(Parameters);
        Parameters = DefinePriors(Parameters);
        Parameters = UpdateParsNoTransfToTransf(Parameters);


        Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);



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
        % Parameters.EInitPropNoise.Estimated = 1;
        % Parameters.IInitPropNoise.Estimated = 1;
        % Parameters.RInitPropNoise.Estimated = 1;
        % Parameters.betainitNoise.Estimated = 1;
        Parameters.InitialCov = 0.15;
        if strcmp(Parameters.DiffusionType,'IBM')
            Parameters.betaderinit.Estimated = 1;
            Parameters.InitialCov = Parameters.InitialCov/2;
            %     Parameters.betaderinitNoise.Estimated = 1;
        end
        Parameters.SigmaRW.Estimated = 1;
        Parameters = DefineEstimatedParametersIndexes(Parameters);
        Parameters = DefineTransfFunctions(Parameters);
        Parameters = DefinePriors(Parameters);
        Parameters = UpdateParsNoTransfToTransf(Parameters);

        Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);


        Parameters.InitialCov = 0.1;
        Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);

        Parameters.InitialCov = 0.05;
        Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);

        Parameters.InitialCov = 0;
        Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);





        Parameters.NoPaths = 0;
        Parameters.PathsToKeep = [1:7]';
        ResultSMC = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);

        betas = mean(squeeze(exp(ResultSMC.CompletePaths(:,6,:))));

        Parameters.ObsNoise = 0.1;
        DataGen = SEIR_CreateData(betas,Parameters,Data,SEIRModel);

        plot(DataGen.Observations(5,:))

        SavePath = 'S:\Results\';
        Name = [SavePath 'SEIR_simsDatas.mat'];
        DataGen.Parameters = Parameters;

        DatasGens{end+1} = DataGen;

        save(Name,'DatasGens')
    end
end

DataGensTot = {};
cpt = 0;
exts = {'','2','3','4'};
for i = 1:length(exts)
    cpt = length(DataGensTot);
    SavePath = 'S:\Results\';
    Name = [SavePath 'SEIR_simsDatas' exts{i} '.mat'];
    load(Name)
    for j = 1:length(DatasGens)
        DataGensTot{end+1} = DatasGens{j};
        plot(DatasGens{j}.RealBetaTraj)
        title(length(DataGensTot))
        pause()
        
    end
end

SavePath = 'S:\Results\';
Name = [SavePath 'SEIR_simsDatas.mat'];
DatasGens =  DataGensTot ;
save(Name,'DatasGens')



load(Name)
for i = 1:length(DatasGens)
    clf
    DatasSim  = SEIR_CreateData(DatasGens{i}.RealBetaTraj,DatasGens{i}.Parameters,Data,SEIRModel);
    subplot(2,1,1)
    plot(DatasGens{i}.Observations(5,:))
    hold on
    plot(DatasSim.Observations(5,:),'g')
    hold off
    subplot(2,1,2)
    plot(DatasGens{i}.RealBetaTraj)
    pause()
end

