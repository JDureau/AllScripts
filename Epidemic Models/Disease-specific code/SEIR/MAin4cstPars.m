% Testing deterministic model with 4 constant beta's

%% Generating data

cd('/Users/dureaujoseph/AllScripts')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR2_2diff_diff'])
SavePath = '/Users/dureaujoseph/Documents/PhD_Data/ResultsMarc/';
load([SavePath '/ParametersSEIR2_2diff_diff.mat']);


load([SavePath '/Temp2_MarcData_StructModel_2diff_diff.mat'])
Data = Res.Data;
Parameters = Res.Parameters;
Parameters.SigmaRW11.Value = 0;
Parameters.SigmaRW22.Value = 0;

SEIRModel.EKF_projection = @SEIR2_2diff_diff_EKF_projection;
SEIRModel.InitializeParameters = @SEIR2_2diff_diff_Initialize;
SEIRModel.SMC_projection = @SEIR2_2diff_diff_SMC_projection;
SEIRModel.LikFunction = 'mvnpdf(log(Variables(:,[9 10])),transpose(log(coeff*Data.Observations([9 10],IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),diag([log(Parameters.SigmaObs.Value^2+1) log(Parameters.SigmaObs.Value^2+1)]))';


Parameters.SigmaRW11.Value = 0.00001;
Parameters.SigmaRW22.Value = 0.00001;
Parameters.beta11init.Value = 2;
Parameters.beta22init.Value = 1;
Parameters.beta12init.Value = 0.1;
Parameters.beta21init = Parameters.beta12init;
Parameters.beta21init.Value = 0.6;
Parameters.SigmaObs.Value = 100;
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.NoPaths = 0;

ResultSMC = EstimationSMCsmoothGen(Data, SEIRModel, Parameters)

clf
subplot(2,1,1)
plot(squeeze(ResultSMC.CompletePaths(1,9,Data.Instants+1)))
subplot(2,1,2)
plot(squeeze(ResultSMC.CompletePaths(1,10,Data.Instants+1)))

Data.Observations(9,:)  = squeeze(ResultSMC.CompletePaths(1,9,Data.Instants+1));
Data.Observations(10,:) = squeeze(ResultSMC.CompletePaths(1,10,Data.Instants+1));


%% Inference

Names = Parameters.Names.All;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.beta11init.Estimated = 1;
Parameters.beta22init.Estimated = 1;
Parameters.beta12init.Estimated = 1;
Parameters.beta21init.Estimated = 1;
Parameters.SigmaObs.Value = 0.1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);

KalHess = Parameters.KalHess;
Cov = (-KalHess)^-1;
Parameters.Correction = 0;
Parameters = KalmOpt(Parameters,Data,SEIRModel,1500);
Parameters.Correction = 1;



Parameters.NbParticules = 2;


Cov = (-KalHess)^-1;
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.AdaptC = 0.999;
Parameters.aim = 0.23;
Parameters.Epsil = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
Parameters.ModelType='SMC';
Parameters.AdaptC = 0.999;
Parameters.AdMet = 0;
Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,10000);
PlotTracePlots(Res)

for i = 1:10
    Cov = cov(Res.TransfThetas');
    Parameters.G = Cov^-1;
    Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,10000);
    save([SavePath '/Temp' num2str(i) '_Model4cstbetas' ],'Res')
end

Cov = cov(Res.TransfThetas');
Parameters.G = Cov^-1;
Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,50000);

save([SavePath '/Temp' num2str(11) '_Model4cstbetas' ],'Res')

Cov = cov(Res.TransfThetas');
Parameters.G = Cov^-1;
Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,100000);

save([SavePath '/Temp' num2str(12) '_Model4cstbetas' ],'Res')

Parameters.beta11init.Value = 2.4;
Parameters.beta22init.Value = 1.3;
Parameters.beta12init.Value = 0.6;
Parameters.beta21init.Value = 0.9;
Parameters = UpdateParsNoTransfToTransf(Parameters);

Cov = cov(Res.TransfThetas');
Parameters.G = Cov^-1;
Parameters.AdaptC = 0;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,1000);

save([SavePath '/TempWrongInit_Model4cstbetas' ],'Res')


load([SavePath '/Temp' num2str(11) '_Model4cstbetas' ])
for i = 1:4
    subplot(2,2,i)
    [fi,xi] = ksdensity(Res.Thetas(i,:));
    plot(xi,fi)
end


