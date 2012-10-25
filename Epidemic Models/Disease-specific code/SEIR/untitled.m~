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


Parameters.SigmaRW11.Value = 0;
Parameters.SigmaRW22.Value = 0;
Parameters.beta11init.Value = 2;
Parameters.beta22init.Value = 1;
Parameters.beta12init.Value = 0.3;
Parameters.adultsmult.Value = 0.6;
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

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.beta11init.Estimated = 1;
Parameters.beta22init.Estimated = 1;
Parameters.beta12init.Estimated = 1;
Parameters.adultsmult.Estimated = 1;


Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);



