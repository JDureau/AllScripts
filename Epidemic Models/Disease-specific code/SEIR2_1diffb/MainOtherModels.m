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

%% Constant \beta model.


SavePath = 'S:\Results\';

Res = load([SavePath 'Marc_FirstEst_WithPaths.mat']);

ResPost = Res.Res3;


Parameters = ResPost.Parameters;
Parameters.Begining.Value = 1.1;
Parameters.Begining.MinLim = 0;
Parameters.Begining.MaxLim = 8;
Parameters.Begining.Min = -10^14;
Parameters.Begining.Max = 10^14;
Parameters.Begining.Estimated = 1;
Parameters.Begining.TransfType = 'Logit';
Parameters.FirstBeta.Value = 4;
Parameters.FirstBeta.Min = -10^14;
Parameters.FirstBeta.Max = 10^14;
Parameters.FirstBeta.Estimated = 1;
Parameters.FirstBeta.TransfType = 'Log';
Parameters.SecondBeta.Value = 4;
Parameters.SecondBeta.Min = -10^14;
Parameters.SecondBeta.Max = 10^14;
Parameters.SecondBeta.Estimated = 1;
Parameters.SecondBeta.TransfType = 'Log';
Parameters.ThirdBeta.Value = 2;
Parameters.ThirdBeta.Min = -10^14;
Parameters.ThirdBeta.Max = 10^14;
Parameters.ThirdBeta.Estimated = 1;
Parameters.ThirdBeta.TransfType = 'Log';
Parameters.FourthBeta.Value =3;
Parameters.FourthBeta.Min = -10^14;
Parameters.FourthBeta.Max = 10^14;
Parameters.FourthBeta.Estimated = 1;
Parameters.FourthBeta.TransfType = 'Log';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)

randind = ceil(rand(1,1)*5000);
Parameters.Pars = (squeeze((ResPost.Paths(randind,6,:))));
aimpath =  (squeeze((ResPost.Paths(randind,5,:))))';
Names = Parameters.Names.Estimated;
NamesPost = ResPost.Parameters.Names.Estimated;
for i = 1:length(NamesPost)
    ind = ResPost.Parameters.(NamesPost{i}).Index;
    Parameters.(NamesPost{i}).TransfValue = ResPost.TransfThetas(ind,randind);
end
Parameters.FirstBeta.TransfValue = Parameters.Pars(1);
Parameters.SecondBeta.TransfValue = mean(Parameters.Pars(2:sum(ResPost.Data.NbComputingSteps(1:8))));
Parameters.ThirddBeta.TransfValue = mean(Parameters.Pars(sum(ResPost.Data.NbComputingSteps(1:8)):sum(ResPost.Data.NbComputingSteps(1:14))));
Parameters.FourthBeta.TransfValue = mean(Parameters.Pars(sum(ResPost.Data.NbComputingSteps(1:14)):sum(ResPost.Data.NbComputingSteps)));
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsTransfToNoTransf(Parameters);
Parameters.LengthPaths = length(Parameters.Pars);
Temp = SEIR_PathLik(ResPost.Data, SEIRModel,Parameters);

subplot(2,1,1)
plot(squeeze(Temp.CompletePaths(1,6,:)))
hold on
plot(Parameters.Pars,'g')
hold off
subplot(2,1,2)
inds = cumsum(ResPost.Data.NbComputingSteps);
plot(aimpath(max(1,inds))','g')
hold on
plot(squeeze(Temp.CompletePaths(1,5,max(1,inds))))
hold off



Cov = zeros(length(Names));
inds = [];
for i = 1:7
    inds(i) = Parameters.(Names{i}).Index;
end

tmpcov = cov(ResPost.TransfThetas(inds,:)');
Cov(1:7,1:7) = tmpcov;
% Cov(1:7,1:7) = eye(7) + diag(diag(tmpcov));
for i = 8:12
    Cov(i,i) = 0.8*abs(Parameters.(Names{i}).TransfValue);
end

Parameters.ModelType = 'MCMC';


SEIRModel = struct();
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
SEIRModel.BuildPath = @SEIR_BuildPathConstantBeta;
SEIRModel.PathLik = @SEIR_PathLik;
Parameters = SEIRModel.InitializeParameters(Parameters);

Parameters.MCMCType = 'Rand';
Parameters.G = Cov^-1;
Parameters.Epsil = 0.1;
Parameters.NoPaths = 0;
Parameters.PathsToKeep = 1:6;
Parameters.NbParticules = 1;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
Res = RunEstimationMethod(Data, SEIRModel,Parameters,Res.TempPar,1000);

    
PlotMarc(Res,8);

