%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DOES IT RECOGNIZE THE MODEL?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
% plot(Weigthed2,'--k','LineWidth',2)
hold off
legend('<1','1-4','5-14','15-24','25-44','45-64','65+')

for i =1:length(Data.Dates)
    disp(i)
    disp(Data.Dates{i})
end


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



Data.Observations = zeros(7,35);
Data.Observations(5,:) = Weigthed*10;
Data.Instants = [1:35]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];


SEIRModel = struct();
SEIRModel.UKF_projection = @SEIR_UKF_projection;
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);

% EKF Optimization
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];



temp = zeros(1,7);
temp(1,5) = 1;
SEIRModel.ObservationJacobian = {};
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.Instants)
    SEIRModel.ObservationJacobian{i} = temp;
    SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs*Data.Observations(5,i))^2;
end




SavePath = 'S:\Results\';
Res = load([SavePath 'Marc_ForPresEst_WithPaths.mat']);

Res = Res.Res3;
Parameters = Res.Parameters;
Parameters.ObsNoise = 0;

betas = mean(squeeze(exp(Res.Paths(:,6,:))));

tmp = mean(Res.Thetas');
Names = Parameters.Names.Estimated;
for i = 1:length(tmp)
    Parameters.(Names{i}).Value = tmp(Parameters.(Names{i}).Index);
end
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters);

Data = SEIR_CreateData(betas,Parameters,Data,SEIRModel);
plot(Data.Observations(5,:))



Res = load([SavePath 'Marc_ForPresEst_WithPaths_IBM.mat'])

Res = Res.Res3;
Parameters = Res.Parameters;

tmp = mean(Res.Thetas');
Names = Parameters.Names.Estimated;
for i = 1:length(tmp)
    Parameters.(Names{i}).Value = tmp(Parameters.(Names{i}).Index);
end
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters);

SavePath = 'S:\Results\';
Name = [SavePath 'DIC_DoesIBMRecognizeAdd.mat'];
Parameters.DiffusionType = 'IBM';

FullSEIRinference(Data,Parameters,Name,1)

SavePath = 'S:\Results\';
Name = [SavePath 'DIC_DoesIBMRecognizeAdd.mat'];
Res = load(Name);

PlotMarc(Res.Res3,8)


