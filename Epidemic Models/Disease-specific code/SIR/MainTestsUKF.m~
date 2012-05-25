% UKF tests on SIR



cd('/Users/dureaujoseph/Dropbox/Taf/tmpsaveKalman/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/SIR'])
addpath([pwd '/Epidemic Models/Disease-specific code/Linear'])


Data = struct();
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc/';
load([SavePath 'ParametersSEIR.mat']);
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).Estimated = 0;
end
Parameters.gammam1.Value = 1;
Parameters.SigmaRW.Value = 0.05;
Parameters.gammam1.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

%% Prevalence data

Data.Observations = zeros(3,20);

mpred = [100000 1 log(3)]';
Parameters.PopSize = 100000;
TStep = 1/100;
betas = 1.5+0.1*(1:20*100)*TStep;
for i = 2:20
    for j = 1:100
        beta = betas((i-1)*100+j);
        mtemp = mpred;
        mpred(1) = mpred(1) + (-beta*mtemp(1)*mtemp(2)/Parameters.PopSize)*TStep;
        mpred(2) = mpred(2) + ( beta*mtemp(1)*mtemp(2)/Parameters.PopSize- Parameters.gammam1.Value^(-1)*mtemp(2))*TStep;
    end
    Data.Observations(2,i) = mpred(2);
end
plot(Data.Observations(2,:))
Data.ObservedVariables = 2*ones(1,20);
Parameters.ComputationTStep = 0.01;
Data.Instants = (1:20)/Parameters.ComputationTStep;
Data.NbComputingSteps = [0 diff(Data.Instants)];
Parameters.NbVariables = 4;

clf
plot(Data.Instants,Data.Observations(2,:))


SIRModel.InitializeParameters = @SIR_Initialize;
SIRModel.UKF_projection = @SIR_UKF_projection;
SIRModel.EKF_projection = @SIR_EKF_projection;
SIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.ObservedVariables)-1
    SIRModel.ObservationMeasurementNoise{i+1} = 0.005*Data.Observations(2,i+1);
end

clf
Result = EstimationEKFGen(Data, SIRModel, Parameters);


% clf
Result = EstimationUKFGen(Data, SIRModel, Parameters);


Parameters.Correction = 0;
Parameters = KalmOpt(Parameters,Data,SIRModel,2000);


%% Incidence data

Data.Observations = zeros(4,20);

mpred = [100000 1 0 log(2)]';
Parameters.PopSize = 100000;
TStep = 1/100;
betas = 1.5+0.1*(1:20*100)*TStep;
record = zeros(4,length(betas));
for i = 2:20
    mpred(3) = 0;
    for j = 1:100
        beta = betas((i-1)*100+j);
        mtemp = mpred;
        mpred(1) = mpred(1) + (-beta*mtemp(1)*mtemp(2)/Parameters.PopSize)*TStep;
        mpred(2) = mpred(2) + ( beta*mtemp(1)*mtemp(2)/Parameters.PopSize- Parameters.gammam1.Value^(-1)*mtemp(2))*TStep;
        mpred(3) = mpred(3) + ( beta*mtemp(1)*mtemp(2)/Parameters.PopSize)*TStep;
        mpred(4) = log(beta);
        record(:,(i-1)*100+j) = mpred;
    end
    Data.Observations(3,i) = mpred(3);
end
plot(Data.Observations(3,:))
Data.ObservedVariables = 3*ones(1,20);
Parameters.ComputationTStep = 0.01;
Data.Instants = (1:20)/Parameters.ComputationTStep;
Data.NbComputingSteps = [0 diff(Data.Instants)];
Parameters.NbVariables = 4;

clf
plot(Data.Instants,Data.Observations(3,:))


SIRModel.InitializeParameters = @SIR_Initialize;
SIRModel.UKF_projection = @SIR_UKF_projection;
SIRModel.EKF_projection = @SIR_EKF_projection;
SIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.ObservedVariables)-1
    SIRModel.ObservationMeasurementNoise{i+1} = (0.01*Data.Observations(3,i+1))^2;
end

clf
Result = EstimationEKFGen(Data, SIRModel, Parameters);


% clf
Result = EstimationUKFGen(Data, SIRModel, Parameters);
plot(100:200,record(3,100:200),'k')

Result = EstimationSRUKFGen(Data, SIRModel, Parameters);


Parameters.Correction = 0;
Parameters = KalmOpt(Parameters,Data,SIRModel,2000);

%% Linear data

Data.Observations = zeros(4,20);

mpred = [1 1 1 1]';
TStep = 1/100;
betas = 1+0.1*(1:20/TStep)*TStep;
for i = 1:20
    for j = 1:100
        mpred(1) = mpred(1)+ (mpred(1)*0.05+mpred(2)*0.06+mpred(3)*0.07+mpred(4)*0.08)*TStep;
        mpred(2) = mpred(2)+ (mpred(1)*0.09+mpred(2)*0.1 +mpred(3)*0.11+mpred(4)*0.12)*TStep;
        mpred(3) = mpred(3)+ (mpred(1)*0.13+mpred(2)*0.14+mpred(3)*0.15+mpred(4)*0.16)*TStep;
        mpred(4) = mpred(4)+ (mpred(1)*0.17+mpred(2)*0.18+mpred(3)*0.19+mpred(4)*0.20)*TStep;
    end
    Data.Observations(3,i) = mpred(3)*( 1+ 0.01*randn(1,1));
end
plot(Data.Observations(3,:))
Data.ObservedVariables = 3*ones(1,20);
Parameters.ComputationTStep = 0.005;
Data.Instants = (1:20)/Parameters.ComputationTStep;
Data.NbComputingSteps = [0 diff(Data.Instants)];
Parameters.NbVariables = 4;

clf
plot(Data.Instants,Data.Observations(3,:))


LinearModel.InitializeParameters = @Linear_Initialize;
LinearModel.UKF_projection = @Linear_UKF_projection;
LinearModel.EKF_projection = @Linear_EKF_projection;
LinearModel.SRUKF_projection = @Linear_SRUKF_projection;
LinearModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.ObservedVariables)-1
    LinearModel.ObservationMeasurementNoise{i+1} = 0.01*Data.Observations(3,i+1);
end

Parameters.SigmaRW.Value = 0.9;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

clf
Result = EstimationEKFGen(Data, LinearModel, Parameters);

% clf
Result = EstimationUKFGen(Data, LinearModel, Parameters);

Result = EstimationSRUKFGen(Data, LinearModel, Parameters);


Parameters.Correction = 0;
Parameters = KalmOpt(Parameters,Data,SIRModel,2000);




%% SEIR

DataPath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/HPA';

A = load([DataPath '/andre_estimates_31_01.txt']);
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


Parameters.ComputationTStep=0.01;
Data.Observations = zeros(7,35);
Data.Observations(5,:) = Weigthed*10;
DiffType = 'Add';
ObsType = 'Estimated';
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc/';
load([SavePath 'ParametersSEIR.mat']);
Data.Instants = [1:size(Data.Observations,2)]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];

SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc';



clf
plot(Data.Instants,Data.Observations(5,:))

Parameters.NbVariables = 6;


SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.UKF_projection = @SEIR_UKF_projection;
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.ObservedVariables)-1
    SEIRModel.ObservationMeasurementNoise{i+1} = 0.5*Data.Observations(5,i+1);
end

clf
Result = EstimationUKFGen(Data, SEIRModel, Parameters);


