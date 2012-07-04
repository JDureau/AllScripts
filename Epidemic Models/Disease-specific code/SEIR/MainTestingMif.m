%% Testing deterministic MIF



cd('/Users/dureaujoseph/Dropbox/AllScriptsGit/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR'])



SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc';

A = load([SavePath '/andre_estimates_31_01.txt']);

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
Students = sum((A(:,1:4)*diag(PopWeigths(1:4).^-1)*diag(PopProps(1:4)/100)*100000)');
Adults = sum((A(:,5:7)*diag(PopWeigths(5:7).^-1)*diag(PopProps(5:7)/100)*100000)');

plot(Students,'g')
hold on
plot(Adults,'b')
hold off

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



load([SavePath '/ParametersSEIR.mat']);

Parameters.MIFL = 5;
Parameters.MIFNbIts = 150;
Parameters.MIFa = 0.99;
Parameters.MIFb = 1;

Parameters.SigmaObs.Value = 0.1;
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.SMC_projection = @SEIR_SMC_projection;
SEIRModel.LikFunction = 'normpdf(log(Variables(:,5)),transpose(log(coeff*Data.Observations(5,IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),sqrt(log(Parameters.SigmaObs.Value^2+1)))';%Parameters.SigmaObs.Value)';


Res = RunDetermMIF(Data, SEIRModel, Parameters)

k = ceil(sqrt(length(Parameters.Names.Estimated)));
subplot(k+1,k,1:k)
plot(Res.LogLiks)
for i = 1:length(Parameters.Names.Estimated)
    subplot(k+1,k,k+i)
    plot(Res.Thetas(i,:))
end







