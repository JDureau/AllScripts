function [] = MarcSensitivityAnalysis(IndPar,IndShift)

% IndShift : +10 +20 -10 -20



cd('/users/ecologie/dureau/src/AllScripts/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR'])


SavePath = '/users/ecologie/dureau/src/AllData/ResultsMarc/';
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR2_1diff'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR2_1diffb'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR2_2diff'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR2_2diffb'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR2_3diff'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR2_4diff'])



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

PopWeigths = [667600,2461800,5904100,6862500,14417400,12847800,7929300];
PopProps = PopWeigths/sum(PopWeigths)*100;
Weighted = sum((A(:,1:7)*diag(PopWeigths.^-1)*diag(PopProps/100)*100000)');

PropStudents = sum(PopWeigths(1:3))/sum(PopWeigths);
PropAdults   = sum(PopWeigths(4:7))/sum(PopWeigths);


Students = sum((A(:,1:3)*diag(PopWeigths(1:3).^-1)*diag(PopProps(1:3)/100)*100000)');
Adults   = sum((A(:,4:7)*diag(PopWeigths(4:7).^-1)*diag(PopProps(4:7)/100)*100000)');


Data.Observations = zeros(12,35);
Data.Observations(5,:) = Weighted*10;

plot(Students,'g')
hold on
plot(Adults)
plot(Weighted,'k')
plot(Adults+Students,'r')
hold off

DiffType = 'Add';
ObsType = 'Estimated';
Name = '';


load([SavePath 'ParametersSEIR.mat']);
pars = {'RInitProp','km1','gammam1'};
shift = [0 0.1 0.2 -0.1 -0.2];
m = Parameters.(pars{IndPar}).Min;
M = Parameters.(pars{IndPar}).Max;
delta = M-m;
Parameters.(pars{IndPar}).Min = max(m+shift(IndShift+1)*delta,Parameters.(pars{IndPar}).MinLim);
Parameters.(pars{IndPar}).Max = min(M+shift(IndShift+1)*delta,Parameters.(pars{IndPar}).MaxLim);
Parameters.(pars{IndPar}).Value = (Parameters.(pars{IndPar}).Min+Parameters.(pars{IndPar}).Max)/2;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

save([SavePath 'ParametersSEIR_' num2str(IndPar) '_' num2str(IndShift) '.mat'],'Parameters');

Data.TwistingPriors = 1;
Data.NameTwistedPriors =  ['/ParametersSEIR_' num2str(IndPar) '_' num2str(IndShift) '.mat'];
Data.NameToSave = ['MarcData_Sensitivity_' pars{IndPar} '_' num2str(IndShift) '.mat'];

FullSEIRinference(Data,DiffType,ObsType,Name,1)







