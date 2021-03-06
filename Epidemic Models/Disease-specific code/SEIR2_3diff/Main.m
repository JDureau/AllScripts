

cd('/Users/dureaujoseph/Dropbox/Taf/AllScripts/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR2_1diff'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR2_2diff'])


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

PropStudents = sum(PopWeigths(1:4))/sum(PopWeigths);
PropAdults = sum(PopWeigths(5:7))/sum(PopWeigths);

tmp = (A(:,1:4)*diag(PopWeigths(1:4).^-1)*diag(PopWeigths(1:4)))/sum(PopWeigths(1:4))*PropStudents*100000;
Students = sum(tmp');
tmp = (A(:,5:7)*diag(PopWeigths(5:7).^-1)*diag(PopWeigths(5:7)))/sum(PopWeigths(5:7))*PropStudents*100000;
Adults = sum(tmp');



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





Data.Observations = zeros(12,35);
Data.Observations(9,:) = Students*10;
Data.Observations(10,:) = Adults*10;

DiffType = 'Add';
ObsType = 'Fixed';
Parameters.StructModel = 1;
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc';


Name = [SavePath '/MarcData_StructModel2_2diff.mat'];


FullSEIRinference(Data,Difftype,Obstype,Name)

Parameters.SigmaRW11.Value = 0.1;
Parameters.SigmaRW22.Value = 0.1;
Parameters.beta11init.Value = 0.1;
Parameters.beta22init.Value = 0.1;
Parameters.beta12init.Value = 0.1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);






