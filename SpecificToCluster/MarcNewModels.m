function [] = MarcNewModels(IndModel)

IndModel = IndModel + 1;


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


if IndModel == 1
    Data.Observations = zeros(12,35);
    Data.Observations(5,:) = Weighted*10;
else
    Data.Observations = zeros(12,35);
    Data.Observations(9,:) = Students*10;
    Data.Observations(10,:) = Adults*10;
end

plot(Students,'g')
hold on
plot(Adults)
plot(Weighted,'k')
plot(Adults+Students,'r')
hold off

DiffType = 'Add';
ObsType = 'Fixed';
Name = '';

FullSEIRinference(Data,DiffType,ObsType,Name,IndModel)


