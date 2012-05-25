function DataExtracted = ExtractGoogleTSeries(Data,Parameters,region,year)

% This function will extract the data starting July of that year untill
% following July.

DataExtracted.Observations = zeros(6,52);
DataExtracted.ObservationsDates = {};
indregion = find(strcmp(Data.regions,region));
inddate = 1;
notfound = 1;
while notfound
    if str2num(Data.Dates{inddate}(1:4)) == year
        if str2num(Data.Dates{inddate}(6:7)) > 6
            notfound = 0;
        end
    end
    inddate = inddate + 1;
end

for i = 0:51
    DataExtracted.Observations(5,i+1) =  Data.NewCases{i+inddate}(indregion);
    DataExtracted.ObservationsDates{i+1} = Data.Dates{i+inddate};
end
DataExtracted.Instants = [1:52]*7/Parameters.ComputationTStep;
DataExtracted.ObservedVariables = 5*ones(1,length(DataExtracted.Instants));
DataExtracted.NbComputingSteps = [0 diff(DataExtracted.Instants)];
DataExtracted.Dates = Data.Dates;