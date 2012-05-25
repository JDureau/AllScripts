function Parameters = GibbsIGsampler(Name,Data,Model,Parameters)

n = length(Data.Instants);

tmp =  Model.PathLik(Data,Model,Parameters);
xis = tmp.RecordVariables(5,:);
yis = Data.Observations(5,:);
Parameters.Name.TransfValue = gamrnd(n/2+Parameters.Name.IGa,sum((xis-yis).^2)/2+b)^-1;
Parameters = UpdateParsTransfToNoTransf(Parameters);



