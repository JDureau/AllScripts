function Parameters = DefineEstimatedParametersIndexes(Parameters)

% Get all names
names = fieldnames(Parameters);

Names = {};
EstimatedNames = {};
ind = 1;
for i = 1:length(names)
    try
        if Parameters.(names{i}).Estimated
            EstimatedNames{end+1} = names{i};
            Parameters.(names{i}).Index = ind;
            ind = ind+1;
        end
        if isstruct(Parameters.(names{i}))
            Names{end+1} = names{i};
            if not(isfield(Parameters.(names{i}),'PlotInv'))
                Parameters.(names{i}).PlotInv = 0;
            end
            if strcmp(Parameters.(names{i}).TransfType, 'Log');
                Parameters.(names{i}).MinLim = 0;
                Parameters.(names{i}).MaxLim = Inf;
            end
            Parameters.(names{i}).MeanLim = (Parameters.(names{i}).MaxLim+Parameters.(names{i}).MinLim)/2;
            Parameters.(names{i}).StdLim = (Parameters.(names{i}).MaxLim-Parameters.(names{i}).MinLim)/2;
            Parameters.(names{i}).MeanPrior = (Parameters.(names{i}).Max+Parameters.(names{i}).Min)/2;
            Parameters.(names{i}).StdPrior = (Parameters.(names{i}).Max-Parameters.(names{i}).Min)/4;
            Parameters.(names{i}).CLim = max(eps,normcdf(Parameters.(names{i}).MaxLim,Parameters.(names{i}).MeanPrior,Parameters.(names{i}).StdPrior) - normcdf(Parameters.(names{i}).MinLim,Parameters.(names{i}).MeanPrior,Parameters.(names{i}).StdPrior));
            if not(isfield(Parameters.(names{i}),'Init'))
                Parameters.(names{i}).Init = 0;
            end
            try if not(Parameters.(names{i}).Sample)
                    Parameters.(names{i}).Sample = 0;
                end
            catch
                Parameters.(names{i}).Sample = 0;
            end
        end
    end
end
Parameters.Names.All = Names;
Parameters.Names.Estimated = EstimatedNames;
Parameters.NbParsEstimated = ind - 1;


tmpInit = 0;
EstInit = {};
tmpNotInit = 0;
EstNotInit = {};
names = EstimatedNames;
for i = 1:ind-1
%     names{i}
    if  Parameters.(names{i}).Init
        tmpInit = tmpInit + 1;
        EstInit{end+1} = names{i};
        Parameters.(names{i}).IndexInit = tmpInit;
    else
        tmpNotInit = tmpNotInit + 1;
        EstNotInit{end+1} = names{i};
        Parameters.(names{i}).IndexNotInit = tmpNotInit;
    end
end
Parameters.NbParsEstimatedInit = tmpInit;
Parameters.Names.EstimatedInit = EstInit;
Parameters.NbParsEstimatedNotInit = tmpNotInit;
Parameters.Names.EstimatedNotInit = EstNotInit;