function LogPrior = ComputeLogPriorZ_Full(Par)

LogPrior = 0;

% % Z
% LogPrior = LogPrior - 0.5*(Z')*Z;


% PARAMETERS
% we make the assumption here that uniform priors are used for parameters
Names = Par.Names.Estimated;
if Par.Prior
    if Par.GradCorr
        for j = 1:length(Names)
            temp = Par.(Names{j}).Prior(Names{j},Par);
%             disp([Names{j} '   ' num2str(log(temp) + log(Par.(Names{j}).Corr(Names{j},Par)))])
            LogPrior = LogPrior + log(temp) + log(Par.(Names{j}).Corr(Names{j},Par));
        end
    end
end
