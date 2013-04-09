function prior = NormalLogPrior(Name,Parameters,TransfValues,Values)

if nargin == 2
    if Parameters.(Name).StdPrior>10^10
        prior = 1;
    else
        prior = normpdf(Parameters.(Name).Value,Parameters.(Name).MeanPrior,Parameters.(Name).StdPrior)/Parameters.(Name).CLim;     
    end
    if strcmp(Name,'betainit')
        prior = 1;
    end
else 
    prior = normpdf(Values,Parameters.(Name).MeanPrior,Parameters.(Name).StdPrior)/Parameters.(Name).CLim;     
end