function prior = DerNormalLogPrior(Name,Parameters,TransfValues,Values)

if nargin == 2
    if 1%Parameters.(Name).StdPrior>10^10
        prior = 0;
        
    else
        prior = normpdf(Parameters.(Name).Value,Parameters.(Name).MeanPrior,Parameters.(Name).StdPrior)/Parameters.(Name).CLim;     
    end
    if strcmp(Name,'betainit')
        prior = 1;
    end
else 
    prior = normpdf(Values,Parameters.(Name).MeanPrior,Parameters.(Name).StdPrior)/Parameters.(Name).CLim;     
end