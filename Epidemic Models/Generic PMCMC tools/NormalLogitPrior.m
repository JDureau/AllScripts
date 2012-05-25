function prior = NormalLogitPrior(Name,Parameters,TransfValues,Values)


if nargin == 2
    prior = normpdf(Parameters.(Name).Value,Parameters.(Name).MeanPrior,Parameters.(Name).StdPrior)/Parameters.(Name).CLim;     
   
else 
    prior = normpdf(Values,Parameters.(Name).MeanPrior,Parameters.(Name).StdPrior)/Parameters.(Name).CLim;     
end
