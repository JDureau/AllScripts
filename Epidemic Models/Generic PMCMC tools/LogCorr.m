function corr = LogCorr(Name,Parameters,TransfValues,Values)
 
if nargin == 2
    corr = exp(-Parameters.(Name).TransfValue);     
else 
    corr = exp(-TransfValues);     
end