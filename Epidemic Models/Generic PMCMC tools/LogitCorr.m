function corr = LogitCorr(Name,Parameters,TransfValues,Values)
 
if nargin == 2
    m = Parameters.(Name).MinLim;
    M = Parameters.(Name).MaxLim;
    corr = ((M-m)*exp(Parameters.(Name).TransfValue)/((1+exp(Parameters.(Name).TransfValue))^2))^-1;     
    if or(isnan(corr),isinf(corr))
        corr = ((M-m)/((1+exp(-Parameters.(Name).TransfValue))*(1+exp(Parameters.(Name).TransfValue))))^-1;     
    end
else 
    m = Parameters.(Name).MinLim;
    M = Parameters.(Name).MaxLim;  
    corr = ((M-m)*exp(TransfValues)./((ones(size(TransfValues))+exp(TransfValues)).^2)).^-1;     
    inds = find(isnan(corr));
    corr(inds) = ((M-m)*ones(size(TransfValues))./((ones(size(inds))+exp(-TransfValues(inds))).*(ones(size(inds))+exp(TransfValues(inds))))).^-1;     
end