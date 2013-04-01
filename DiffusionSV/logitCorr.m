function corr = logitCorr(Name,Parameters)
 
    m = Parameters.(Name).MinLim;
    M = Parameters.(Name).MaxLim;
    corr = ((M-m)*exp(Parameters.(Name).TransfValue)/((1+exp(Parameters.(Name).TransfValue))^2));     
    if or(isnan(corr),isinf(corr))
        corr = ((M-m)*exp(-Parameters.(Name).TransfValue)/((1+exp(-Parameters.(Name).TransfValue))^2));     
    end
    corr = max(corr,eps);
