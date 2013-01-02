function corr = logitCorr(Name,Parameters)
 
    m = Parameters.(Name).MinLim;
    M = Parameters.(Name).MaxLim;
    corr = ((M-m)*exp(Parameters.(Name).TransfValue)/((1+exp(Parameters.(Name).TransfValue))^2))^-1;     
    if or(isnan(corr),isinf(corr))
        corr = ((M-m)/((1+exp(-Parameters.(Name).TransfValue))*(1+exp(Parameters.(Name).TransfValue))))^-1;     
    end
