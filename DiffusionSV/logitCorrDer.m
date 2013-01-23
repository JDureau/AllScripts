function corr = logitCorrDer(Name,Parameters)
 
    m = Parameters.(Name).MinLim;
    M = Parameters.(Name).MaxLim;
    corr = (M-m)*(1-exp(2*Parameters.(Name).TransfValue))*exp(Parameters.(Name).TransfValue)/((1+exp(Parameters.(Name).TransfValue))^4);
%    corr = (exp(2*Parameters.(Name).TransfValue)-1)/((M-m)*exp(Parameters.(Name).TransfValue));     
%     if or(isnan(corr),isinf(corr))
%  %       corr = (1-exp(-2*Parameters.(Name).TransfValue))/((M-m)*exp(-Parameters.(Name).TransfValue)); 
%     end
