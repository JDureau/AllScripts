function corr = logCorrDer(Name,Parameters)
 
corr = exp(Parameters.(Name).TransfValue);     
