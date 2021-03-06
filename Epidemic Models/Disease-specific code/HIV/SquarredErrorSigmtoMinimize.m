function mLogPost = SquarredErrorSigmtoMinimize(x,NamesEst,Data,Parameters,HIVModel)


for i = 1:length(NamesEst)
    Parameters.(NamesEst{i}).TransfValue = x(i);
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsTransfToNoTransf(Parameters);

TempSim = SimulateSigmoidDetermEpid(Data,Parameters,HIVModel);

if not(isreal(TempSim.Observations))
    TellParsValues(TempSim.Parameters)
    die
end


LogLik = 0;
LogPrior = 0;
LogCorr = 0;

for i = 2:length(Data.ObservedVariables);
%     LogLik = LogLik + max(-700,log(normpdf(Data.Observations(Data.ObservedVariables(i),i),TempSim.Observations(Data.ObservedVariables(i),i),sqrt(Data.Observations(Data.ObservedVariables(i),i)*(100-Data.Observations(Data.ObservedVariables(i),i))/400))));
    LogLik = LogLik + max(-700,log(binopdf(round(425*Data.Observations(Data.ObservedVariables{i}(1),i)),425,TempSim.Observations(Data.ObservedVariables{i}(1),i)/100)));
    if length(Data.ObservedVariables{i})==2
        corr = TempSim.Observations(Data.ObservedVariables{i}(2),i);
        LogLik = LogLik + max(-700,log(binopdf(round(425*Data.Observations(Data.ObservedVariables{i}(2),i)),425,Parameters.Rho.Value*corr)));
    end
end

for i = 1:length(NamesEst)
    tmp = Parameters.(NamesEst{i}).Prior(NamesEst{i},Parameters);
    LogPrior = LogPrior + log(tmp);
    temp = Parameters.(NamesEst{i}).CorrFunct(NamesEst{i},Parameters);
    LogCorr = LogCorr + log(temp);
%     NamesEst{i}
%     'prior'
%     log(tmp)
%     'corr'
%     log(temp)
end
mLogPost = -(LogLik + LogPrior - LogCorr);
title(mLogPost)
pause(0.01)




