function Parameters = DefineTransfFunctions(Parameters)

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    if strcmp(Parameters.(Names{i}).TransfType, 'Logit')
        Parameters.(Names{i}).Transf = 'LogitTransf(Parameters.(Names{i}).Value,Parameters.(Names{i}).MinLim,Parameters.(Names{i}).MaxLim)';
        Parameters.(Names{i}).InvTransf = 'InvLogitTransf(Parameters.(Names{i}).TransfValue,Parameters.(Names{i}).MinLim,Parameters.(Names{i}).MaxLim)';
        Parameters.(Names{i}).CorrFunct = @LogitCorr;
    elseif strcmp(Parameters.(Names{i}).TransfType, 'Log')
        Parameters.(Names{i}).Transf = 'log(Parameters.(Names{i}).Value)';
        Parameters.(Names{i}).InvTransf = 'exp(Parameters.(Names{i}).TransfValue)';
        Parameters.(Names{i}).CorrFunct = @LogCorr;
    else
        disp('problem, unknown transformation')
    end
end
    