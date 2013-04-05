function Par = DefineIndexes(Par)


Names = Par.Names.All;

tmp = {};
cpt = 1;
for i = 1:length(Names)
    if Par.(Names{i}).Estimated     
        Par.(Names{i}).Index = cpt;
        tmp{cpt} = Names{i};
        cpt = cpt+1;
    end
end
Par.Names.Estimated = tmp;


Names = Par.Names.Estimated;
for i = 1:length(Names)
    Names{i}
    if strcmp(Par.(Names{i}).TransfType, 'Logit')
%         if Par.(Names{i}).Max == 10^14
            Par.(Names{i}).Prior = @UnifLogitPrior;
            Par.(Names{i}).DerPrior = @DerUnifLogitPrior;
            Par.(Names{i}).CorrFunct = Par.(Names{i}).Corr;
%         else
%             Par.(Names{i}).Prior = @NormalLogitPrior;
%         end
    elseif strcmp(Par.(Names{i}).TransfType, 'Log')
        Par.(Names{i}).Prior = @NormalLogPrior;    
        Par.(Names{i}).DerPrior = @DerNormalLogPrior;
        Par.(Names{i}).CorrFunct = @logCorr;
        Par.(Names{i}).MinLim = 0;
        Par.(Names{i}).MaxLim = 0;
    elseif strcmp(Par.(Names{i}).TransfType, 'Id')
        Par.(Names{i}).Prior = @NormalLogPrior;    
        Par.(Names{i}).DerPrior = @DerNormalLogPrior;
        Par.(Names{i}).CorrFunct = @idCorr;
        Par.(Names{i}).MinLim = 0;
        Par.(Names{i}).MaxLim = 0;
    end
%     if strcmp(Names{i}, 'rho')
%         Par.(Names{i}).Prior = @RhoPrior;
%         Par.(Names{i}).DerPrior = @DerRhoPrior;
%         Par.(Names{i}).CorrFunct = Par.(Names{i}).Corr;
%     end
    if strcmp(Names{i}, 'sigma_X')
        Par.(Names{i}).Prior = @SigmaPrior2;
        Par.(Names{i}).DerPrior = @DerSigmaPrior2;
        Par.(Names{i}).CorrFunct = Par.(Names{i}).Corr;
    end
%     if strcmp(Names{i}, 'kappa')
%         Par.(Names{i}).Prior = @KappaPrior;
%         Par.(Names{i}).DerPrior = @DerKappaPrior;
%         Par.(Names{i}).CorrFunct = Par.(Names{i}).Corr;
%     end
end
