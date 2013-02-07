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
    end
    if strcmp(Names{i}, 'rho')
        Par.(Names{i}).Prior = @RhoPrior;
        Par.(Names{i}).DerPrior = @DerRhoPrior;
        Par.(Names{i}).CorrFunct = Par.(Names{i}).Corr;
    end
end
