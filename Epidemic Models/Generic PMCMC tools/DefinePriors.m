function Parameters = DefinePriors(Parameters)

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    if strcmp(Parameters.(Names{i}).TransfType, 'Logit')
        if Parameters.(Names{i}).Max == 10^14
            Parameters.(Names{i}).Prior = @UnifLogitPrior;
        else
            Parameters.(Names{i}).Prior = @NormalLogitPrior;
        end
    elseif strcmp(Parameters.(Names{i}).TransfType, 'Log')
        Parameters.(Names{i}).Prior = @NormalLogPrior;        
    end
end

% Names = Parameters.Names.Estimated;
% for i = 1:length(Names)
%     if strcmp(Parameters.(Names{i}).TransfType, 'Logit')
%         if Parameters.(Names{i}).Max == 10^14
%             Parameters.(Names{i}).Prior = 'exp(ParametersForPriors.(Names{i}).TransfValue)/((1+exp(ParametersForPriors.(Names{i}).TransfValue))^2)*1/(Parameters.(Names{i}).MaxLim-Parameters.(Names{i}).MinLim)';
%             Parameters.(Names{i}).Prior2 = '1/((1+exp(-ParametersForPriors.(Names{i}).TransfValue))^2)*1/(Parameters.(Names{i}).MaxLim-Parameters.(Names{i}).MinLim)';
%         else
%             Parameters.(Names{i}).Prior = 'exp(ParametersForPriors.(Names{i}).TransfValue)/((1+exp(ParametersForPriors.(Names{i}).TransfValue))^2)*normpdf(ParametersForPriors.(Names{i}).Value,ParametersForPriors.(Names{i}).MeanPrior,ParametersForPriors.(Names{i}).StdPrior)/ParametersForPriors.(Names{i}).CLim';
%             Parameters.(Names{i}).Prior2 = '1/((1+exp(-ParametersForPriors.(Names{i}).TransfValue))^2)*normpdf(ParametersForPriors.(Names{i}).Value,ParametersForPriors.(Names{i}).MeanPrior,ParametersForPriors.(Names{i}).StdPrior)/ParametersForPriors.(Names{i}).CLim';
%         end
%     elseif strcmp(Parameters.(Names{i}).TransfType, 'Log')
%         Parameters.(Names{i}).Prior = 'exp(ParametersForPriors.(Names{i}).TransfValue)*normpdf(ParametersForPriors.(Names{i}).Value,ParametersForPriors.(Names{i}).MeanPrior,ParametersForPriors.(Names{i}).StdPrior)/ParametersForPriors.(Names{i}).CLim';        
%     end
% end

