function Parameters = SampleParameters(Parameters)


Names = Parameters.Names.Estimated;

for i = 1:length(Names)
    if Parameters.(Names{i}).Sample
        if strcmp(Parameters.(Names{i}).TransfType,'Logit')
            if Parameters.(Names{i}).MaxLim - Parameters.(Names{i}).MinLim < Parameters.(Names{i}).Max - Parameters.(Names{i}).Min
                m = Parameters.(Names{i}).MinLim;
                M = Parameters.(Names{i}).MaxLim;                
                mean = (m+M)/2;
                sig = (M-m)/4;
                tmp = mean+randn(1,1)*sig;
                Parameters.(Names{i}).Value = m+ 0.05*(M-m) + 0.9*rand(1,1)*(M-m);
            else
                m = Parameters.(Names{i}).MinLim;
                M = Parameters.(Names{i}).MaxLim;
                sig = (M-m)/4;
                tmp = Parameters.(Names{i}).Min + rand(1,1)*(Parameters.(Names{i}).Max - Parameters.(Names{i}).Min);
                Parameters.(Names{i}).Value = max(m+0.001*sig,min(M-0.001*sig,tmp));
                Parameters.(Names{i}).Value = m+ 0.05*(M-m) + 0.9*rand(1,1)*(M-m);
            end
        elseif strcmp(Parameters.(Names{i}).TransfType,'Log')
%             tmp = Parameters.(Names{i}).MinLim + rand(1,1)*(Parameters.(Names{i}).MaxLim - Parameters.(Names{i}).MinLim);
%             Parameters.(Names{i}).Value = tmp;
        end
    end
    
%     test = 0;
%     while not(test)
%         if Parameters.(Names{i}).Max == 10^14
%             if or(isinf(Parameters.(Names{i}).MaxLim),isinf(Parameters.(Names{i}).MinLim))
%                 tmp = Parameters.(Names{i}).Value;
%                 test = 1;
%             else
%                 tmp = Parameters.(Names{i}).MinLim + rand(1,1)*(Parameters.(Names{i}).MaxLim - Parameters.(Names{i}).MinLim);
%                 test = 1;
%             end
%         else    
%             ampl = Parameters.(Names{i}).Max - Parameters.(Names{i}).Min;
%             mean = ( Parameters.(Names{i}).Max + Parameters.(Names{i}).Min )/2;
%             tmp = randn(1,1)*ampl/2 + mean;
%             if and(tmp > Parameters.(Names{i}).MinLim, tmp < Parameters.(Names{i}).MaxLim) 
%                 test = 1;
%             end
%         end
%     end
%     Parameters.(Names{i}).Value = tmp;
end

Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
