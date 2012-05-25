function [] = TellParsValues(Parameters)

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    try
        if strcmp(Parameters.(Names{i}).TransfType,'Logit')
            if Parameters.(Names{i}).MaxLim - Parameters.(Names{i}).MinLim < Parameters.(Names{i}).Max - Parameters.(Names{i}).Min
                disp([Names{i} ': ' num2str(Parameters.(Names{i}).Value) ' (' num2str((Parameters.(Names{i}).Value-Parameters.(Names{i}).MinLim)/(Parameters.(Names{i}).MaxLim-Parameters.(Names{i}).MinLim)*100,2) '%)'])
            else
                disp([Names{i} ': ' num2str(Parameters.(Names{i}).Value) ' (' num2str((Parameters.(Names{i}).Value-Parameters.(Names{i}).Min)/(Parameters.(Names{i}).Max-Parameters.(Names{i}).Min)*100,2) '%)'])
            end
        else
            disp([Names{i} ': ' num2str(Parameters.(Names{i}).Value) ' (' num2str((Parameters.(Names{i}).Value-Parameters.(Names{i}).Min)/(Parameters.(Names{i}).Max-Parameters.(Names{i}).Min)*100,2) '%)'])
        end
    catch
        disp([Names{i} ': ' num2str(Parameters.(Names{i}).Value) ' (' num2str((Parameters.(Names{i}).Value-Parameters.(Names{i}).MinLim)/(Parameters.(Names{i}).MaxLim-Parameters.(Names{i}).MinLim)*100,2) '%)'])
    end
end