function []  = TellParValues(Par)

Names = Par.Names.Estimated;

for i = 1:length(Names)
    disp([Names{i} ': ' num2str(Par.(Names{i}).Value)]);
end