function Par = TransfToNoTransf(Par)

Names = Par.Names.Estimated;

for i = 1:length(Names)
    Par.(Names{i}).Value = max(Par.(Names{i}).MinLim + 0.001, min(Par.(Names{i}).MaxLim - 0.001, Par.(Names{i}).InvTransf(Par.(Names{i}).TransfValue,Par.(Names{i}).MinLim,Par.(Names{i}).MaxLim)));
    Par.(Names{i}).TransfValue = Par.(Names{i}).Transf(Par.(Names{i}).Value,Par.(Names{i}).MinLim,Par.(Names{i}).MaxLim);
end