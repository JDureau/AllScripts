function Par = NoTransfToTransf(Par)

Names = Par.Names.Estimated;

for i = 1:length(Names)
    Par.(Names{i}).TransfValue = Par.(Names{i}).Transf(Par.(Names{i}).Value,Par.(Names{i}).MinLim,Par.(Names{i}).MaxLim);
end