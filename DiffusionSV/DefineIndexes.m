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