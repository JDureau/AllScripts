function Parameters = UpdateParsNoTransfToTransf(Parameters)

Names = Parameters.Names.Estimated;

for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = eval(Parameters.(Names{i}).Transf);
end

CheckParametersGen(Parameters)