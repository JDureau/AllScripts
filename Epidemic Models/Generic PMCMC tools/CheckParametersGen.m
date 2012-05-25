function [] = CheckParametersGen(Parameters)

% check if logbeta = beta etc...
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    if not(Parameters.(Names{i}).TransfValue == eval(Parameters.(Names{i}).Transf))
        if not(Parameters.(Names{i}).Value == eval(Parameters.(Names{i}).InvTransf))
            disp('Problem with parameters transfs')
            Names{i}
            die
        end
    end
end
