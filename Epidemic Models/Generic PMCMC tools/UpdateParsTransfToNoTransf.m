function Output = UpdateParsTransfToNoTransf(Parameters,Thetas,Names)


if nargin == 1
    Names = Parameters.Names.Estimated;
    for i = 1:length(Names)
        try
            Parameters.(Names{i}).Value = eval(Parameters.(Names{i}).InvTransf);
        catch
            disp('you')
        end
    end
    CheckParametersGen(Parameters)
    Output = Parameters;
elseif nargin == 3
    Output = zeros(size(Thetas));
    if not(size(Thetas,2) == length(Names))
        disp('Format Problem')
        die
    end
    for i = 1:length(Names)
        if Parameters.(Names{i}).Init
            ind = Parameters.(Names{i}).IndexInit;
        else
            ind = Parameters.(Names{i}).IndexNotInit;
        end
        try
            if strcmp(Parameters.(Names{i}).TransfType,'Log')
                Output(:,i) = exp(Thetas(:,i));
            elseif strcmp(Parameters.(Names{i}).TransfType,'Logit')
                Output(:,i) = InvLogitTransf(Thetas(:,i),Parameters.(Names{i}).MinLim,Parameters.(Names{i}).MaxLim);
            end
        catch
            disp('Problem')
        end
    end
end
