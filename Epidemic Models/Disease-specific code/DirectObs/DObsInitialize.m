function Parameters = DObsInitialize(Parameters)


try
    n = Parameters.NbParticules;
catch
    n = 1000;
    disp('NbParticules temporarily set to 1000')
end

InitialState(1:2,1) = Parameters.betainit.Value;

InitialStates = repmat(InitialState,1,n);
Parameters.InitialStates = InitialStates;
Parameters.InitialState = InitialState;
Parameters.NbVariables = length(InitialState);