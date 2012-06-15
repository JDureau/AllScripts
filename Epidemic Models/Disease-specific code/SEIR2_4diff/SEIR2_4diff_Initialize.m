function Parameters = SEIR2_3diff_Initialize(Parameters)






Parameters.E1InitProp.Value = min(1,Parameters.E1InitProp.Value);
Parameters.I1InitProp.Value = min(1-Parameters.E1InitProp.Value,Parameters.I1InitProp.Value);
Parameters.R1InitProp.Value = min(1-Parameters.E1InitProp.Value-Parameters.I1InitProp.Value,Parameters.R1InitProp.Value);

Parameters.E1InitProp.Value = max(0,Parameters.E1InitProp.Value);
Parameters.I1InitProp.Value = max(0,Parameters.I1InitProp.Value);
Parameters.R1InitProp.Value = max(0,Parameters.R1InitProp.Value);

Parameters.E2InitProp.Value = min(2,Parameters.E2InitProp.Value);
Parameters.I2InitProp.Value = min(2-Parameters.E2InitProp.Value,Parameters.I2InitProp.Value);
Parameters.R2InitProp.Value = min(2-Parameters.E2InitProp.Value-Parameters.I2InitProp.Value,Parameters.R2InitProp.Value);

Parameters.E2InitProp.Value = max(0,Parameters.E2InitProp.Value);
Parameters.I2InitProp.Value = max(0,Parameters.I2InitProp.Value);
Parameters.R2InitProp.Value = max(0,Parameters.R2InitProp.Value);

Parameters = UpdateParsNoTransfToTransf(Parameters);

try
    n = Parameters.NbParticules;
catch
    n = 1000;
%     disp('NbParticules temporarily set to 1000')
end


InitialState(3) = Parameters.TotalPopulation1*(Parameters.E1InitProp.Value);
InitialState(4) = Parameters.TotalPopulation2*(Parameters.E2InitProp.Value);
InitialState(5) = Parameters.TotalPopulation1*(Parameters.I1InitProp.Value);
InitialState(6) = Parameters.TotalPopulation2*(Parameters.I2InitProp.Value);
InitialState(7) = Parameters.TotalPopulation1*(Parameters.R1InitProp.Value);
InitialState(8) = Parameters.TotalPopulation2*(Parameters.R2InitProp.Value);
InitialState(1) = max(0,Parameters.TotalPopulation1 -InitialState(3)-InitialState(5)-InitialState(7));
InitialState(2) = max(0,Parameters.TotalPopulation2 -InitialState(4)-InitialState(6)-InitialState(8));
InitialState(9) = 0;
InitialState(10) = 0;
InitialState(11) = log(max(0,Parameters.beta11init.Value));
InitialState(12) = log(max(0,Parameters.beta22init.Value));
InitialState(13) = log(max(0,Parameters.beta12init.Value));
InitialState(14) = log(max(0,Parameters.beta21init.Value));



InitialCov = zeros(Parameters.NbVariables,Parameters.NbVariables);

InitialStates = repmat(InitialState,n,1)';


% try
%     if strcmp(Parameters.DiffusionType,'IBM')
%         InitialStates(7,:) = Parameters.betaderinit.Value*(ones(1,n));
%         InitialState(7) = Parameters.betaderinit.Value;
%     elseif strcmp(Parameters.DiffusionType,'SVO')
%         InitialStates(7,:) = Parameters.VolInit.Value;
%     else
%         InitialStates(7,:) = zeros(1,n);
%         InitialState(7) = 0;
%     end
% catch
%     InitialStates(7,:) = zeros(1,n);
%     InitialState(7) = 0;
% end
Parameters.InitialStates = InitialStates;
Parameters.InitialState = InitialState';
Parameters.InitialCov = InitialCov;