function Parameters = SEIR2_1diff_Initialize(Parameters)






Parameters.EInitProp.Value = min(1,Parameters.EInitProp.Value);
Parameters.IInitProp.Value = min(1-Parameters.EInitProp.Value,Parameters.IInitProp.Value);
Parameters.RInitProp.Value = min(1-Parameters.EInitProp.Value-Parameters.IInitProp.Value,Parameters.RInitProp.Value);

Parameters.EInitProp.Value = max(0,Parameters.EInitProp.Value);
Parameters.IInitProp.Value = max(0,Parameters.IInitProp.Value);
Parameters.RInitProp.Value = max(0,Parameters.RInitProp.Value);

Parameters = UpdateParsNoTransfToTransf(Parameters);

try
    n = Parameters.NbParticules;
catch
    n = 1000;
    disp('NbParticules temporarily set to 1000')
end


InitialState(3) = Parameters.TotalPopulation1*(Parameters.EInitProp.Value);
InitialState(4) = Parameters.TotalPopulation2*(Parameters.EInitProp.Value);
InitialState(5) = Parameters.TotalPopulation1*(Parameters.IInitProp.Value);
InitialState(6) = Parameters.TotalPopulation2*(Parameters.IInitProp.Value);
InitialState(7) = Parameters.TotalPopulation1*(Parameters.RInitProp.Value);
InitialState(8) = Parameters.TotalPopulation2*(Parameters.RInitProp.Value);
InitialState(1) = max(0,Parameters.TotalPopulation1 -InitialState(3)-InitialState(5)-InitialState(7));
InitialState(2) = max(0,Parameters.TotalPopulation2 -InitialState(4)-InitialState(6)-InitialState(8));
InitialState(9) = 0;
InitialState(10) = 0;
InitialState(11) = log(max(0,Parameters.beta11init.Value));
InitialState(12) = 0;



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