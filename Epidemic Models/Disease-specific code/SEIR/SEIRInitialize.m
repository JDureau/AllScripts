function Parameters = SEIRInitialize(Parameters)






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
%     disp('NbParticules temporarily set to 1000')
end


InitialState(2) = Parameters.TotalPopulation*(Parameters.EInitProp.Value);
InitialState(3) = Parameters.TotalPopulation*(Parameters.IInitProp.Value);
InitialState(4) = Parameters.TotalPopulation*(Parameters.RInitProp.Value);
InitialState(1) = max(0,Parameters.TotalPopulation -InitialState(2)-InitialState(3)-InitialState(4));
InitialState(5) = 0;
InitialState(6) = log(max(0,Parameters.betainit.Value));
InitialState(7) = 0;

InitialCov = zeros(Parameters.NbVariables,Parameters.NbVariables);


InitialStates(2,:) = Parameters.TotalPopulation*((ones(1,n))*Parameters.EInitProp.Value);
InitialStates(3,:) = Parameters.TotalPopulation*((ones(1,n))*Parameters.IInitProp.Value);
InitialStates(4,:) = Parameters.TotalPopulation*((ones(1,n))*Parameters.RInitProp.Value);
InitialStates(1,:) = max(0,Parameters.TotalPopulation*ones(1,n) -InitialStates(2,:)-InitialStates(3,:)-InitialStates(4,:));
InitialStates(5,:) = zeros(1,n);
InitialStates(6,:) = log(max(0,(ones(1,n))*Parameters.betainit.Value));
InitialStates(7,:) = zeros(1,n);

InitialStates(2,:) = min(Parameters.TotalPopulation,InitialStates(2,:));
InitialStates(3,:) = min(Parameters.TotalPopulation-InitialStates(2,:),InitialStates(3,:));
InitialStates(4,:) = min(Parameters.TotalPopulation-InitialStates(2,:)-InitialStates(3,:),InitialStates(4,:));
InitialStates(2,:) = max(InitialStates(2,:),0);
InitialStates(3,:) = max(InitialStates(3,:),0);
InitialStates(4,:) = max(InitialStates(4,:),0);

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