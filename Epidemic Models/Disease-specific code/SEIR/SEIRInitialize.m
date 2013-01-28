function Parameters = fullvolInitialize(Parameters)


 
try
    n = Parameters.NbParticules;
catch
    n = 1000;
%     disp('NbParticules temporarily set to 1000')
end


InitialState(1) = Parameters.X0.Value;
InitialState(2) = 0;
InitialState(3) = 0;

InitialCov = zeros(Parameters.NbVariables,Parameters.NbVariables);

InitialStates(1,:) = Parameters.X0.Value*ones(1,n);
InitialStates(2,:) = zeros(1,n);
InitialStates(3,:) = zeros(1,n);

% InitialStates(2,:) = min(Parameters.TotalPopulation,InitialStates(2,:));
% InitialStates(3,:) = min(Parameters.TotalPopulation-InitialStates(2,:),InitialStates(3,:));
% InitialStates(4,:) = min(Parameters.TotalPopulation-InitialStates(2,:)-InitialStates(3,:),InitialStates(4,:));


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