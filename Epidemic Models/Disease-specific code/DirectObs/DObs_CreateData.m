function Data = DObs_CreateData(Betas,Parameters,Data,Model)
    
    Parameters = Model.InitializeParameters(Parameters);
 
    CheckParametersGen(Parameters)

    Variables = Parameters.InitialState;


    TempVariables = Variables;


    ComputationTStep = Parameters.ComputationTStep;
    
    BuiltObs = 0; 
    

    
    IndBeta = 2;
    Record = [];
    for IndIt = 2:length(Data.Instants)
        NbIts = Data.NbComputingSteps(IndIt);

        for IndDiscr = 1:NbIts

            % Variables
            beta = Betas(IndBeta);
%             beta = R0s(IndBeta)*Parameters.gamma.Value*Parameters.TotalPopulation/TempVariables(1);
            IndBeta = IndBeta + 1;
            TempVariables(1) = beta ;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
            Variables = TempVariables;
            Record(IndBeta,:) = Variables;

        end
        BuiltObs(IndIt) = Variables(1)*1+Parameters.ObsNoise*randn(1,1);
    end    

%     disp(Record(2,1)/Parameters.TotalPopulation)
    Data.Observations = zeros(1,length(BuiltObs));
    Data.Observations(1,:) = BuiltObs;
    Data.RealBetaTraj = Betas;
%     plot(BuiltObs)