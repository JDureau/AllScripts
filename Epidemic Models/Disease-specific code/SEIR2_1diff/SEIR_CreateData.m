function Data = SEIR_CreateData(Betas,Parameters,Data,Model)
    

    
    Parameters = Model.InitializeParameters(Parameters);
 
    CheckParametersGen(Parameters)

    Variables = Parameters.InitialState;

    
    gamma = Parameters.gammam1.Value^-1;
    k = Parameters.km1.Value^-1;
    TotPop = Parameters.TotalPopulation;
    TempVariables = Variables;
    

%     gamma = 0.94963^-1;
%     k = 1.5647^-1;
%     TotPop = 100000;
    
    
    TempVariables = Variables;


    ComputationTStep = Parameters.ComputationTStep;
    
    BuiltObs = 0; 
    

    
    IndBeta = 1;
    Record = zeros((length(Data.Instants)-1)*Data.NbComputingSteps(2),7);
    for IndIt = 2:length(Data.Instants)
        NbIts = Data.NbComputingSteps(IndIt);
        Variables (5) = 0;
        TempVariables (5) = 0;
        for IndDiscr = 1:NbIts

            % Variables
            beta = Betas(IndBeta);
%             beta = R0s(IndBeta)*Parameters.gamma.Value*Parameters.TotalPopulation/TempVariables(1);
            IndBeta = IndBeta + 1;
            TempVariables(1) = TempVariables(1) + (-beta.*Variables(1).*Variables(3)/TotPop)*ComputationTStep ;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
            TempVariables(2) = TempVariables(2) + ( beta.*Variables(1).*Variables(3)/TotPop-k*Variables(2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands(:,IndDiscr,2);
            TempVariables(3) = TempVariables(3) + (-gamma*Variables(3) + k*Variables(2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
            TempVariables(4) = TempVariables(4) + ( gamma*Variables(3))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
            TempVariables(5) = TempVariables(5) + ( k*Variables(2))*ComputationTStep ;
            TempVariables(1) = max(TempVariables(1),0);
            TempVariables(2) = max(TempVariables(2),0);
            TempVariables(3) = max(TempVariables(3),0);
            Variables = TempVariables;
            Record(IndBeta,:) = Variables;
            if TempVariables(5) ==0
                disp('obs to zero... problematic!')
            end
        end
        BuiltObs(IndIt) = exp(log(Variables(5))+Parameters.SigmaObs.Value*randn(1,1));
    end    

%     disp(Record(2,1)/Parameters.TotalPopulation)
    Data.Observations = zeros(6,length(BuiltObs));
    Data.Observations(5,:) = BuiltObs;
    Data.RealBetaTraj = Betas;
%     plot(BuiltObs)