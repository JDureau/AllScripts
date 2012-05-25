function Res = DObs_SMC_projection(Variables,IndIt,NbIts,Data,Parameters,Path)

    SigmaLambda = Parameters.SigmaRW.Value;
    TempVariables = Variables;

    rands = randn(Parameters.NbParticules ,NbIts);

    ComputationTStep = Parameters.ComputationTStep;
    
    for IndDiscr = 1:NbIts
        TempVariables(:,1) = TempVariables(:,1) + sqrt(ComputationTStep)*SigmaLambda*rands(:,IndDiscr);
        
        try
            % This is for GIBBS PMCMC
             if Parameters.ForceTraj
                 TempVariables(1,1) = Parameters.ForcedTraj(1,sum(Data.NbComputingSteps(1:IndIt-1)) + IndDiscr);
            end
        end
        Variables = TempVariables;
        if not(Parameters.NoPaths)
            Path(:,:,sum(Data.NbComputingSteps(1:IndIt-1)) + IndDiscr ) = TempVariables(:,Parameters.PathsToKeep);
        end 
    end
    

    
    Res.Paths = Path;
    Res.Variables = Variables; 