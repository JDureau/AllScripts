function Res = SEIRPlusObs_SMC_projection(Variables,IndIt,NbIts,Data,Parameters,Path)

    

    gamma = Parameters.gammam1.Value^-1;
    k = Parameters.km1.Value^-1;
    SigmaLambda = Parameters.SigmaRW.Value;
    SigmaObs = Parameters.SigmaRWObs.Value;
    TotPop = Parameters.TotalPopulation;
    TempVariables = Variables;
    TempVariables(:,5) = zeros(size(TempVariables(:,5)));
    rands = randn(Parameters.NbParticules ,NbIts,2);

    ComputationTStep = Parameters.ComputationTStep;
    
    for IndDiscr = 1:NbIts
        
        % Variables
        beta = exp(Variables(:,6));
        PropObs = exp(Variables(:,7));
        TempVariables(:,1) = TempVariables(:,1) + (-beta.*Variables(:,1).*Variables(:,3)/TotPop)*ComputationTStep ;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
        TempVariables(:,2) = TempVariables(:,2) + ( beta.*Variables(:,1).*Variables(:,3)/TotPop-k*Variables(:,2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands(:,IndDiscr,2);
        TempVariables(:,3) = TempVariables(:,3) + (-gamma*Variables(:,3) + k*Variables(:,2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
        TempVariables(:,4) = TempVariables(:,4) + ( gamma*Variables(:,3))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
        TempVariables(:,5) = TempVariables(:,5) + ( k*Variables(:,2).*PropObs)*ComputationTStep ;
        TempVariables(:,1) = max(TempVariables(:,1),0);
        TempVariables(:,2) = max(TempVariables(:,2),0);
        TempVariables(:,3) = max(TempVariables(:,3),0);
        if Parameters.DiffusionType =='OUD'
            TempVariables(:,6) = TempVariables(:,6) + (Mu-Variables(:,6))*Kappa*ComputationTStep + sqrt(ComputationTStep)*SigmaLambda*rands(:,IndDiscr);            
        elseif Parameters.DiffusionType =='Add'
            TempVariables(:,6) = TempVariables(:,6) + sqrt(ComputationTStep)*SigmaLambda*rands(:,IndDiscr,1);
            TempVariables(:,7) = min(0,TempVariables(:,7) + sqrt(ComputationTStep)*SigmaObs*rands(:,IndDiscr,2));
        else
            disp('Unknown diffusion type')
            die
        end
        
            
        Variables = TempVariables;
        if not(Parameters.NoPaths)
            Path(:,:,sum(Data.NbComputingSteps(1:IndIt-1)) + IndDiscr) = TempVariables(:,Parameters.PathsToKeep);
        end 
    end
    

    
    Res.Paths = Path;
    Res.Variables = Variables; 