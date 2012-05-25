function Res = SIR_SMC_projection(Variables,IndIt,NbIts,Data,Parameters,Path)

    Beta = Parameters.Beta.Value;
    Gamma = Parameters.Gamma.Value;
    SigmaLambda = Parameters.SigmaRW.Value;

    LambdaNm1 = Variables(:,4);
    TempVariables = Variables;
    rands = randn(Parameters.NbParticules ,NbIts);

    ComputationTStep = Parameters.ComputationTStep;
    
    for IndDiscr = 1:NbIts
        
        % Variables
        TempVariables(:,1) = TempVariables(:,1) + (-Variables(:,4).*Variables(:,1).*Variables(:,2)+Beta*Variables(:,3))*ComputationTStep ;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
        TempVariables(:,2) = TempVariables(:,2) + ( Variables(:,4).*Variables(:,1).*Variables(:,2)-Gamma*Variables(:,2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands(:,IndDiscr,2);
        TempVariables(:,3) = TempVariables(:,3) + (-Beta*Variables(:,3) + Gamma*Variables(:,2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
        if Parameters.DiffusionType =='OUD'
            TempVariables(:,4) = TempVariables(:,4) + (Mu-Variables(:,4))*Kappa*ComputationTStep + sqrt(ComputationTStep)*SigmaLambda*rands(:,IndDiscr);            
        elseif Parameters.DiffusionType =='Add'
            TempVariables(:,4) = TempVariables(:,4) + sqrt(ComputationTStep)*SigmaLambda*rands(:,IndDiscr);
        else
            disp('Unknown diffusion type')
            die
        end
        
            
        Variables = TempVariables;
        if not(Parameters.NoPaths)
            Path(:,:,sum(Data.NbComputingSteps(1:IndIt-1)) + IndDiscr) = TempVariables;
        end 
    end
    

    
    Res.Paths = Path;
    Res.Variables = Variables; 