function Res = fullvol_SMC_projection(Variables,IndIt,NbIts,Data,Parameters,Path)

    kappa = Parameters.kappa.Value;
    mu_X  = Parameters.mu_X.Value;
    sigma_X  = Parameters.sigma_X.Value;
    mu_Y  = Parameters.mu_Y.Value;
    rho  = Parameters.rho.Value;

    TempVariables = Variables;
    TempVariables(:,5) = zeros(size(TempVariables(:,5)));
    rands = randn(Parameters.NbParticules ,NbIts);

    ComputationTStep = Parameters.ComputationTStep;
    
    for IndDiscr = 1:NbIts
        
        % Variables
        
        %X
        TempVariables(:,1) = TempVariables(:,1) + kappa*(mu_X-Variables(:,1))*ComputationTStep + sigma_X*sqrt(ComputationTStep)*rands(1,IndDiscr);% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
        
        % drift
        TempVariables(:,2) = TempVariables(:,2) + (mu_Y - Vol(Variables(:,1))^2/2)*ComputationTStep + rho*Vol(Variables(:,1))*rands(1,IndDiscr);%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands(:,IndDiscr,2);
        
        % vol
        TempVariables(:,3) = TempVariables(:,3) + (1-rho^2)*Vol(Variables(:,1))^2*ComputationTStep ;

        Variables = TempVariables;
        if not(Parameters.NoPaths)
            Path(:,:,sum(Data.NbComputingSteps(1:IndIt-1)) + IndDiscr) = TempVariables(:,Parameters.PathsToKeep);
        end 
    end
    

    
    Res.Paths = Path;
    Res.Variables = Variables; 