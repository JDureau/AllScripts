function Res = fullvol_SMC_projection(Variables,IndIt,NbIts,Data,Parameters,Path)

    kappa = Parameters.kappa.Value;
    mu_X  = Parameters.mu_X.Value;
    sigma_X  = Parameters.sigma_X.Value;
    mu_Y  = Parameters.mu_Y.Value;
    rho  = Parameters.rho.Value;
    Vol = Parameters.Vol;

    TempVariables = Variables;
    TempVariables(:,2) = zeros(size(TempVariables(:,2)));
    TempVariables(:,3) = zeros(size(TempVariables(:,3)));
    rands = randn(Parameters.NbParticules ,NbIts);

    ComputationTStep = Parameters.ComputationTStep;
    
    for IndDiscr = 1:NbIts
        
        % Variables
        vols = Vol(Variables(:,1));
        
        %X
        TempVariables(:,1) = TempVariables(:,1) + kappa*(mu_X-Variables(:,1))*ComputationTStep + sigma_X*sqrt(ComputationTStep)*rands(:,IndDiscr);% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
        
        % drift
        TempVariables(:,2) = TempVariables(:,2) + (mu_Y - vols.^2/2)*ComputationTStep + rho*vols*sqrt(ComputationTStep).*rands(:,IndDiscr);%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands(:,IndDiscr,2);
        
        % vol
        TempVariables(:,3) = TempVariables(:,3) + (1-rho^2)*vols.^2*ComputationTStep ;

        Variables = TempVariables;
        if not(Parameters.NoPaths)
            Path(:,:,sum(Data.NbComputingSteps(1:IndIt-1)) + IndDiscr) = TempVariables(:,Parameters.PathsToKeep);
        end 
    end
    

    
    Res.Paths = Path;
    Res.Variables = Variables; 