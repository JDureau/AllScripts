function Res = HIVECP_SMC_projection(Variables,IndIt,NbIts,Data,Parameters,Path)



    TempVariables = Variables;
    rands = randn(Parameters.NbParticules ,NbIts);

    TStep = Parameters.ComputationTStep;
    
    WentOutOrNot = ones(Parameters.NbParticules,1);
    
    
    for IndDiscr = 1:NbIts
        % Variables
        
        TempVariables(:,1) = TempVariables(:,1) + ( -Parameters.Beta.Value*Variables(:,1).*Variables(:,2) + Parameters.Gamma.Value*(Variables(:,2)+Variables(:,3)))*TStep;
        TempVariables(:,2) = TempVariables(:,2) + ( Parameters.Beta.Value*Variables(:,1).*Variables(:,2) - (Parameters.Alpha.Value+Parameters.Gamma.Value)*(Variables(:,2)))*TStep;
        TempVariables(:,3) = TempVariables(:,3) + ( Parameters.Alpha.Value*Variables(:,2) - Parameters.Gamma.Value*Variables(:,3))*TStep;
        TempVariables(:,4) = TempVariables(:,2) + TempVariables(:,3);
            
        Variables = TempVariables;
        if not(Parameters.NoPaths)
            Path(:,:,sum(Data.NbComputingSteps(1:IndIt-1)) + IndDiscr) = TempVariables(:,Parameters.PathsToKeep);
        end 

    end
    
    if sum(sum(isnan(Variables)))>0
        disp('There''s a problem...')
        die
    end

    Res.WentOutOrNot = WentOutOrNot;
    Res.Paths = Path;
    Res.Variables = Variables; 