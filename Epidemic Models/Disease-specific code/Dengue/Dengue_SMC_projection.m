 function Res = Dengue_SMC_projection(Variables,IndIt,NbIts,Data,Parameters,Path)



    TempVariables = Variables;
    rands = randn(Parameters.NbParticules ,NbIts);

    TStep = Parameters.ComputationTStep;
    
 
    TempVariables(:,11) = zeros(size(TempVariables(:,11)));
    
    for IndDiscr = 1:NbIts
        % Variables
            

        
        e = Parameters.e.Value;
        d = Parameters.d.Value;
        seas = e*(1+sin(2*pi*((Data.Instants(IndIt)+IndDiscr )*Parameters.ComputationTStep/12 + d)));
        MuB = Parameters.MuB.Value;
        MuD = Parameters.MuD.Value;
        Tot = Parameters.PopulationSize.Value;
        r0  = exp(Variables(:,10));
        v = Parameters.vm1.Value;
        ade = Parameters.ade.Value;
        eta = Parameters.eta.Value;
        q = Parameters.qm1.Value;
        
        TempVariables(:,1) = TempVariables(:,1) + (MuB*Tot - MuD*Variables(:,1) - r0/Tot*v*seas.*(Variables(:,2)+ade*Variables(:,8)+eta).*Variables(:,1) - r0/Tot*v*seas.*(Variables(:,3)+ade*Variables(:,9)+eta).*Variables(:,1))*TStep;
        TempVariables(:,2) = TempVariables(:,2) + (- MuD*Variables(:,2) + r0/Tot*v*seas.*(Variables(:,2)+ade*Variables(:,8)+eta).*Variables(:,1)-Variables(:,2)*v)*TStep;
        TempVariables(:,3) = TempVariables(:,3) + (- MuD*Variables(:,3) + r0/Tot*v*seas.*(Variables(:,3)+ade*Variables(:,9)+eta).*Variables(:,1)-Variables(:,3)*v)*TStep;
        TempVariables(:,4) = TempVariables(:,4) + (- MuD*Variables(:,4) + Variables(:,2)*v - Variables(:,4)*q)*TStep;
        TempVariables(:,5) = TempVariables(:,5) + (- MuD*Variables(:,5) + Variables(:,3)*v - Variables(:,5)*q)*TStep;
        TempVariables(:,6) = TempVariables(:,6) + (- MuD*Variables(:,6) + Variables(:,4)*q  - r0/Tot*v*seas.*(Variables(:,3)+ade*Variables(:,9)+eta).*Variables(:,6))*TStep;
        TempVariables(:,7) = TempVariables(:,7) + (- MuD*Variables(:,7) + Variables(:,5)*q  - r0/Tot*v*seas.*(Variables(:,2)+ade*Variables(:,8)+eta).*Variables(:,7))*TStep;
        TempVariables(:,8) = TempVariables(:,8) + (- MuD*Variables(:,8) + r0/Tot*v*seas.*(Variables(:,3)+ade*Variables(:,9)+eta).*Variables(:,6)-Variables(:,8)*v)*TStep;
        TempVariables(:,9) = TempVariables(:,9) + (- MuD*Variables(:,9) + r0/Tot*v*seas.*(Variables(:,2)+ade*Variables(:,8)+eta).*Variables(:,7)-Variables(:,9)*v)*TStep;
        TempVariables(:,11) = TempVariables(:,11) + (r0/Tot*v*seas.*(Variables(:,3)+ade*Variables(:,9)+eta).*Variables(:,6)+r0/Tot*v*seas.*(Variables(:,2)+ade*Variables(:,8)+eta).*Variables(:,7))*TStep;
               

        
        TempVariables(:,10) = TempVariables(:,10) + sqrt(TStep)*Parameters.SigmaRW.Value*rands(:,IndDiscr);            
        TempVariables(:,1:9) = max(0,TempVariables(:,1:9));
        
        if  (sum(sum(isnan(TempVariables)))>0)
            disp('There''s a problem...')
            die
        end
            
        Variables = TempVariables;
        if not(Parameters.NoPaths)
            Path(:,:,sum(Data.NbComputingSteps(1:IndIt-1)) + IndDiscr) = TempVariables(:,Parameters.PathsToKeep);
        end 

    end
    
    if sum(sum(isnan(Variables)))>0
        disp('There''s a problem...')
        die
    end

    Res.Paths = Path;
    Res.Variables = Variables; 