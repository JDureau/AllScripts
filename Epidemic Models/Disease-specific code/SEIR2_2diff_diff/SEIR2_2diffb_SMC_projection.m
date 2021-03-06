function Res = SEIR2_2diffb_SMC_projection(Variables,IndIt,NbIts,Data,Parameters,Path)

    

    gamma = Parameters.gammam1.Value^-1;
    k = Parameters.km1.Value^-1;
    SigmaRW11 = Parameters.SigmaRW11.Value;
    SigmaRW22 = Parameters.SigmaRW22.Value;
    TotPop1 = Parameters.TotalPopulation1;
    TotPop2 = Parameters.TotalPopulation2;
    TempVariables = Variables;
    TempVariables(:,9) = zeros(size(TempVariables(:,9)));
    TempVariables(:,10) = zeros(size(TempVariables(:,10)));
    rands = randn(2,Parameters.NbParticules ,NbIts);

    ComputationTStep = Parameters.ComputationTStep;
    
    for IndDiscr = 1:NbIts
        
        % Variables
        beta11 = exp(TempVariables(:,11));
        beta12 = Parameters.kidsadd.Value + Parameters.kidsmult.Value*exp(TempVariables(:,11));
        beta21 = Parameters.adultsadd.Value + Parameters.adultsmult.Value*exp(TempVariables(:,11));
        beta22 = exp(TempVariables(:,12));
        %S
        TempVariables(:,1) = TempVariables(:,1) + (-beta11.*Variables(:,1).*Variables(:,5)/TotPop1 -beta12.*Variables(:,1).*Variables(:,6)/TotPop2)*ComputationTStep ;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
        TempVariables(:,2) = TempVariables(:,2) + (-beta22.*Variables(:,2).*Variables(:,6)/TotPop2 -beta21.*Variables(:,2).*Variables(:,5)/TotPop1)*ComputationTStep ;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
        %E
        TempVariables(:,3) = TempVariables(:,3) + ( beta11.*Variables(:,1).*Variables(:,5)/TotPop1 + beta12.*Variables(:,1).*Variables(:,6)/TotPop2 -k*Variables(:,3))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands(:,IndDiscr,2);
        TempVariables(:,4) = TempVariables(:,4) + ( beta22.*Variables(:,2).*Variables(:,6)/TotPop2 + beta21.*Variables(:,2).*Variables(:,5)/TotPop1 -k*Variables(:,4))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands(:,IndDiscr,2);
        %I
        TempVariables(:,5) = TempVariables(:,5) + (-gamma*Variables(:,5) + k*Variables(:,3))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
        TempVariables(:,6) = TempVariables(:,6) + (-gamma*Variables(:,6) + k*Variables(:,4))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
        %R
        TempVariables(:,7) = TempVariables(:,7) + ( gamma*Variables(:,5))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
        TempVariables(:,8) = TempVariables(:,8) + ( gamma*Variables(:,6))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
        %Inc
        TempVariables(:,9)  = TempVariables(:,9)  + ( k*Variables(:,3))*ComputationTStep ;
        TempVariables(:,10) = TempVariables(:,10) + ( k*Variables(:,4))*ComputationTStep ;
        
        
        
        TempVariables(:,1) = max(TempVariables(:,1),0);
        TempVariables(:,2) = max(TempVariables(:,2),0);
        TempVariables(:,3) = max(TempVariables(:,3),0);
        TempVariables(:,4) = max(TempVariables(:,4),0);
        TempVariables(:,5) = max(TempVariables(:,5),0);
        TempVariables(:,6) = max(TempVariables(:,6),0);
        TempVariables(:,1) = min(TempVariables(:,1),Parameters.TotalPopulation1);
        TempVariables(:,2) = min(TempVariables(:,2),Parameters.TotalPopulation2);
        TempVariables(:,3) = min(TempVariables(:,3),Parameters.TotalPopulation1);
        TempVariables(:,4) = min(TempVariables(:,4),Parameters.TotalPopulation2);
        TempVariables(:,5) = min(TempVariables(:,5),Parameters.TotalPopulation1);
        TempVariables(:,6) = min(TempVariables(:,6),Parameters.TotalPopulation2);
        
     
        TempVariables(:,11) = TempVariables(:,11) + sqrt(ComputationTStep)*SigmaRW11*squeeze(rands(1,:,IndDiscr))';
        TempVariables(:,12) = TempVariables(:,12) + sqrt(ComputationTStep)*SigmaRW22*squeeze(rands(2,:,IndDiscr))';
        
%         if IndDiscr == 2
%             (cov(TempVariables))
%             die
%         end            

        Variables = TempVariables;
        if not(Parameters.NoPaths)
            Path(:,:,sum(Data.NbComputingSteps(1:IndIt-1)) + IndDiscr) = TempVariables(:,Parameters.PathsToKeep);
        end 
    end
    

    
    Res.Paths = Path;
    Res.Variables = Variables; 