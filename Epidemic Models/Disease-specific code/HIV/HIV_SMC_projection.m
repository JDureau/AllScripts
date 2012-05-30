 function Res = HIV_SMC_projection(Variables,IndIt,NbIts,Data,Parameters,Path)



    TempVariables = Variables;
    rands = randn(Parameters.NbParticules ,NbIts);

    TStep = Parameters.ComputationTStep;
    
    WentOutOrNot = ones(Parameters.NbParticules,1);
    
    if or(strcmp(Parameters.DiffusionType,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
        m = Parameters.BRmm1.Value + 1;
        mu = Parameters.BRmu.Value;
        k = Parameters.k;
    elseif strcmp(Parameters.DiffusionType,'Sigmoid')
        % beta(t) = a + b/c(1+x_t)
        rate = Parameters.Sigmrate.Value;
        base = Parameters.Sigmbase.Value;
        mu = Parameters.Sigmmu.Value;
        tinfl = Parameters.Sigmtinfl.Value;
        
        c = 1/(1+exp(tinfl/rate));
        b = (mu-base)*c/(1-c);
        a = base - b;
    end
    
    Crash = 0;
    if Crash 
        NbIts = 1;
    end
    for IndDiscr = 1:NbIts
        % Variables
        if or(strcmp(Parameters.DiffusionType,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
            beta = ((1-m)*Variables(:,9)+mu^(1-m)).^(1/(1-m));
            if Parameters.BRbase.Value>m^(1/(1-m))*mu
                Crash = 1;
                beta = Parameters.BRbase.Value;
            end
        elseif strcmp(Parameters.DiffusionType,'Sigmoid')
            Variables(:,9) = max(0.000001,Variables(:,9));
            beta = a + b./(c*(1+min(10^6,Variables(:,9))));
        elseif or(strcmp(Parameters.DiffusionType,'Add'),strcmp(Parameters.DiffusionType,'AddConstr'))
            beta = min(1,max(0,exp(Variables(:,9))./(1+exp(Variables(:,9)))));
        end    
        TotF1 = Parameters.TotF1.Value;
        TotF2 = Parameters.TotF2.Value;
        TotM  = Parameters.TotM.Value;
        TempVariables(:,1) = TempVariables(:,1) + ( Parameters.MuFm1.Value^-1*(TotF1) + Parameters.Alpham1.Value^-1*Variables(:,2) - Parameters.BetaMF.Value*Parameters.CF1.Value*(1-Parameters.eHIV.Value*beta).*Variables(:,1).*Variables(:,6)./TotM-Parameters.MuFm1.Value^-1*Variables(:,1))*TStep;
        TempVariables(:,2) = TempVariables(:,2) + ( Parameters.BetaMF.Value*Parameters.CF1.Value*(1-Parameters.eHIV.Value*beta).*Variables(:,1).*Variables(:,6)./TotM - (Parameters.MuFm1.Value^-1 + Parameters.Alpham1.Value^-1)*Variables(:,2))*TStep;
        TempVariables(:,3) = TempVariables(:,3) + ( Parameters.MuFm1.Value^-1*(TotF2) + Parameters.Alpham1.Value^-1*Variables(:,4) - Parameters.BetaMF.Value*Parameters.CF2.Value*(1-Parameters.eHIV.Value*beta).*Variables(:,3).*Variables(:,6)./TotM-Parameters.MuFm1.Value^-1*Variables(:,3))*TStep;
        TempVariables(:,4) = TempVariables(:,4) + ( Parameters.BetaMF.Value*Parameters.CF2.Value*(1-Parameters.eHIV.Value*beta).*Variables(:,3).*Variables(:,6)./TotM - (Parameters.MuFm1.Value^-1 + Parameters.Alpham1.Value^-1)*Variables(:,4))*TStep;
        TempVariables(:,5) = TempVariables(:,5) + ( Parameters.MuMm1.Value^-1*(TotM) + Parameters.Alpham1.Value^-1*Variables(:,6) - Parameters.BetaFM.Value*Parameters.NbContactsForMen*(1-Parameters.eHIV.Value*beta).*Variables(:,5).*(Parameters.CF1.Value*TotF1./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(:,2)./TotF1+Parameters.CF2.Value*TotF2./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(:,4)./TotF2)-Parameters.MuMm1.Value^-1*Variables(:,5))*TStep;
        TempVariables(:,6) = TempVariables(:,6) + ( Parameters.BetaFM.Value*Parameters.NbContactsForMen*(1-Parameters.eHIV.Value*beta).*Variables(:,5).*(Parameters.CF1.Value*TotF1./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(:,2)./TotF1+Parameters.CF2.Value*TotF2./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(:,4)./TotF2) - (Parameters.MuMm1.Value^-1+Parameters.Alpham1.Value^-1)*Variables(:,6))*TStep;
        TempVariables(:,7) = (TempVariables(:,2) + TempVariables(:,4))/(TotF1+TotF2)*100; 
        TempVariables(:,8) = (TempVariables(:,6))/(TotM)*100;         

%         if or(not(Parameters.StableCUseConstraint),sum(Data.NbComputingSteps(1:IndIt-1)) + IndDiscr>204/Parameters.ComputationTStep)
%             if Parameters.DiffusionType =='OUD'
%                 if IndDiscr<204/Parameters.ComputationTStep
%                     Mu = Parameters.MuOU1.TransfValue;
%                 else
%                     Mu = Parameters.MuOU2.TransfValue;
%                 end
%                 TransfFt = TransfFt + (Mu-TransfFt)*Parameters.KappaOU.Value*TStep + sqrt(TStep)*Parameters.SigmaOU.Value*rands(:,IndDiscr);            
%             elseif Parameters.DiffusionType =='Add'
%                 TransfFt = TransfFt + sqrt(TStep)*Parameters.SigmaRW.Value*rands(:,IndDiscr);
%             else
%                 disp('Unknown diffusion type')
%                 die
%             end
%             try
%                 TempVariables(:,9) = InvLogitTransf(TransfFt,0,1);
%             catch
%                 die
%             end
%         end
        if 1
            if or(strcmp(Parameters.DiffusionType,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
                try
                    TempVariables(:,9) = TempVariables(:,9) - k*TempVariables(:,9)*TStep + sqrt(TStep)*Parameters.BRsigma.Value*TempVariables(:,9).*rands(:,IndDiscr);            
                    WentOutOrNot = WentOutOrNot.*(TempVariables(:,9)<0);
                    tmpinds = find(not(WentOutOrNot));
                    TempVariables(tmpinds,9) = -0.00001;
                catch
                    'problem'
                end
            elseif strcmp(Parameters.DiffusionType,'Sigmoid')
                try
                    TempVariables(:,9) = TempVariables(:,9) - 1/rate*TempVariables(:,9)*TStep + sqrt(TStep)*Parameters.Sigmsigma.Value*TempVariables(:,9).*rands(:,IndDiscr);            
                    WentOutOrNot = WentOutOrNot.*(TempVariables(:,9)>0);
                    tmpinds = find(not(WentOutOrNot));
                    TempVariables(tmpinds,9) = 0.00001;
                catch
                    'problem'
                end
                    %                 TempVariables(:,9) = TempVariables(:,9) + Parameters.CUsteepness.Value*Parameters.k*TStep + sqrt(TStep)*Parameters.SigmaRW.Value/(mu*exp(Variables(:,9))./(1+exp(Variables(:,9))^2))*rands(:,IndDiscr);            
            elseif or(strcmp(Parameters.DiffusionType,'Add'),strcmp(Parameters.DiffusionType,'AddConstr'))
                TempVariables(:,9) = TempVariables(:,9) + sqrt(TStep)*Parameters.SigmaRW.Value*rands(:,IndDiscr);            
            
            else
                disp('Unknown diffusion type')
                die
            end


        end
        TempVariables(:,1:8) = max(0,TempVariables(:,1:8));
        
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

    Res.WentOutOrNot = WentOutOrNot;
    Res.Paths = Path;
    Res.Variables = Variables; 
    Res.Crash = Crash;