function Res = HIV_EKF_projection(Data,Model,m,Cov,NbIts,IndTime,Parameters)

    mpred = m;
    TStep = Parameters.ComputationTStep;

    
    deltabetas = zeros(1,NbIts);
    record = zeros(9,NbIts);
    
    if or(strcmp(Parameters.DiffusionType,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
        m = Parameters.BRmm1.Value + 1;
        mu = Parameters.BRmu.Value;
        k = Parameters.k;
    elseif strcmp(Parameters.DiffusionType,'Sigmoid')
        rate = Parameters.Sigmrate.Value;
        base = Parameters.Sigmbase.Value;
        mu = Parameters.Sigmmu.Value;
    else
%         disp('Unknown Diffusion')
%         die
    end
        
    MuFm1 = Parameters.MuFm1.Value;
    Alpham1 = Parameters.Alpham1.Value;
    BetaMF = Parameters.BetaMF.Value;
    CF1 = Parameters.CF1.Value;
    BetaFM = Parameters.BetaFM.Value;
    CF2 = Parameters.CF2.Value;
    eHIV = Parameters.eHIV.Value;
    NbContactsForMen = Parameters.NbContactsForMen;
    MuMm1 = Parameters.MuMm1.Value;
    Crash = 0;
    difftype = Parameters.DiffusionType;
    for IndDiscr = 1:NbIts
        
%         if sum(mpred(1:8)<0)
%             'stop'
%         end
        
%         mpred(1:8) = max(0,mpred(1:8));
        
       
        mtemp = mpred;
        if or(strcmp(difftype,'Bertallanfy'),strcmp(difftype,'BertallanfyConstr'))
            mtemp(9) = max(-10^6,min(-0.01,mtemp(9)));
            beta = ((1-m)*mtemp(9)+mu^(1-m))^(1/(1-m));
            if or(Parameters.BRbase.Value>m^(1/(1-m))*mu,not(isreal(k)))
                Crash = 1;
                'BR forbidden beta0'
                k = real(k);
                beta = Parameters.BRbase.Value;
            end
        elseif strcmp(difftype,'Sigmoid')
            mtemp(9) = min(10^6,max(0.000001,mtemp(9)));
            beta = base + (mu-base)/(1+mtemp(9));
        else
            beta = exp(mtemp(9))/(1+exp(mtemp(9)));
        end
        TotF1 = mtemp(1) + mtemp(2);
        TotF2 = mtemp(3) + mtemp(4);
        TotM  = mtemp(5) + mtemp(6);
        mpred(1) = mpred(1) + ( MuFm1^-1*(TotF1) + Alpham1^-1*mtemp(2) - BetaMF*CF1*(1-eHIV*beta)*mtemp(1)*mtemp(6)/TotM-MuFm1^-1*mtemp(1))*TStep;
        mpred(2) = mpred(2) + ( BetaMF*CF1*(1-eHIV*beta)*mtemp(1)*mtemp(6)/TotM - (MuFm1^-1 + Alpham1^-1)*mtemp(2))*TStep;
        mpred(3) = mpred(3) + ( MuFm1^-1*(TotF2) + Alpham1^-1*mtemp(4) - BetaMF*CF2*(1-eHIV*beta)*mtemp(3)*mtemp(6)/TotM-MuFm1^-1*mtemp(3))*TStep;
        mpred(4) = mpred(4) + ( BetaMF*CF2*(1-eHIV*beta)*mtemp(3)*mtemp(6)/TotM - (MuFm1^-1 + Alpham1^-1)*mtemp(4))*TStep;
        mpred(5) = mpred(5) + ( MuMm1^-1*(TotM) + Alpham1^-1*mtemp(6) - BetaFM*NbContactsForMen*(1-eHIV*beta)*mtemp(5)*(CF1*TotF1/(CF2*TotF2+CF1*TotF1)*mtemp(2)./TotF1+CF2*TotF2/(CF2*TotF2+CF1*TotF1)*mtemp(4)/TotF2)-MuMm1^-1*mtemp(5))*TStep;
        mpred(6) = mpred(6) + ( BetaFM*NbContactsForMen*(1-eHIV*beta)*mtemp(5)*(CF1*TotF1/(CF2*TotF2+CF1*TotF1)*mtemp(2)/TotF1+CF2*TotF2/(CF2*TotF2+CF1*TotF1)*mtemp(4)/TotF2) - (MuMm1^-1+Alpham1^-1)*mtemp(6))*TStep;
        mpred(7) = ((mpred(2) + mpred(4))/(TotF1+TotF2))*100; 
        mpred(8) = ((mpred(6))/(TotM))*100; 
        
     
        if or(strcmp(difftype,'Bertallanfy'),strcmp(difftype,'BertallanfyConstr'))
            mpred(9) = mtemp(9) - k*mtemp(9)*TStep;
            if Crash
                mpred(9) = mtemp(9);
            end
        elseif strcmp(difftype,'Sigmoid')
            mpred(9) = mtemp(9) - 1/rate*mtemp(9)*TStep;
        elseif or(strcmp(difftype,'Add'),strcmp(difftype,'AddConstr'))
            mpred(9) = mtemp(9);
        else
            disp('Unknown Diffusion')
            die
        end
        
        if not(isreal(mpred))
            '....'
        end

        
        if or(strcmp(difftype,'Bertallanfy'),strcmp(difftype,'BertallanfyConstr'))
            betader = ((1-m)*mpred(9)+mu^(1-m))^(m/(1-m));
            if Crash
                betader = 0;
            end
        elseif strcmp(difftype,'Sigmoid')
            betader = -(mu-base)/((1+mpred(9))^2);
        else
            betader = exp(mpred(9))/((1+exp(mpred(9)))^2);
        end
        Jacobian = zeros(9,9);
        Jacobian(1,1) =   - BetaMF*CF1*(1-eHIV*beta)*mtemp(6)/TotM-MuFm1^-1;
        Jacobian(1,2) =     Alpham1^-1 ;
        Jacobian(1,3) = 0;
        Jacobian(1,4) = 0;
        Jacobian(1,5) = 0;
        Jacobian(1,6) =   - BetaMF*CF1*(1-eHIV*beta)*mtemp(1)/TotM;
        Jacobian(1,7) = 0;
        Jacobian(1,8) = 0;
        Jacobian(1,9) =   - BetaMF*CF1*(-eHIV*betader)*mtemp(1)*mtemp(6)./TotM;
        Jacobian(2,1) =     BetaMF*CF1*(1-eHIV*beta)*mtemp(6)./TotM ;
        Jacobian(2,2) =   - (MuFm1^-1 + Alpham1^-1);
        Jacobian(2,3) = 0;
        Jacobian(2,4) = 0;
        Jacobian(2,5) = 0;
        Jacobian(2,6) =     BetaMF*CF1*(1-eHIV*beta)*mtemp(1)./TotM ;
        Jacobian(2,7) = 0;
        Jacobian(2,8) = 0;
        Jacobian(2,9) =     BetaMF*CF1*(-eHIV*betader)*mtemp(1)*mtemp(6)./TotM ;
        Jacobian(3,1) = 0;
        Jacobian(3,2) = 0;
        Jacobian(3,3) =  - BetaMF*CF2*(1-eHIV*beta)*mtemp(6)./TotM-MuFm1^-1;
        Jacobian(3,4) =    Alpham1^-1;
        Jacobian(3,5) = 0;
        Jacobian(3,6) =  - BetaMF*CF2*(1-eHIV*beta)*mtemp(3)./TotM;
        Jacobian(3,7) = 0;
        Jacobian(3,8) = 0;
        Jacobian(3,9) =  - BetaMF*CF2*(-eHIV*betader)*mtemp(3)*mtemp(6)./TotM;
        Jacobian(4,1) =  0;
        Jacobian(4,2) =  0;
        Jacobian(4,3) =  BetaMF*CF2*(1-eHIV*beta)*mtemp(6)./TotM ;
        Jacobian(4,4) =  - (MuFm1^-1 + Alpham1^-1);
        Jacobian(4,5) = 0;
        Jacobian(4,6) =   BetaMF*CF2*(1-eHIV*beta)*mtemp(3)./TotM ;
        Jacobian(4,7) = 0;
        Jacobian(4,8) = 0;
        Jacobian(4,9) =    BetaMF*CF2*(-eHIV*betader)*mtemp(3)*mtemp(6)./TotM ;        
        Jacobian(5,1) = 0;
        Jacobian(5,2) =  - BetaFM*NbContactsForMen*(1-eHIV*beta)*mtemp(5)*(CF1*TotF1./(CF2*TotF2+CF1*TotF1)./TotF1);
        Jacobian(5,3) = 0;
        Jacobian(5,4) =  - BetaFM*NbContactsForMen*(1-eHIV*beta)*mtemp(5)*(CF2*TotF2./(CF2*TotF2+CF1*TotF1)./TotF2);
        Jacobian(5,5) =  - BetaFM*NbContactsForMen*(1-eHIV*beta)*(CF1*TotF1./(CF2*TotF2+CF1*TotF1)*mtemp(2)./TotF1+CF2*TotF2./(CF2*TotF2+CF1*TotF1)*mtemp(4)./TotF2)-MuMm1^-1;
        Jacobian(5,6) =    Alpham1^-1;
        Jacobian(5,7) = 0;
        Jacobian(5,8) = 0;
        Jacobian(5,9) =  - BetaFM*NbContactsForMen*(-eHIV*betader)*mtemp(5)*(CF1*TotF1./(CF2*TotF2+CF1*TotF1)*mtemp(2)./TotF1+CF2*TotF2./(CF2*TotF2+CF1*TotF1)*mtemp(4)./TotF2);
        Jacobian(6,1) = 0;
        Jacobian(6,2) =  BetaFM*NbContactsForMen*(1-eHIV*beta)*mtemp(5)*(CF1*TotF1./(CF2*TotF2+CF1*TotF1)./TotF1) ;
        Jacobian(6,3) = 0;
        Jacobian(6,4) =  BetaFM*NbContactsForMen*(1-eHIV*beta)*mtemp(5)*(CF2*TotF2./(CF2*TotF2+CF1*TotF1)./TotF2) ;
        Jacobian(6,5) =  BetaFM*NbContactsForMen*(1-eHIV*beta)*(CF1*TotF1./(CF2*TotF2+CF1*TotF1)*mtemp(2)./TotF1+CF2*TotF2./(CF2*TotF2+CF1*TotF1)*mtemp(4)./TotF2) ;
        Jacobian(6,6) =  - (MuMm1^-1+Alpham1^-1);
        Jacobian(6,7) = 0;
        Jacobian(6,8) = 0;
        Jacobian(6,9) =  BetaFM*NbContactsForMen*(-eHIV*betader)*mtemp(5)*(CF1*TotF1./(CF2*TotF2+CF1*TotF1)*mtemp(2)./TotF1+CF2*TotF2./(CF2*TotF2+CF1*TotF1)*mtemp(4)./TotF2);
        Jacobian(7,:) = ((Jacobian(2,:)+Jacobian(4,:))/(TotF1+TotF2))*100;
        Jacobian(8,:) = (Jacobian(6,:)/(TotM))*100;

        
        if strcmp(difftype,'OUD')
            Jacobian(9,9) = -Parameters.KappaOU.Value;
        elseif strcmp(difftype,'IBM')
            Jacobian(9,10) = 1;
        elseif or(strcmp(difftype,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
            Jacobian(9,9) = - k;
        elseif strcmp(difftype,'Sigmoid')
            Jacobian(9,9) = - 1/rate;
        end
        
        Q = zeros(9,9);
        if or(strcmp(difftype,'Bertallanfy'),strcmp(difftype,'BertallanfyConstr'))
            Q(9,9) = (Parameters.BRsigma.Value*mpred(9))^2;
            if isinf(Q(9,9))
                Q(9,9) = (Parameters.BRsigma.Value)^2;
                Crash = 1;
            end
        elseif strcmp(difftype,'Sigmoid')
            Q(9,9) = (Parameters.Sigmsigma.Value*mpred(9))^2;
            if isinf(Q(9,9))
                Q(9,9) = (Parameters.Sigmsigma.Value)^2;
                Crash = 1;
                'crash'
            end
%             Q(9,9) = (Parameters.SigmaRW.Value/(mu*exp(mtemp(9))/(1+exp(mtemp(9)))^2))^2;
        elseif or(strcmp(difftype,'Add'),strcmp(difftype,'AddConstr'))
            Q(9,9) = (Parameters.SigmaRW.Value)^2;
        else
            die
        end
        Cov = Cov + (Jacobian*Cov+Cov*Jacobian'+Q)*TStep;
        
        try
            if not(Crash)
                [V,D] = eig(Cov);
                temp = real(eig(Cov));
                temp = max(temp,0);
                D = diag(temp);
                temp = ( V*D*V' + (V*D*V')')/2;
                Cov = temp;
            end
        catch
            Crash = 1;
        end
        record(:,IndDiscr)= mtemp;
        
%         mpred(1:8) = max(0,mpred(1:8));
    end
    
% try
%     Res.deltabetas = deltabetas;
% end
    Res.m = mpred;
    Res.Cov = Cov;
    Res.Model = Model;
    Res.Crash = Crash;