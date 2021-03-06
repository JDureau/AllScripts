function Res = SEIRPlusObs_EKF_projection(m,Cov,NbIts,IndTime,Parameters)

    mpred(1:5,1) = max(m(1:5),0);
    mpred(1:5,1) = min(m(1:5),Parameters.TotalPopulation);
    mpred(6,1) = m(6);
    mpred(7,1) = min(m(7),0);
    mpred(5) = 0;
    Cov(5,:) = 0;
    Cov(:,5) = 0;
    TStep = Parameters.ComputationTStep;
    TotPop = Parameters.TotalPopulation;
    for IndDiscr = 1:NbIts
        mpred(1:5,1) = max(mpred(1:5,1),0);
        mpred(1:5,1) = min(mpred(1:5,1),Parameters.TotalPopulation);
        mtemp = mpred;
        if sum(mpred(1:5)<0)
            disp('pb signe')
            die
        end
        beta = exp(mtemp(6));
        PropObs = exp(mtemp(7));
        mpred(1) = mpred(1) + (-beta*mtemp(1)*mtemp(3)/TotPop)*TStep;
        mpred(2) = mpred(2) + ( beta*mtemp(1)*mtemp(3)/TotPop- Parameters.km1.Value^-1*mtemp(2))*TStep;
        mpred(3) = mpred(3) + ( Parameters.km1.Value^-1*mtemp(2) - Parameters.gammam1.Value^-1*mtemp(3))*TStep;
        mpred(4) = mpred(4) + ( Parameters.gammam1.Value^-1*mtemp(3))*TStep;
        mpred(5) = mpred(5) + ( Parameters.km1.Value^-1*mtemp(2)*PropObs)*TStep;
        mpred(6) = mtemp(6);
        mpred(7) = mtemp(7);
       
        Jacobian = zeros(7,7);
        Jacobian(1,1) = -beta*mpred(3)/TotPop;
        Jacobian(1,3) = -beta*mpred(1)/TotPop;
        Jacobian(1,6) = -beta*mpred(3)*mpred(1)/TotPop;
        Jacobian(2,1) =  beta*mpred(3)/TotPop;
        Jacobian(2,2) = -Parameters.km1.Value^-1;
        Jacobian(2,3) =  beta*mpred(1)/TotPop;
        Jacobian(2,6) =  beta*mpred(3)*mpred(1)/TotPop;
        Jacobian(3,2) =  Parameters.km1.Value^-1;
        Jacobian(3,3) = -Parameters.gammam1.Value^-1;
        Jacobian(4,3) =  Parameters.gammam1.Value^-1;
        Jacobian(5,2) =  Parameters.km1.Value^-1*PropObs;
        Jacobian(5,7) =  Parameters.km1.Value^-1*mtemp(2)*PropObs;
%         if strcmp(Parameters.DiffusionType,'OUD')
%             Jacobian(6,6) = -Parameters.KappaOU.Value;
%         end
        Q = zeros(7,7);
        if strcmp(Parameters.DiffusionType,'Add')
            Q(6,6) = (Parameters.SigmaRW.Value)^2;
            Q(7,7) = (Parameters.SigmaRWObs.Value)^2;
        elseif strcmp(Parameters.DiffusionType,'OUD')
            Q(6,6) = (Parameters.SigmaOU.Value)^2;
        end
        Cov2 = Cov + (Jacobian*Cov+Cov*Jacobian'+Q)*TStep;
                
        
        if sum(not(isreal(Cov2)))
            disp('pb')
        else
            Cov = Cov2;
        end
    end
%      for i=1:6
%             subplot(6,1,i)
%             hold on
%             plot(IndDiscr,mpred(i),'.')
%             plot(IndDiscr,mpred(i)+sqrt(Cov(i,i)),'.r')
%             plot(IndDiscr,mpred(i)-sqrt(Cov(i,i)),'.r')
%         end
    
    Res.m = mpred;
    Res.Cov = Cov;