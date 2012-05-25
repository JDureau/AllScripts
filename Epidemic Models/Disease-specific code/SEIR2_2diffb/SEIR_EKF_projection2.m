function Res = SEIR_EKF_projection2(m,Cov,NbIts,IndTime,Parameters)

    mpred = m;
    TStep = Parameters.ComputationTStep;
    TotPop = Parameters.TotalPopulation;
    mpred(5) = log(eps);
    Cov(5,5) = 0;
    Cov(:,5) = 0;
    Cov(5,5) = 0;
    for IndDiscr = 1:NbIts;
        mtemp = mpred;
        beta = exp(mtemp(6));
        expmtemp = exp(mtemp);
        mpred(1) = mpred(1) + (-beta*expmtemp(3)/TotPop)*TStep;
        mpred(2) = mpred(2) + ( beta*expmtemp(1)*expmtemp(3)/(TotPop*expmtemp(2))- Parameters.k.Value)*TStep;
        mpred(3) = mpred(3) + ( Parameters.k.Value*expmtemp(2)/expmtemp(3) - Parameters.gamma.Value)*TStep;
        mpred(4) = mpred(4) + ( Parameters.gamma.Value*expmtemp(3)/expmtemp(4))*TStep;
        mpred(5) = mpred(5) + ( Parameters.k.Value*expmtemp(2)/expmtemp(5))*TStep;
        mpred(6) = mtemp(6);
       
        expmpred = exp(mpred);
        Jacobian = zeros(6,6);
        Jacobian(1,1) = -beta*expmpred(3)/TotPop*expmpred(1);
        Jacobian(1,3) = -beta*expmpred(1)/TotPop*expmpred(3);
        Jacobian(1,6) = -beta*expmpred(3)*mpred(1)/TotPop*expmpred(6);
        Jacobian(2,1) =  beta*expmpred(3)/TotPop*expmpred(1);
        Jacobian(2,2) = -Parameters.k.Value*expmpred(2);
        Jacobian(2,3) =  beta*expmpred(1)/TotPop*expmpred(3);
        Jacobian(2,6) =  beta*expmpred(3)*expmpred(1)/TotPop*expmpred(6);
        Jacobian(3,2) =  Parameters.k.Value*expmpred(2);
        Jacobian(3,3) = -Parameters.gamma.Value*expmpred(3);
        Jacobian(4,3) =  Parameters.gamma.Value*expmpred(3);
        Jacobian(5,2) =  Parameters.k.Value*expmpred(2);
%         if strcmp(Parameters.DiffusionType,'OUD')
%             Jacobian(6,6) = -Parameters.KappaOU.Value;
%         end
        Q = zeros(6,6);
        if strcmp(Parameters.DiffusionType,'Add')
            Q(6,6) = (Parameters.SigmaRW.Value)^2;
        elseif strcmp(Parameters.DiffusionType,'OUD')
            Q(6,6) = (Parameters.SigmaOU.Value)^2;
        end
        Cov = Cov + (Jacobian*Cov+Cov*Jacobian'+Q)*TStep;
        
       
        try sum(eig(Cov)<-10^-3)
            i = 1;
        catch
            disp('pb')
%             die
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