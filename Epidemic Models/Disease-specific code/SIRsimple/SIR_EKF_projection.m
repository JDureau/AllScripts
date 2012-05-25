function Res = SIR_EKF_projection(Data,Model,m,Cov,NbIts,IndTime,Parameters)

mpred = m;
    mpred(3) = 0;
    Cov(3,:) = 0;
    Cov(:,3) = 0;
    Cov(3,3) = 0.0001;

    TStep = Parameters.ComputationTStep;
    
    temp = zeros(1,4);
    temp(1,3) = 1;
    Model.ObservationJacobian = {};
%     Model.ObservationMeasurementNoise = {};
   
    
    for i = 1:length(Data.Instants)
        Model.ObservationJacobian{i} = temp;
%         Model.ObservationMeasurementNoise{i} = Model.ObservationMeasurementNoise{i};
    end
    
    plind = 3;
    plot((IndTime-1)*NbIts,m(plind),'ob')

%     plot(IndTime*NbIts + 1,exp(m(plind)),'ob')
    hold on
    plot((IndTime-1)*NbIts ,mpred(plind)+sqrt(Cov(plind,plind)),'+r')
    plot((IndTime-1)*NbIts ,mpred(plind)-sqrt(Cov(plind,plind)),'+r')
    
    
    for IndDiscr = 1:NbIts;
        mtemp = mpred;
        if not(mean(mean(Cov == Cov')))
            [IndTime IndDiscr 'notsym']
        end
        
        beta = exp(mtemp(4));

        
        mpred(1) = mpred(1) + (-beta*mtemp(1)*mtemp(2)/Parameters.PopSize)*TStep;
        mpred(2) = mpred(2) + ( beta*mtemp(1)*mtemp(2)/Parameters.PopSize- Parameters.gammam1.Value^(-1)*mtemp(2))*TStep;
        mpred(3) = mpred(3) + ( beta*mtemp(1)*mtemp(2)/Parameters.PopSize)*TStep;
        mpred(4) = mpred(4) ;
        
%         mpred(1)
%         mpred(2)
%         mpred(3)
%         mpred(4)
%         die
      
        betader = beta;
        Jacobian = zeros(4,4);
        Jacobian(1,1) = -mtemp(4)*mtemp(2)/Parameters.PopSize;
        Jacobian(1,2) = -mtemp(4)*mtemp(1)/Parameters.PopSize;
        Jacobian(1,4) = -betader*mtemp(1)*mtemp(2)/Parameters.PopSize;
        Jacobian(2,1) =  mtemp(4)*mtemp(2)/Parameters.PopSize;
        Jacobian(2,2) =  mtemp(4)*mtemp(1)/Parameters.PopSize- Parameters.gammam1.Value^(-1);
        Jacobian(2,4) =  betader*mtemp(1)*mtemp(2)/Parameters.PopSize;
        Jacobian(3,1) =  mtemp(4)*mtemp(2)/Parameters.PopSize;
        Jacobian(3,2) =  mtemp(4)*mtemp(1)/Parameters.PopSize;
        Jacobian(3,4) =  betader*mtemp(1)*mtemp(2)/Parameters.PopSize;
        

       
% exp(mpred(plind))
%         plot(IndTime*NbIts + IndDiscr,exp(mpred(plind)),'ob')
%         plot(IndTime*NbIts + IndDiscr,exp(mpred(plind)+sqrt(Cov(plind,plind))),'+r')
%         plot(IndTime*NbIts + IndDiscr,exp(mpred(plind)-sqrt(Cov(plind,plind))),'+r') 
%         eig(ftmp*BigW*(xtmp') + xtmp*BigW*(ftmp'))
        
        Q = zeros(4,4);
        Q(4,4) = (Parameters.SigmaRW.Value)^2;
        Cov = Cov + (Jacobian*Cov+Cov*Jacobian'+Q)*TStep;
        
        plot((IndTime-1)*NbIts + IndDiscr,mpred(plind),'ob')
        plot((IndTime-1)*NbIts + IndDiscr,mpred(plind)+sqrt(Cov(plind,plind)),'+r')
        plot((IndTime-1)*NbIts + IndDiscr,mpred(plind)-sqrt(Cov(plind,plind)),'+r')

%         Cov(3,3)
      
%         IndDiscr
%         mpred
%         Cov
%         die
        
    end
%     mpred
%     Cov
%     die
    
    Res.m = mpred;
    Res.Cov = Cov;
    Res.Model = Model;
    Res.Crash = 0;
    