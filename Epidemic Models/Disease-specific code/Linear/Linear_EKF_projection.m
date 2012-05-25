function Res = Linear_EKF_projection(Data,Model,m,Cov,NbIts,IndTime,Parameters)

mpred = m;
 
    TStep = Parameters.ComputationTStep;
    
    temp = zeros(1,4);
    temp(1,3) = 1;
    Model.ObservationJacobian = {};
%     Model.ObservationMeasurementNoise = {};
   
    
    for i = 1:length(Data.Instants)
        Model.ObservationJacobian{i} = temp;
%         Model.ObservationMeasurementNoise{i} = Model.ObservationMeasurementNoise{i};
    end
    
    plind = 1;
    plot(IndTime*NbIts,m(plind),'ob')

%     plot(IndTime*NbIts + 1,exp(m(plind)),'ob')
    hold on
    plot(IndTime*NbIts ,mpred(plind)+sqrt(Cov(plind,plind)),'+r')
    plot(IndTime*NbIts ,mpred(plind)-sqrt(Cov(plind,plind)),'+r')
    
    
    for IndDiscr = 1:NbIts;
        mtemp = mpred;
        if not(mean(mean(Cov == Cov')))
            [IndTime IndDiscr 'notsym']
        end
        

        mpred(1) = mpred(1) + (0.05*mtemp(1)+0.06*mtemp(2)+0.07*mtemp(3)+0.08*mtemp(4))*TStep;
        mpred(2) = mpred(2) + (0.09*mtemp(1)+0.10*mtemp(2)+0.11*mtemp(3)+0.12*mtemp(4))*TStep;
        mpred(3) = mpred(3) + (0.13*mtemp(1)+0.14*mtemp(2)+0.15*mtemp(3)+0.16*mtemp(4))*TStep;
        mpred(4) = mpred(4) + (0.17*mtemp(1)+0.18*mtemp(2)+0.19*mtemp(3)+0.20*mtemp(4))*TStep;

       
      
        Jacobian = zeros(4,4);
        Jacobian(1,1) = 0.05;
        Jacobian(1,2) = 0.06;
        Jacobian(1,3) = 0.07;
        Jacobian(1,4) = 0.08;
        Jacobian(2,1) = 0.09;
        Jacobian(2,2) = 0.10;
        Jacobian(2,3) = 0.11;
        Jacobian(2,4) = 0.12;
        Jacobian(3,1) = 0.13;
        Jacobian(3,2) = 0.14;
        Jacobian(3,3) = 0.15;
        Jacobian(3,4) = 0.16;
        Jacobian(4,1) = 0.17;
        Jacobian(4,2) = 0.18;
        Jacobian(4,3) = 0.19;
        Jacobian(4,4) = 0.20;
        

       
% exp(mpred(plind))
%         plot(IndTime*NbIts + IndDiscr,exp(mpred(plind)),'ob')
%         plot(IndTime*NbIts + IndDiscr,exp(mpred(plind)+sqrt(Cov(plind,plind))),'+r')
%         plot(IndTime*NbIts + IndDiscr,exp(mpred(plind)-sqrt(Cov(plind,plind))),'+r') 
%         eig(ftmp*BigW*(xtmp') + xtmp*BigW*(ftmp'))
        
        Q = zeros(4,4);
        Q(1,1) = (Parameters.SigmaRW.Value)^2;
        Q = (Parameters.SigmaRW.Value)^2*eye(4,4);

        Cov = Cov + (Jacobian*Cov+Cov*(Jacobian')+Q)*TStep;
       
        plot(IndTime*NbIts + IndDiscr,mpred(plind),'ob')
        plot(IndTime*NbIts + IndDiscr,mpred(plind)+sqrt(Cov(plind,plind)),'+r')
        plot(IndTime*NbIts + IndDiscr,mpred(plind)-sqrt(Cov(plind,plind)),'+r')
%  mpred(1)
%         die
%         Cov(3,3)
      
%         IndDiscr
%         mpred
%         Cov
%         die
    end
%     die
   
    
    Res.m = mpred;
    Res.Cov = Cov;
    Res.Model = Model;
    Res.Crash = 0;
    