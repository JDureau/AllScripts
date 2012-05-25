function Res = SIRsimf_EKF_projection(Data,Model,m,Cov,NbIts,IndTime,Parameters)

    mpred = m;
    mpred(7) = 0;
    mpred(8) = 0;
    Cov(7,:) = 0;
    Cov(:,7) = 0;
    Cov(8,:) = 0;
    Cov(:,8) = 0;
    
    TStep = Parameters.ComputationTStep;
    rep2 = Parameters.rep2.Value;
    rep1 = Parameters.rep1;
    phi = Parameters.phi.Value;
    iota = Parameters.iota.Value;
    e = Parameters.e.Value;
    d = Parameters.d.Value;
    
    p_t = Parameters.p_t;
    mu_b = Parameters.mu_b;
    mu_d = Parameters.mu_d;
    v = (Parameters.vm1.Value/7)^(-1);
    
    for IndDiscr = 1:NbIts;
        mtemp = mpred;
        
        beta = exp(mtemp(9));
        t = sum(Data.NbComputingSteps(1:IndTime-1)) + IndDiscr-1;
        a_year = 364/Parameters.ComputationTStep;
        seas = (1.0+ e*sin(2.0*pi*( t/a_year )+d*2.0*pi ));
        
        mtemp(1:8) = max(0,mtemp(1:8));
        mtemp(1:8) = min(p_t,mtemp(1:8));
        
        mpred(1) = mpred(1) + (-beta*v/p_t*seas*mtemp(1)*(mtemp(3)+mtemp(5)+iota) + mu_b*p_t - mu_d*mtemp(1))*TStep; 
        mpred(2) = mpred(2) + (-beta*v/p_t*seas*mtemp(2)*(mtemp(4)+mtemp(6)+iota) + mu_b*p_t - mu_d*mtemp(2))*TStep;
        mpred(3) = mpred(3) + ( beta*v/p_t*seas*mtemp(1)*(mtemp(3)+mtemp(5)+iota) - 2*v * mpred(3) - mu_d*mtemp(3))*TStep; 
        mpred(4) = mpred(4) + ( beta*v/p_t*seas*mtemp(2)*(mtemp(4)+mtemp(6)+iota) - 2*v * mpred(4) - mu_d*mtemp(4))*TStep;
        mpred(5) = mpred(5) + (  2*v * mpred(3) - 2*v * mpred(5) - mu_d*mtemp(5)) *TStep;
        mpred(6) = mpred(6) + (  2*v * mpred(4) - 2*v * mpred(6) - mu_d*mtemp(6)) *TStep;
        mpred(7) = mpred(7) + ( beta*v/p_t*seas*mtemp(1)*(mtemp(3)+mtemp(5)+iota)  ) *TStep;
        mpred(8) = mpred(8) + ( mu_d*mtemp(4) + mu_d*mtemp(6) + 2*v * mpred(6)) *TStep;
        
        
        
      
        betader = beta;
        Jacobian = zeros(9,9);
        Jacobian(1,1) = -beta*v/p_t*seas*(mtemp(3)+mtemp(5)+iota) - mu_d ;
        Jacobian(1,3) = -beta*v/p_t*seas*mtemp(1) ;
        Jacobian(1,5) = -beta*v/p_t*seas*mtemp(1) ;
        Jacobian(1,9) = -betader*v/p_t*seas*mtemp(1)*(mtemp(3)+mtemp(5)+iota) ;
        Jacobian(2,2) = -beta*v/p_t*seas*(mtemp(4)+mtemp(6)+iota) - mu_d;
        Jacobian(2,4) = -beta*v/p_t*seas*mtemp(2) ;
        Jacobian(2,6) = -beta*v/p_t*seas*mtemp(2) ;
        Jacobian(2,9) = -betader*v/p_t*seas*mtemp(2)*(mtemp(4)+mtemp(6)+iota) ;
        Jacobian(3,1) = beta*v/p_t*seas*(mtemp(3)+mtemp(5)+iota)  ;
        Jacobian(3,3) = beta*v/p_t*seas*mtemp(1)  - 2*v - mu_d;
        Jacobian(3,5) = beta*v/p_t*seas*mtemp(1)  ;
        Jacobian(3,9) = betader*v/p_t*seas*mtemp(1)*(mtemp(3)+mtemp(5)+iota) ;
        Jacobian(4,2) = beta*v/p_t*seas*(mtemp(4)+mtemp(6)+iota)  ;
        Jacobian(4,4) = beta*v/p_t*seas*mtemp(2)  - 2*v - mu_d;
        Jacobian(4,6) = beta*v/p_t*seas*mtemp(2)  ;
        Jacobian(4,9) = betader*v/p_t*seas*mtemp(2)*(mtemp(4)+mtemp(6)+iota) ;
        Jacobian(5,3) =   2*v;
        Jacobian(5,5) = - 2*v;
        Jacobian(6,4) =   2*v;
        Jacobian(6,6) = - 2*v;
        Jacobian(7,1) = beta*v/p_t*seas*(mtemp(3)+mtemp(5)+iota)  ;
        Jacobian(7,3) = beta*v/p_t*seas*mtemp(1)  ;
        Jacobian(7,5) = beta*v/p_t*seas*mtemp(1)  ;
        Jacobian(7,9) = betader*v/p_t*seas*mtemp(1)*(mtemp(3)+mtemp(5)+iota)  ;
        Jacobian(8,4) = mu_d ;
        Jacobian(8,6) = mu_d + 2*v ;
        
        
        Q = zeros(9,9);
        Q(9,9) = (Parameters.vol__r0.Value)^2;
                
        Cov = Cov + (Jacobian*Cov+Cov*Jacobian'+Q)*TStep;
        Cov = (Cov + Cov')/2;
        try
            [V,D] = eig(Cov);
            temp = real(eig(Cov));
            temp = max(temp,0);
            D = diag(temp);
            temp = ( V*D*V' + (V*D*V')')/2;
            Cov = temp;
        catch
           'yo' 
        end
    end


    tmp = [];
    if Data.ObservedVariables(1,IndTime) 
        tmp(end+1) = rep2*(1-rep2)*rep1*mpred(7) + (rep2*phi*rep1*mpred(7))^2;
    end
    if Data.ObservedVariables(2,IndTime) 
        tmp(end+1) = rep2*(1-rep2)*rep1*mpred(8) + (rep2*phi*rep1*mpred(8))^2;
    end
    Model.ObservationMeasurementNoise{IndTime} = diag(tmp);
    
    tmp = [];
    if Data.ObservedVariables(1,IndTime) 
        tmp2 = zeros(9,1);
        tmp2(7,1) = rep2*(1-rep2)*rep1;
        tmp(:,end+1) = tmp2;
    end
    if Data.ObservedVariables(2,IndTime) 
        tmp2 = zeros(9,1);
        tmp2(8,1) = rep2*(1-rep2)*rep1;
        tmp(:,end+1) = tmp2;
    end
    Model.ObservationJacobian{IndTime} = tmp';

    Res.m = mpred;
    Res.Cov = Cov;
    Res.Model = Model;
    Res.Crash = 0;
    