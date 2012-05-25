function Res = SIRsimple_EKF_projection(Data,Model,m,Cov,NbIts,IndTime,Parameters)

   
    mpred = m;
    mpred(3) = 0;
    Cov(3,:) = 0;
    Cov(:,3) = 0;
    
    TStep = Parameters.ComputationTStep;
    rep2 = Parameters.rep2.Value;
    rep1 = Parameters.rep1;
    phi = Parameters.phi.Value;
    
    
    
    p_t = Parameters.p_t;
    v = (Parameters.vm1.Value/7)^(-1);
    
    for IndDiscr = 1:NbIts;
        
%         for i = 1:4
%             disp(mpred(i));
%         end
%         disp(['t=' num2str(IndTime + IndDiscr*TStep)])

        mtemp = mpred;
        
        r0 = exp(mtemp(4));
        
        mtemp(1:3) = max(0,mtemp(1:3));
        mtemp(1:3) = min(p_t,mtemp(1:3));
        
%         disp(['S: ' num2str(mtemp(1))])
%         disp(['I: ' num2str(mtemp(2))])
%         disp(['r0: ' num2str(r0)])
%         disp(['v: ' num2str(v,15)])
%         disp(['p_t: ' num2str(p_t)])
%         disp(num2str(TStep,15))
%         
        mpred(1) = mpred(1) + ( -r0*v/p_t*mtemp(1)*mtemp(2) )*TStep; 
        mpred(2) = mpred(2) + (  r0*v/p_t*mtemp(1)*mtemp(2) - v * mtemp(2) )*TStep; 
        mpred(3) = mpred(3) + (  r0*v/p_t*mtemp(1)*mtemp(2) )*TStep;
          
        
        
      
        r0der = r0;
        Jacobian = zeros(4,4);
        Jacobian(1,1) = -r0*v/p_t*mtemp(2) ;
        Jacobian(1,2) = -r0*v/p_t*mtemp(1) ;
        Jacobian(1,4) = -r0der*v/p_t*mtemp(1)*mtemp(2) ;     
        Jacobian(2,1) =  r0*v/p_t*mtemp(2) ;
        Jacobian(2,2) =  r0*v/p_t*mtemp(1) - v;
        Jacobian(2,4) =  r0der*v/p_t*mtemp(1)*mtemp(2) ;   
        Jacobian(3,1) =  r0*v/p_t*mtemp(2) ;
        Jacobian(3,2) =  r0*v/p_t*mtemp(1) ;
        Jacobian(3,4) =  r0der*v/p_t*mtemp(1)*mtemp(2) ;   
        
        Q = zeros(4,4);
        Q(4,4) = (Parameters.vol__r0.Value)^2;
                
        Cov = Cov + (Jacobian*Cov+Cov*Jacobian'+Q)*TStep;
%         Cov = (Cov + Cov')/2;
%         disp(Cov)
%         try
%             [V,D] = eig(Cov);
%             temp = (eig(Cov));
%             if sum(temp<0)
%                 'to0'
%             end
%             temp = max(temp,0);
%             D = diag(temp);
%             temp = ( V*D*V' + (V*D*V')')/2;
%             Cov = temp;
%         catch
%            'yo' 
%         end
    end



    Model.ObservationMeasurementNoise{IndTime} = rep2*(1-rep2)*rep1*mpred(3) + (rep2*phi*rep1*mpred(3))^2;
    
    
    tmp = zeros(4,1);
    tmp(3,1) = rep2*rep1;
    Model.ObservationJacobian{IndTime} = tmp';

    Res.m = mpred;
    Res.Cov = Cov;
    Res.Model = Model;
    Res.Crash = 0;
    