function Res = Dengue_EKF_projection(Data,Model,m,Cov,NbIts,IndTime,Parameters)

    TempVariables = m;
    TempVariables(11)=0;
    Variables = TempVariables; 
    
    TStep = Parameters.ComputationTStep;

   
    
    rep1 = Parameters.rep1.Value;
    rep2 = Parameters.rep2.Value;
    phi = Parameters.phi.Value;
    record = zeros(11,NbIts);
    
  
        
    for IndDiscr = 1:NbIts;
     
           
        MuB = Parameters.MuB.Value;
        MuD = Parameters.MuD.Value;
        Tot = Parameters.PopulationSize.Value;
        r0  = exp(TempVariables(10));
        v = Parameters.vm1.Value^-1;
        e = Parameters.e.Value;
        d = Parameters.d.Value;
        seas = e*(1+sin(2*pi*((Data.Instants(IndTime)+IndDiscr)*Parameters.ComputationTStep/12+d )));
        ade = Parameters.ade.Value;
        eta = Parameters.eta.Value;
        q = Parameters.qm1.Value^-1;
        
        
        
        TempVariables(1) = TempVariables(1) + (MuB*Tot - MuD*Variables(1) - r0/Tot*v*seas.*(Variables(2)+ade*Variables(8)+eta).*Variables(1) - r0/Tot*v*seas.*(Variables(3)+ade*Variables(9)+eta).*Variables(1))*TStep;
        TempVariables(2) = TempVariables(2) + (- MuD*Variables(2) + r0/Tot*v*seas.*(Variables(2)+ade*Variables(8)+eta).*Variables(1)-Variables(2)*v)*TStep;
        TempVariables(3) = TempVariables(3) + (- MuD*Variables(3) + r0/Tot*v*seas.*(Variables(3)+ade*Variables(9)+eta).*Variables(1)-Variables(3)*v)*TStep;
        TempVariables(4) = TempVariables(4) + (- MuD*Variables(4) + Variables(2)*v - Variables(4)*q)*TStep;
        TempVariables(5) = TempVariables(5) + (- MuD*Variables(5) + Variables(3)*v - Variables(5)*q)*TStep;
        TempVariables(6) = TempVariables(6) + (- MuD*Variables(6) + Variables(4)*q  - r0/Tot*v*seas.*(Variables(3)+ade*Variables(9)+eta).*Variables(6))*TStep;
        TempVariables(7) = TempVariables(7) + (- MuD*Variables(7) + Variables(5)*q  - r0/Tot*v*seas.*(Variables(2)+ade*Variables(8)+eta).*Variables(7))*TStep;
        TempVariables(8) = TempVariables(8) + (- MuD*Variables(8) + r0/Tot*v*seas.*(Variables(3)+ade*Variables(9)+eta).*Variables(6)-Variables(8)*v)*TStep;
        TempVariables(9) = TempVariables(9) + (- MuD*Variables(9) + r0/Tot*v*seas.*(Variables(2)+ade*Variables(8)+eta).*Variables(7)-Variables(9)*v)*TStep;
        TempVariables(11) = TempVariables(11) + (r0/Tot*v*seas.*(Variables(3)+ade*Variables(9)+eta).*Variables(6)+r0/Tot*v*seas.*(Variables(2)+ade*Variables(8)+eta).*Variables(7))*TStep;
               
    
        TempVariables(1:9) = max(TempVariables(1:9),eps);
        TempVariables(11) = max(TempVariables(11),eps);
        
        Variables = TempVariables;  
        
 
        r0der = exp(TempVariables(10));
        
        
        Jacobian = zeros(11,11);
        Jacobian(1,1) =    - MuD - r0/Tot*v*seas*(Variables(2)+ade*Variables(8)+eta) - r0/Tot*v*seas*(Variables(3)+ade*Variables(9)+eta);
        Jacobian(1,2) =    - r0/Tot*v*seas*Variables(1) ;
        Jacobian(1,3) =    - r0/Tot*v*seas*Variables(1);
        Jacobian(1,4) =    0;
        Jacobian(1,5) =    0;
        Jacobian(1,6) =    0;
        Jacobian(1,7) =    0;
        Jacobian(1,8) =    - r0/Tot*v*seas*ade*Variables(1);
        Jacobian(1,9) =    - r0/Tot*v*seas*ade*Variables(1);
        Jacobian(1,10) =   - r0der/Tot*v*seas*(Variables(2)+ade*Variables(8)+eta).*Variables(1) - r0der/Tot*v*seas*(Variables(3)+ade*Variables(9)+eta).*Variables(1);
        Jacobian(2,1) =   r0/Tot*v*seas*(Variables(2)+ade*Variables(8)+eta);
        Jacobian(2,2) =   - MuD + r0/Tot*v*seas*Variables(1)-v;
        Jacobian(2,8) =   r0/Tot*v*seas*ade*Variables(1);
        Jacobian(2,10) =  r0der/Tot*v*seas*(Variables(2)+ade*Variables(8)+eta).*Variables(1);
        Jacobian(3,1) =   r0/Tot*v*seas*(Variables(3)+ade*Variables(9)+eta);
        Jacobian(3,3) =   - MuD + r0/Tot*v*seas*Variables(1)-v;
        Jacobian(3,9) =   r0/Tot*v*seas*ade*Variables(1);
        Jacobian(3,10) =  r0der/Tot*v*seas*(Variables(3)+ade*Variables(9)+eta).*Variables(1);
        Jacobian(4,2) =  v;
        Jacobian(4,4) =  - MuD -q;
        Jacobian(5,3) =  v;
        Jacobian(5,5) =  - MuD -q;
        Jacobian(6,3) = - r0/Tot*v*seas.*Variables(6);
        Jacobian(6,4) = q;
        Jacobian(6,6) = - MuD  - r0/Tot*v*seas*(Variables(3)+ade*Variables(9)+eta);
        Jacobian(6,9) = - r0/Tot*v*seas*ade*Variables(6);
        Jacobian(6,10) = - r0der/Tot*v*seas*(Variables(3)+ade*Variables(9)+eta).*Variables(6);
        Jacobian(7,2) = - r0/Tot*v*seas.*Variables(7);
        Jacobian(7,5) = q;
        Jacobian(7,7) = - MuD  - r0/Tot*v*seas*(Variables(2)+ade*Variables(8)+eta);
        Jacobian(7,8) = - r0/Tot*v*seas*ade*Variables(7);
        Jacobian(7,10) = - r0der/Tot*v*seas*(Variables(2)+ade*Variables(8)+eta).*Variables(7);
        Jacobian(8,3) =  r0/Tot*v*seas*Variables(6);
        Jacobian(8,6) =  r0/Tot*v*seas*(Variables(3)+ade*Variables(9)+eta);
        Jacobian(8,8) =  - MuD -v;
        Jacobian(8,9) =  r0/Tot*v*seas*ade*Variables(6);
        Jacobian(8,10) =  r0der/Tot*v*seas*(Variables(3)+ade*Variables(9)+eta).*Variables(6);
        Jacobian(9,2) =  r0/Tot*v*seas*Variables(7);
        Jacobian(9,7) =  r0/Tot*v*seas*(Variables(2)+ade*Variables(8)+eta);
        Jacobian(9,9) =  - MuD -v;
        Jacobian(9,10) =  r0der/Tot*v*seas*(Variables(2)+ade*Variables(8)+eta).*Variables(7);
        Jacobian(11,2) = r0/Tot*v*seas.*Variables(7);
        Jacobian(11,3) = r0/Tot*v*seas.*Variables(6);
        Jacobian(11,6) = r0/Tot*v*seas.*(Variables(3)+ade*Variables(9)+eta);
        Jacobian(11,7) = r0/Tot*v*seas.*(Variables(2)+ade*Variables(8)+eta);  
        Jacobian(11,8) = r0/Tot*v*seas.*ade.*Variables(7);
        Jacobian(11,9) = r0/Tot*v*seas.*ade.*Variables(6);
        Jacobian(11,10) = r0der/Tot*v*seas.*(Variables(3)+ade*Variables(9)+eta).*Variables(6)+r0der/Tot*v*seas.*(Variables(2)+ade*Variables(8)+eta).*Variables(7);
        Jacobian(11,11) = 0;
        
        Q = zeros(11,11);
        
        Q(10,10) = (Parameters.SigmaRW.Value)^2;
        Cov = Cov + (Jacobian*Cov+Cov*Jacobian'+Q)*TStep;
        record(:,IndDiscr)= Variables;
        
%         mpred(1:8) = max(0,mpred(1:8));
    end
    
    
    Model.ObservationMeasurementNoise{IndTime} = rep2*(1.0-rep2)*rep1*(Variables(11)) + (rep2*phi*rep1*Variables(11))^2;

    
% try
%     Res.deltabetas = deltabetas;
% end
    Res.m = Variables;
    Res.Cov = Cov;
    Res.Model = Model;
    Res.Crash = 0;