function Res = SEIR_EKF_projection(Data,Model,m,Cov,NbIts,IndTime,Parameters)




    
   

    if Parameters.RunningMif 
        mpred = m; 
    end
    mpred(1:5,1) = max(m(1:5),0);
    mpred(1:5,1) = min(m(1:5),Parameters.TotalPopulation);
    mpred(6,1) = m(6);
    mpred(5) = 0;
    if length(m) == 6
        m(7) = 0;
    end       
    mpred(7,1) = m(7);
    
    Cov(5,:) = 0;
    Cov(:,5) = 0;
    TStep = Parameters.ComputationTStep;
    TotPop = Parameters.TotalPopulation;
    for IndDiscr = 1:NbIts
        mpred(1:5,1) = max(mpred(1:5,1),0);
        mpred(1:5,1) = min(mpred(1:5,1),Parameters.TotalPopulation);
        mtemp = mpred;
%         if sum(mpred(1:5)<0)
%             disp('pb signe')
%             die
%         end
        beta = exp(mtemp(6));
        mpred(1) = mpred(1) + (-beta*mtemp(1)*mtemp(3)/TotPop)*TStep;
        mpred(2) = mpred(2) + ( beta*mtemp(1)*mtemp(3)/TotPop- Parameters.km1.Value^-1*mtemp(2))*TStep;
        mpred(3) = mpred(3) + ( Parameters.km1.Value^-1*mtemp(2) - Parameters.gammam1.Value^-1*mtemp(3))*TStep;
        mpred(4) = mpred(4) + ( Parameters.gammam1.Value^-1*mtemp(3))*TStep;
        mpred(5) = mpred(5) + ( Parameters.km1.Value^-1*mtemp(2))*TStep;
        
        if strcmp(Parameters.DiffusionType,'Add')
            mpred(6) = mtemp(6);
%         elseif strcmp(Parameters.DiffusionType,'IBM')
%             mpred(6) = mtemp(6) + mtemp(7)*TStep;
%             mpred(7) = mtemp(7);
%         elseif strcmp(Parameters.DiffusionType,'SVO')
%             mpred(6) = mtemp(6);
%             mpred(7) = mtemp(7);
        end
            
        if Parameters.RunningMif 
            names = Parameters.Names.Estimated;
            tempn = 7+length(names); 
            Jacobian = zeros(tempn,tempn);
        else
            Jacobian = zeros(7,7);
        end
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
        Jacobian(5,2) =  Parameters.km1.Value^-1;
        if strcmp(Parameters.DiffusionType,'IBM')
            Jacobian(6,7) =  1;
        end    
        if Parameters.RunningMif 
            if Parameters.km1.Estimated 
                Jacobian(2,7+Parameters.km1.Index) =  Parameters.km1.Value^-2*mtemp(2);
                Jacobian(3,7+Parameters.km1.Index) = -Parameters.km1.Value^-2*mtemp(2);
                Jacobian(5,7+Parameters.km1.Index) = -Parameters.km1.Value^-2*mtemp(2);
            end
            if Parameters.gammam1.Estimated 
                Jacobian(3,7+Parameters.gammam1.Index) =  Parameters.gammam1.Value^-2*mtemp(3);
                Jacobian(4,7+Parameters.gammam1.Index) = -Parameters.gammam1.Value^-2*mtemp(3);
            end            
        end
        
%         if strcmp(Parameters.DiffusionType,'OUD')
%             Jacobian(6,6) = -Parameters.KappaOU.Value;
%         end
        Q = zeros(7,7);
        if strcmp(Parameters.DiffusionType,'Add')
            Q(6,6) = (Parameters.SigmaRW.Value)^2;
%         elseif strcmp(Parameters.DiffusionType,'OUD')
%             Q(6,6) = (Parameters.SigmaOU.Value)^2;
%         elseif strcmp(Parameters.DiffusionType,'IBM')
%             Q(7,7) = (Parameters.SigmaRW.Value)^2;
%         elseif strcmp(Parameters.DiffusionType,'IBM')
%             Q(6,6) = (mpred(7))^2;
%             Q(7,7) = (Parameters.SigmaRW.Value)^2;
        end
        if Parameters.RunningMif 
            names = Parameters.Names.Estimated;
            tempn = 7+length(names);
            Q = zeros(tempn,tempn);
            for i = 8:7+length(names)
                if Parameters.(names{i-7}).Init
                    Q(i,i) = Parameters.MIFVarInit;
                else
                    Q(i,i) = Parameters.MIFVarNotInit;
                end
            end
        end
        Cov2 = Cov + (Jacobian*Cov+Cov*Jacobian'+Q)*TStep;
                
        
        if sum(not(isreal(Cov2)))
%             disp('pb')
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
    if Parameters.RunningMif 
        temp = zeros(1,7+length(Parameters.Names.Estimated));
    else
        temp = zeros(1,7);
    end
    temp(1,5) = 1;
    Model.ObservationJacobian = {};
    Model.ObservationMeasurementNoise = {};
    try
        coeff = Parameters.MultCoeff.Value/10;
    catch
        coeff = 1;
    end
    for i = 1:length(Data.Instants)
        Model.ObservationJacobian{i} = temp;
        Model.ObservationMeasurementNoise{i} = (Parameters.SigmaObs.Value*Data.Observations(5,i))^2; % Additional approximation here due to kalman gaussian approximation
    end

    Res.m = mpred;
    Res.Cov = Cov;
    Res.Model = Model;
    Res.Crash = 0;