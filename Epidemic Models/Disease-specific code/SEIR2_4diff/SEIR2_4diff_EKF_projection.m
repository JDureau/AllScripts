function Res = SEIR2_3diff_EKF_projection(Data,Model,m,Cov,NbIts,IndTime,Parameters)


    temp = zeros(2,14);
    temp(1,9) = 1;
    temp(2,10) = 1;
    Model.ObservationJacobian = {};
    Model.ObservationMeasurementNoise = {};
    try
        coeff = Parameters.MultCoeff.Value/10;
    catch
        coeff = 1;
    end
    for i = 1:length(Data.Instants)
        Model.ObservationJacobian{i} = temp;
        Model.ObservationMeasurementNoise{i} = diag([(Parameters.SigmaObs.Value*Data.Observations(9,i))^2    (Parameters.SigmaObs.Value*Data.Observations(10,i))^2]);
    end

    inds1 = [1 3 5 7];
    inds2 = [2 4 6 8];
    mpred(1:8,1) = max(m(1:8),0);
    mpred(inds1,1) = min(m(inds1),Parameters.TotalPopulation1);
    mpred(inds2,1) = min(m(inds2),Parameters.TotalPopulation2);
    mpred(11,1) = m(11);
    mpred(12,1) = m(12);
    mpred(13,1) = m(13);
    mpred(14,1) = m(14);
    mpred(9) = 0;
    mpred(10) = 0;
    Cov(9,:) = 0;
    Cov(:,9) = 0;
    Cov(10,:) = 0;
    Cov(:,10) = 0;
    TStep = Parameters.ComputationTStep;
    TotPop1 = Parameters.TotalPopulation1;
    TotPop2 = Parameters.TotalPopulation2;
    for IndDiscr = 1:NbIts
        
        mpred(1:8,1) = max(mpred(1:8,1),0);
        mpred(inds1,1) = min(mpred(inds1,1),Parameters.TotalPopulation1);
        mpred(inds2,1) = min(mpred(inds2,1),Parameters.TotalPopulation2);
       
        mtemp = mpred;
        beta11 = exp(mtemp(11));
        beta12 = exp(mtemp(13));
        beta22 = exp(mtemp(12));
        beta21 = exp(mtemp(14));
        
        % S
        mpred(1) = mpred(1) + (-beta11*mtemp(1)*mtemp(5)/TotPop1  -beta12*mtemp(1)*mtemp(6)/TotPop2)*TStep;
        mpred(2) = mpred(2) + (-beta22*mtemp(2)*mtemp(6)/TotPop2  -beta21*mtemp(2)*mtemp(5)/TotPop1)*TStep;
        % E
        mpred(3) = mpred(3) + ( beta11*mtemp(1)*mtemp(5)/TotPop1 + beta12*mtemp(1)*mtemp(6)/TotPop2  - Parameters.km1.Value^(-1)*mtemp(3))*TStep;
        mpred(4) = mpred(4) + ( beta22*mtemp(2)*mtemp(6)/TotPop2 + beta21*mtemp(2)*mtemp(5)/TotPop1  - Parameters.km1.Value^(-1)*mtemp(4))*TStep;
        % I
        mpred(5) = mpred(5) + ( Parameters.km1.Value^(-1)*mtemp(3) - Parameters.gammam1.Value^(-1)*mtemp(5))*TStep;
        mpred(6) = mpred(6) + ( Parameters.km1.Value^(-1)*mtemp(4) - Parameters.gammam1.Value^(-1)*mtemp(6))*TStep;
        % R
        mpred(7) = mpred(7) + ( Parameters.gammam1.Value^(-1)*mtemp(5))*TStep;
        mpred(8) = mpred(8) + ( Parameters.gammam1.Value^(-1)*mtemp(6))*TStep;
        %Inc
        mpred(9)  = mpred(9)  + ( Parameters.km1.Value^(-1)*mtemp(3))*TStep;
        mpred(10) = mpred(10) + ( Parameters.km1.Value^(-1)*mtemp(4))*TStep;
        
        
        mpred(1) = max(mpred(1),0);
        mpred(2) = max(mpred(2),0);
        mpred(3) = max(mpred(3),0);
        mpred(4) = max(mpred(4),0);
        mpred(5) = max(mpred(5),0);
        mpred(6) = max(mpred(6),0);
        mpred(1) = min(mpred(1),Parameters.TotalPopulation1);
        mpred(2) = min(mpred(2),Parameters.TotalPopulation2);
        mpred(3) = min(mpred(3),Parameters.TotalPopulation1);
        mpred(4) = min(mpred(4),Parameters.TotalPopulation2);
        mpred(5) = min(mpred(5),Parameters.TotalPopulation1);
        mpred(6) = min(mpred(6),Parameters.TotalPopulation2);
        
            
        
            
        Jacobian = zeros(14,14);
        Jacobian(1,1)  = -beta11*mpred(5)/TotPop1 -beta12*mpred(6)/TotPop2;
        Jacobian(1,5)  = -beta11*mpred(1)/TotPop1;
        Jacobian(1,6)  = -beta12*mpred(1)/TotPop2;
        Jacobian(1,11) = -beta11*mpred(1)*mpred(5)/TotPop1;
        Jacobian(1,13) = -beta12*mtemp(1)*mtemp(6)/TotPop2;
        
        Jacobian(2,2)  = -beta22*mpred(6)/TotPop2 -beta21*mpred(5)/TotPop1;
        Jacobian(2,6)  = -beta22*mpred(2)/TotPop2;
        Jacobian(2,5)  = -beta21*mpred(2)/TotPop1; 
        Jacobian(2,12) = -beta22*mpred(2)*mpred(6)/TotPop2;        
        Jacobian(2,14) = -beta21*mtemp(2)*mtemp(5)/TotPop1;
        
        Jacobian(3,1)  =  beta11*mpred(5)/TotPop1 + beta12*mpred(6)/TotPop2;
        Jacobian(3,3)  = -Parameters.km1.Value^(-1);
        Jacobian(3,5)  =  beta11*mpred(1)/TotPop1;
        Jacobian(3,6)  =  beta12*mpred(1)/TotPop2;
        Jacobian(3,11) =  beta11*mpred(1)*mpred(5)/TotPop1;
        Jacobian(3,13) =  beta12*mtemp(1)*mtemp(6)/TotPop2;
        
        Jacobian(4,2)  =  beta22*mpred(6)/TotPop2 + beta21*mpred(5)/TotPop1;
        Jacobian(4,4)  = -Parameters.km1.Value^(-1);
        Jacobian(4,6)  =  beta22*mpred(2)/TotPop2;
        Jacobian(4,5)  =  beta21*mpred(2)/TotPop1;
        Jacobian(4,12) =  beta22*mpred(2)*mpred(6)/TotPop2;
        Jacobian(4,14) =  beta21*mtemp(2)*mtemp(5)/TotPop1;
        
        Jacobian(5,3)  =  Parameters.km1.Value^(-1);
        Jacobian(5,5)  =  -Parameters.gammam1.Value^(-1);
        
        Jacobian(6,4)  =  Parameters.km1.Value^-1;
        Jacobian(6,6)  =  -Parameters.gammam1.Value^-1;
        
        Jacobian(7,5)  =  Parameters.gammam1.Value^(-1);
        
        Jacobian(8,6)  =  Parameters.gammam1.Value^(-1);
        
        Jacobian(9,3)  =  Parameters.km1.Value^-1;
        
        Jacobian(10,4) =  Parameters.km1.Value^-1;
        
           
        
        
%         if strcmp(Parameters.DiffusionType,'OUD')
%             Jacobian(6,6) = -Parameters.KappaOU.Value;
%         end
        Q = zeros(14,14);
        if strcmp(Parameters.DiffusionType,'Add')
            Q(11,11) = (Parameters.SigmaRW11.Value)^2;
            Q(12,12) = (Parameters.SigmaRW22.Value)^2;
            Q(13,13) = (Parameters.SigmaRW12.Value)^2;
            Q(14,14) = (Parameters.SigmaRW21.Value)^2;
        end
        Cov2 = Cov + (Jacobian*Cov+Cov*Jacobian'+Q)*TStep;
                
%         if IndDiscr == 2
%             (Cov2)
%             die
%         end
%         
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
    Res.Model = Model;
    Res.Crash = 0;