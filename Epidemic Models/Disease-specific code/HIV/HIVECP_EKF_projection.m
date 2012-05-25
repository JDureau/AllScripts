function Res = HIVECP_EKF_projection(m,Cov,NbIts,IndTime,Parameters)

    mpred = m;
    TStep = Parameters.ComputationTStep;

    
    
%     record = zeros(8,NbIts);
    for IndDiscr = 1:NbIts;
        mtemp = mpred;
        mpred(1) = mpred(1) + ( -Parameters.Beta.Value*mpred(1)*mpred(2) + Parameters.Gamma.Value*(mtemp(2)+mtemp(3)))*TStep;
        mpred(2) = mpred(2) + ( Parameters.Beta.Value*mpred(1)*mpred(2) - (Parameters.Alpha.Value+Parameters.Gamma.Value)*(mtemp(2)))*TStep;
        mpred(3) = mpred(3) + ( Parameters.Alpha.Value*(mtemp(2)) - Parameters.Gamma.Value*mtemp(3))*TStep;
        mpred(4) = mpred(2) + mpred(3);
    end

    disp([num2str(IndTime) num2str(mpred(4))])
    
    Res.m = mpred;
    Res.Cov = Cov;