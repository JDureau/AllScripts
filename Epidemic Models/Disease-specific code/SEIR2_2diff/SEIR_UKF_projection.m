function Res = SEIR_UKF_projection(Data,Model,m,Cov,NbIts,IndTime,Parameters)

    TStep = Parameters.ComputationTStep;

    mpred = m;
    
    
    gamma = Parameters.gammam1.Value^-1;
    
    L = length(m);
    UKFalpha = Parameters.UKFalpha.Value;
    UKFkappa = Parameters.UKFkappa.Value;
    UKFbeta = Parameters.UKFbeta.Value;
    c = UKFalpha^2*(L+UKFkappa);
    lambda = c-L;
    
    record = zeros(L,NbIts);

    Wms = [lambda/(L+lambda) repmat(1/(2*(L+lambda)),1,2*L)];
    Wcs = [lambda/(L+lambda)+(1-UKFalpha^2+UKFbeta) repmat(1/(2*(L+lambda)),1,2*L)];
    
    BigW = (eye(2*L+1)-repmat(Wms',1,2*L+1))*diag(Wcs)*((eye(2*L+1)-repmat(Wms',1,2*L+1))');
    
%     clf
    plind = 5;
    plot(IndTime*NbIts + 1,m(plind),'ob')
    hold on
    
    mpred(5) = 0;
%     Cov(5,:) = 0.1;
%     Cov(:,5) = 0.1;
    TStep = Parameters.ComputationTStep;
    TotPop = Parameters.TotalPopulation;
    
    for IndDiscr = 1:NbIts
        IndDiscr
        xtmp = repmat(mpred,1,2*L+1);

        choltmp = chol(Cov);
        
        xtmp = xtmp + sqrt(c)*[zeros(size(mpred)) choltmp -choltmp];
        
        
        xtmp(1:5,:) = max(xtmp(1:5,:),0);
        xtmp(1:5,:) = min(xtmp(1:5,:),Parameters.TotalPopulation);
        mtemp = mpred;
%         if sum(mpred(1:5)<0)
%             disp('pb signe')
%             die
%         end
        beta = exp(xtmp(6,:));
        
        ftmp = zeros(size(xtmp));
        ftmp(1,:) =  -beta.*xtmp(1,:).*xtmp(3,:)/TotPop;
        ftmp(2,:) =   beta.*xtmp(1,:).*xtmp(3,:)/TotPop - Parameters.km1.Value^-1*xtmp(2,:);
        ftmp(3,:) =   Parameters.km1.Value^-1*xtmp(2,:) - Parameters.gammam1.Value^-1*xtmp(3,:);
        ftmp(4,:) =   Parameters.gammam1.Value^-1*xtmp(3,:);
        ftmp(5,:) =   Parameters.km1.Value^-1*xtmp(2,:);
        ftmp(6,:) =   0;
        
        mpred = mpred + ftmp*(Wms')*TStep;
        
        Q = zeros(L,L);
        Q(6,6) = (Parameters.SigmaRW.Value)^2;
        
        
        xpreds = xtmp + ftmp*TStep;
        
%         Cov = Cov + (xpreds-repmat(mpred,1,2*L+1))*diag(Wcs)*((xpreds-repmat(mpred,1,2*L+1))');

        
        Cov = Cov + (ftmp*BigW*(xtmp') + xtmp*BigW*(ftmp')+Q)*TStep;
        


        plot(IndTime*NbIts + IndDiscr,mpred(plind),'ob')
        plot(IndTime*NbIts + IndDiscr,mpred(plind)+sqrt(Cov(plind,plind)),'+r')
        plot(IndTime*NbIts + IndDiscr,mpred(plind)-sqrt(Cov(plind,plind)),'+r')
%         eig(ftmp*BigW*(xtmp') + xtmp*BigW*(ftmp'))
        
        record(:,IndDiscr)= mpred;
        
    end
    

    Res.m = mpred;
    Res.Cov = 1/2*(Cov+Cov');
    Res.Model = Model;
