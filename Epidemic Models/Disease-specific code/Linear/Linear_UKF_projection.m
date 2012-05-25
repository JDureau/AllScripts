function Res = Linear_UKF_projection(Data,Model,m,Cov,NbIts,IndTime,Parameters)

    TStep = Parameters.ComputationTStep;

    if size(m,2)>1
        L = size(m,1);
        kappa = Parameters.UKFkappa.Value;
        W = zeros(2*L+1,1);
        W(1) = kappa / (L+kappa);
        for k = 1:2*L
            W(k+1) = 1/(2*(L+kappa));
        end
        [m frfrCov] = UT(m,W);
    end
        
    mpred = m;
    
   
     temp = zeros(1,2);
    temp(1,1) = 1;
    Model.ObservationJacobian = {};
    for i = 1:length(Data.Instants)
        Model.ObservationJacobian{i} = temp;
    end
    
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
    plind = 1;
    plot(IndTime*NbIts,m(plind),'og')
    hold on
    plot(IndTime*NbIts ,mpred(plind)+sqrt(Cov(plind,plind)),'+y')
    plot(IndTime*NbIts ,mpred(plind)-sqrt(Cov(plind,plind)),'+y')
    
    
    for IndDiscr = 1:NbIts
%         IndDiscr
        if not(mean(mean(Cov == Cov'))==1)
            [num2str(IndTime) ' ' num2str(IndDiscr) ' ' 'notsym']
            Cov-Cov'
        end
         if sum(eig(Cov)<0)
            'not def pos'
        end
      
%         [V,D] = eig(Cov);
%         temp = real(eig(Cov));
%         temp = max(temp,eps);
%         D = diag(temp);
%         temp = ( V*D*V' + (V*D*V')')/2;
%         Cov = temp;
        
        
%         xtmp = xtmp + sqrt(c)*[zeros(size(mpred)) choltmp -choltmp];
        kappa = Parameters.UKFkappa.Value;
        [xtmp,W] = SigmaPoints(mpred,Cov,kappa);
      
        [xtmpmean xtmpCov] = UT(xtmp,W);
   
        
        ftmp = zeros(size(xtmp));
        ftmp(1,:) =  0.05*xtmp(1,:)+0.06*xtmp(2,:)+0.07*xtmp(3,:)+0.08*xtmp(4,:);
        ftmp(2,:) =  0.09*xtmp(1,:)+0.10*xtmp(2,:)+0.11*xtmp(3,:)+0.12*xtmp(4,:);
        ftmp(3,:) =  0.13*xtmp(1,:)+0.14*xtmp(2,:)+0.15*xtmp(3,:)+0.16*xtmp(4,:);
        ftmp(4,:) =  0.17*xtmp(1,:)+0.18*xtmp(2,:)+0.19*xtmp(3,:)+0.20*xtmp(4,:);
            
        [fmean fCov] = UT(ftmp,W);
        
        
        
        
        mpred = mpred + fmean*TStep;
        
        
        
        Q = zeros(L,L);
        Q(1,1) = (Parameters.SigmaRW.Value)^2;
        Q = (Parameters.SigmaRW.Value)^2*eye(L,L);
        
     
        
%         xpreds = xtmp + ftmp*TStep;
        
%         Cov = Cov + (xpreds-repmat(mpred,1,2*L+1))*diag(Wcs)*((xpreds-repmat(mpred,1,2*L+1))');
       
        n = length(m);
        dCov = zeros(n,n);
        for k = 1:2*n+1
            tmp = (ftmp(:,k)-fmean)*((xtmp(:,k) - xtmpmean)');
            dCov = dCov +   W(k)*(tmp+tmp');
        end
%         tmp = cov(ftmp',xtmp');
        dCov = dCov + Q ;
       
%         dCov
%         Cov
        chol(Cov+dCov)
        die
        
        Cov = Cov + dCov*TStep;
        
        
       
        
%         [V,D] = eig(Cov);
%         temp = real(eig(Cov));
%         temp = max(temp,0.001);
%         D = diag(temp);
%         vp = temp;
%         temp = ( V*D*V' + (V*D*V')')/2;
%         Cov = temp;

        if sum(eig(Cov)<0)
            'not def pos'
        end
      
%         Cov = (Cov+Cov')/2;
%         
%         [IndTime Cov(3,3)]

        plot(IndTime*NbIts + IndDiscr,mpred(plind),'og')
        plot(IndTime*NbIts + IndDiscr,mpred(plind)+sqrt(Cov(plind,plind)),'+y')
        plot(IndTime*NbIts + IndDiscr,mpred(plind)-sqrt(Cov(plind,plind)),'+y')
%         mpred(1)
%        
%         die
%         eig(ftmp*BigW*(xtmp') + xtmp*BigW*(ftmp'))
        
        record(:,IndDiscr)= mpred;
        
%         IndDiscr 
%         xpreds
%         Cov
%         die
        
    end
%     die
   
    [msigm,W] = SigmaPoints(mpred,Cov,kappa);
    Res.msigm = msigm;
    Res.m = mpred;
    Res.Cov = 1/2*(Cov+Cov');
    Res.Model = Model;
