function Res = SIR_UKF_projection(Data,Model,m,Cov,NbIts,IndTime,Parameters)

    TStep = Parameters.ComputationTStep;

    L = length(Parameters.InitialState);
    Parameters.UKFalpha.Value = 1;%1;
    Parameters.UKFkappa.Value = 2;
    Parameters.UKFbeta.Value = 0;
    
    mpred = m;
    mpred(3) = 0;
    Cov(3,:) = 0;
    Cov(:,3) = 0;
    Cov(3,3) = 0.0001;
    
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
    plind = 3;
    plot((IndTime-1)*NbIts,m(plind),'og')
    hold on
    plot((IndTime-1)*NbIts ,mpred(plind)+sqrt(Cov(plind,plind)),'+y')
    plot((IndTime-1)*NbIts ,mpred(plind)-sqrt(Cov(plind,plind)),'+y')
    
    
    for IndDiscr = 1:NbIts
%         IndDiscr
%         if not(mean(mean(Cov == Cov'))==1)
%             [num2str(IndTime) ' ' num2str(IndDiscr) ' ' 'notsym']
%             Cov-Cov'
%         end
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
        
        beta = exp(xtmp(4,:));
        
        ftmp = zeros(size(xtmp));
        ftmp(1,:) =  -beta.*xtmp(1,:).*xtmp(2,:)/Parameters.PopSize;
        ftmp(2,:) =   beta.*xtmp(1,:).*xtmp(2,:)/Parameters.PopSize -  gamma*xtmp(2,:);
        ftmp(3,:) =   beta.*xtmp(1,:).*xtmp(2,:)/Parameters.PopSize ;
        ftmp(4,:) =   zeros(size(ftmp(3,:)));
            
        [fmean fCov] = UT(ftmp,W);
        
       
        mpred = mpred + fmean*TStep;
        
%         mpred(1)
%         mpred(2)
%         mpred(3)
%         mpred(4)
%         die
        
        Q = zeros(L,L);
        Q(4,4) = (Parameters.SigmaRW.Value)^2;
        Q = (Parameters.SigmaRW.Value)^2*eye(L,L)*10;

        
%         xpreds = xtmp + ftmp*TStep;
        
%         Cov = Cov + (xpreds-repmat(mpred,1,2*L+1))*diag(Wcs)*((xpreds-repmat(mpred,1,2*L+1))');

        n = length(m);
        dCov = zeros(n,n);
        for k = 1:2*n+1
            tmp = (ftmp(:,k)-fmean)*((xtmp(:,k) - xtmpmean)');
            dCov = dCov +   W(k)*(tmp+tmp');
        end
        dCov = dCov +  Q;
        
        Cov = Cov + dCov*TStep;
%         
%         [V,D] = eig(Cov);
%         temp = real(eig(Cov));
%         temp = max(temp,0.001);
%         D = diag(temp);
%         vp = temp;
%         temp = ( V*D*V' + (V*D*V')')/2;
%         Cov = temp;
        
        
%         Cov = (Cov+Cov')/2;
%         [V,D] = eig(Cov);
%         temp = real(eig(Cov));
%         temp = max(temp,0.001);
%         D = diag(temp);
%         vp = temp;
%         temp = ( V*D*V' + (V*D*V')')/2;
%         Cov = temp;
%         [IndTime Cov(3,3)]

        plot((IndTime-1)*NbIts + IndDiscr,mpred(plind),'og')
        plot((IndTime-1)*NbIts + IndDiscr,mpred(plind)+sqrt(Cov(plind,plind)),'+y')
        plot((IndTime-1)*NbIts + IndDiscr,mpred(plind)-sqrt(Cov(plind,plind)),'+y')
%         eig(ftmp*BigW*(xtmp') + xtmp*BigW*(ftmp'))
        
        record(:,IndDiscr)= mpred;
        
%         IndDiscr 
%         xpreds
%         Cov
%         die
        
    end
    mpred
    Cov
%     die

    Res.m = mpred;
    Res.Cov = Cov;%1/2*(Cov+Cov');
    Res.Model = Model;
