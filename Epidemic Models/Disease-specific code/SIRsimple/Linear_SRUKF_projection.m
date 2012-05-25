function Res = Linear_SRUKF_projection(Data,Model,msigm,Cov,NbIts,IndTime,Parameters)


    TStep = Parameters.ComputationTStep;

    
    mpred = mean(msigm);
    
   
     temp = zeros(1,2);
    temp(1,1) = 1;
    Model.ObservationJacobian = {};
    for i = 1:length(Data.Instants)
        Model.ObservationJacobian{i} = temp;
    end
    
    L = size(Cov,1);
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
    plot(IndTime*NbIts,mpred(plind),'oc')
    hold on
    plot(IndTime*NbIts ,mpred(plind)+sqrt(Cov(plind,plind)),'+o')
    plot(IndTime*NbIts ,mpred(plind)-sqrt(Cov(plind,plind)),'+o')
    
    
    for IndDiscr = 1:NbIts
%         IndDiscr
        if not(mean(mean(Cov == Cov'))==1)
            [num2str(IndTime) ' ' num2str(IndDiscr) ' ' 'notsym']
            Cov-Cov'
        end
        if sum(eig(Cov)<0)
            'not def pos'
        end
         
        kappa = Parameters.UKFkappa.Value;
        W = zeros(2*L+1,1);
        W(1) = kappa / (L+kappa);
        for k = 1:2*L
            W(k+1) = 1/(2*(L+kappa));
        end
        
        [xm Cov] = UT(msigm,W);
        try
            A = chol(Cov);
        catch
            [IndTime IndDiscr]
            die
        end
        
        ftmp = zeros(size(msigm));
        xtmp = msigm;
        ftmp(1,:) =  0.05*xtmp(1,:)+0.06*xtmp(2,:)+0.07*xtmp(3,:)+0.08*xtmp(4,:);
        ftmp(2,:) =  0.09*xtmp(1,:)+0.10*xtmp(2,:)+0.11*xtmp(3,:)+0.12*xtmp(4,:);
        ftmp(3,:) =  0.13*xtmp(1,:)+0.14*xtmp(2,:)+0.15*xtmp(3,:)+0.16*xtmp(4,:);
        ftmp(4,:) =  0.17*xtmp(1,:)+0.18*xtmp(2,:)+0.19*xtmp(3,:)+0.20*xtmp(4,:);
            
        [fmean fCov] = UT(ftmp,W);
        
        Q = zeros(L,L);
        Q(1,1) = (Parameters.SigmaRW.Value)^2;
        Q = (Parameters.SigmaRW.Value)^2*eye(L,L);

        dCov = zeros(L,L);
        for k = 1:2*L+1
            tmp = (ftmp(:,k)-fmean)*((msigm(:,k) - xm)');
            dCov = dCov +   W(k)*(tmp+tmp');
        end    
        dCov = dCov +  Q;
        M = A^(-1)*(dCov)*((A^(-1))');
        
        
        phi = zeros(size(M));
        for i = 1:L
            for j = i:L
                if i == j
                    phi(i,j) = 0.5*M(i,j);
                elseif i>j
                    phi(i,j) = M(i,j);
                end
            end
        end
        tmp = A*phi;
        msigm = msigm + (ftmp + (L+kappa)*[zeros(L,1) tmp -tmp])*TStep;
        
        [xm Cov] = UT(msigm,W);
      
%         Cov = (Cov+Cov')/2;
%         
%         [IndTime Cov(3,3)]

        mpred = xm;
        
        plot(IndTime*NbIts + IndDiscr,mpred(plind),'oc')
        plot(IndTime*NbIts + IndDiscr,mpred(plind)+sqrt(Cov(plind,plind)),'+o')
        plot(IndTime*NbIts + IndDiscr,mpred(plind)-sqrt(Cov(plind,plind)),'+o')
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
   
    
    Res.m = mpred;
    Res.msigm = msigm;
    Res.Cov = 1/2*(Cov+Cov');
    Res.Model = Model;
