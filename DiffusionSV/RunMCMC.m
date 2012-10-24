%Create data and likelihood components objects etc
obstep=T/npoints;%must be integer
step=obstep/(m+1);
nbpoints_discr = npoints*(m+1);
Q=V/sigma; %set R to be the true simulated path (initial run, need to change)
Qb=Q;     %also the code assumes that V(1)=Q(1)=0.
QZ=Q-Q(1);
MU=-kappa*Q; %drift values
dMU=-kappa; %drift derivative
lG=logGirs(MU,dMU,Q(end),kappa,step); %log girsanov value
%lG=logGirs2(Q,MU,step); %log girsanov value
lGb=lG;
lF=zeros(1,npoints); %values of logF for each pair of data
for i=1:npoints
    lF(i)=logF(X(i+1)-X(i),Q((i-1)*(m+1)+1:i*(m+1)),sigma,step);
end
lF1=lF;

%create objects to store MCMC output
out_kappa=zeros(loop,1);
%out_sigma=zeros(loop,1);
out_Q=zeros(min(loop,500),npoints*(m+1)+1);
out_loglike=zeros(loop,1);
out_Qm=zeros(loop,npoints);
midind1=m+2;
midind2=npoints*(m+1)+1;


for iter=1:loop %mcmc loop
    out_kappa(iter)=kappa;
    %out_sigma(iter)=sigma;
    if iter<501
        out_Q(iter,:)=Q;
    end
    out_Qm(iter,:)=Q(midind1:m+1:midind2);    
    out_loglike(iter)=sum(lF)+lG;
    
    %display progress
    if (rem(iter,500) == 0)
        disp(iter)
    end
    
    %update Q : independence sampler - random walk
    if strcmp(Qsampler,'ind_sampler') || strcmp(Qsampler,'rw_diff')
        QZ1 = BrownianMotion(0,T,npoints*(m+1)); % This is part of the proposal from Browian motion
        if strcmp(Qsampler,'rw_diff')
            QZ1 = r*QZ1 + sqrt(1-(r^2))*QZ; 
        end
        Q1=QZ1+Q(1); %for the first point only. the endpoint is QZ
        MU1 = -kappa*Q1;
        lG1 = logGirs(MU1,-kappa,Q1(end),kappa,step);
        %lG1 = logGirs2(Q1,MU1,step);
        for i=1:npoints
            lF1(i)=logF(X(i+1)-X(i),Q1((i-1)*(m+1)+1:i*(m+1)),sigma,step);
        end
        
        alpha = lG1 - lG + sum(lF1) - sum(lF);
        if ( alpha >= log(rand) )
            lG = lG1;
            lF = lF1;
            Q = Q1;
            QZ = QZ1;
        end
    end
    
    if strcmp(Qsampler,'MALA') || strcmp(Qsampler,'HybridMC')
        VDeriv=VarDer(X,Q,kappa,sigma,step,npoints,m);
        u=BrownianMotion(0,T,npoints*(m+1));
        [Q1,QZ1,u1]=ImplicitHybridMC2(Q,QZ,u,h(1),VDeriv,a2,b2,c2,...
            X,kappa,sigma,step,npoints,m,T,nsteps);
        MU1 = -kappa*Q1;
        lG1 = logGirs(MU1,-kappa,Q1(end),kappa,step);
        %lG1 = logGirs2(Q1,MU1,step);
        for i=1:npoints
            lF1(i) =logF(X(i+1)-X(i),Q1((i-1)*(m+1)+1:i*(m+1)),sigma,step);
        end
        logPropRatio=logBMotiondensity(Q1,step)+logBMotiondensity(u1,step)-...
            logBMotiondensity(Q,step)-logBMotiondensity(u,step);
        
        alpha = lG1 - lG + sum(lF1) - sum(lF)+ logPropRatio;
        
        if ( alpha >= log(rand) )
            lG = lG1;
            lF = lF1;
            Q = Q1;
            QZ = QZ1;
        end
    end
    
    %MU=-kappa*Q; %drift values
    %lG=logGirs2(Q,MU,step); %log girsanov value
    
    %update kappa
    if strcmp(theta_sampler,'rw')
        kappa1=kappa*exp(pstd_kappa*randn(1,1));
        MU1 = -kappa1*Q;
        dMU1 = -kappa1;
        lG1 = logGirs(MU1,dMU1,Q(end),kappa1,step);

        alpha =  lG1 - lG;% + log(kappa1) - log(kappa) + kappa - kappa1 ;

        if ( alpha >= log(rand) )
            kappa = kappa1;
            lG = lG1;
        end
    end

    %update kappa Gibbs
    if strcmp(theta_sampler,'Gibbs')
        I=sum(Q.^2)*step;
        %kmu=(T-2-Q(end)^2)/(2*I); %exponetial(1) prior
        kmu=(T-Q(end)^2)/(2*I); %prior propto 1
        kM=((kmu*0.04)+(0.03/I))/(0.04+(1/I)); %trunc normal prior with mean 0.03 and std 0.2 
        %kstd = 1/sqrt(I);
        kstd=sqrt((0.04/I)/((1/I)+0.04));
        %simulate from truncated normal
        kappa=randraw('normaltrunc', [0, 1e5, kM, kstd], 1);
        MU=-kappa*Q; %drift values
        dMU=-kappa; %drift derivative
        lG=logGirs(MU,dMU,Q(end),kappa,step); %log girsanov value
        %lG=logGirs2(Q,MU,step); %log girsanov value
    end
    

%     %update sigma
%     if iter<0
%         sigma1=sigma*exp(pstd_sigma*randn(1,1));
%         for i=1:npoints
%             lF1(i)=logF(X(i+1)-X(i),Q((i-1)*(m+1)+1:i*(m+1)),sigma1,step);
%         end
% 
%         alpha =  sum(lF1) - sum(lF) + log(sigma1) - log(sigma) + sigma - sigma1 ;
% 
%         if ( alpha >= log(rand) )
%             sigma = sigma1;
%             lF = lF1;
%         end
%     end
% 
%     %periodically save output
%     if (rem(iter,10000000) == 0)
%         save(output_file,'out*')
%     end
    
end %MCMC loop

