function Res = RunJointMCMC_Full(Data,Par)

Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative

Par.GradCorr = 1;
Par.Prior = 1;

Accepted = [];
AcceptedZ = [];
AcceptedPar = [];

if strcmp(Par.theta_sampler,'JointHMC')
    h = Par.h;    
elseif strcmp(Par.theta_sampler,'GibbsHMC')
    hZ = Par.hZ;
    hP = Par.hP;
elseif strcmp(Par.theta_sampler,'GibbsRW')
    hZ = Par.hZ;
    hH = Par.hH;
    hsig = Par.hsig;
    hmu = Par.hmu;
    hrho = Par.hrho;
    hkappa = Par.hkappa;
elseif strcmp(Par.theta_sampler,'Blocks')
    d = Par.d;
end
  

TellParValues(Par)

loop = Par.loop;
nsteps = Par.nsteps;


obsstep = Data.obsstep;
N = Data.N;
Y = Data.Y;
nobs = length(Y);
npoints = N/(nobs-1);

try
    Names = Par.Names.Estimated;
    for i = 1:length(Names)
        Par.(Names{i}).TransfValue = (1+0.00000000001*randn)*Data.ParTrue.(Names{i}).TransfValue;
    end
    Par = TransfToNoTransf(Par);
end



Z = Data.Z; % initialise Z with true path


step = Data.step;
nobs = Data.nobs;

LogLik = ComputeLogLikZ_Full(Z,Y,Vol,Par);
LogPost = LogLik;
LogPriorTheta = ComputeLogPriorZ_Full(Par);
LogPriorZ = - 0.5*(Z')*Z;
% LogPrior = -Inf;

n2 = length(obsstep:obsstep:2*(N-1));
n1 = length(obsstep:obsstep:(N-1));

out_Ls = zeros(loop,1); % Logliks
out_Lposts = zeros(loop,1); % Logliks
out_Zs = zeros(loop,n2);
out_Bhs = zeros(loop,n1);
out_Xs = zeros(loop,n1);
out_Lpriorthetas= zeros(loop,1);
out_LpriorZ= zeros(loop,1);
out_ZpriorGain = zeros(loop,1);
out_LogLikGain = zeros(loop,1);

Par.thetafixed = 0;
Par.Zfixed = 0;
Thetas = zeros(length(Par.Names.Estimated),loop);
TransfThetas = zeros(length(Par.Names.Estimated),loop);

for iter=1:loop %mcmc loop
    disp(iter)

    

    if strcmp(Par.theta_sampler,'JointHMC2')
        % Update Ze : HMC / MALA
        Ve=(randn(1,2*(N)+length(Par.Names.Estimated)))';
        Vestar = Ve;
        Zestar = Z;
        Names = Par.Names.Estimated;
        for i = 1:length(Names)
            Zestar(length(Z) + Par.(Names{i}).Index) = Par.(Names{i}).TransfValue;
        end
        Ze = Zestar;
        ParStar = Par;

        for i = 1:nsteps
            
            Zstar = Zestar(1:end-length(Names));
            for k = 1:length(Names)
                ParStar.(Names{k}).TransfValue = Zestar(length(Zstar) + ParStar.(Names{k}).Index);
            end
            try
                ParStar = TransfToNoTransf(ParStar);
            catch
                'stop';
            end        
            Grad = -ComputeScore_Full(Zstar,Y,Vol,VolDer,ParStar);
            Zestar_h = 4/(4+h^2)*(Zestar + h*Vestar - h^2/4*Zestar - h^2/2*Grad);
            Vestar_hd2 =  Vestar - h/2 * (Zestar + Zestar_h)/2 - h/2 *Grad;
            
            ParStar_h = ParStar;
            Zstar_h = Zestar_h(1:end-length(Names));
            for k = 1:length(Names)
                ParStar_h.(Names{k}).TransfValue = Zestar_h(length(Zstar) + ParStar.(Names{k}).Index);
                if isnan(ParStar_h.(Names{k}).TransfValue)
                    'stop';
                end
            end
            try
                ParStar_h = TransfToNoTransf(ParStar_h);
            catch
                'stop';
            end
            
            Grad= -ComputeScore_Full(Zstar_h,Y,Vol,VolDer,ParStar_h); % gradient
            Vestar_h = Vestar_hd2 - h/2*(Zestar + Zestar_h)/2 - h/2 * Grad;
            Zestar = Zestar_h;
            Vestar = Vestar_h;
            Zstar = Zstar_h;
            ParStar = ParStar_h;

        end

        % Accept / reject Z
        LogLikStar = ComputeLogLikZ_Full(Zstar,Y,Vol,ParStar);
%         LogPriorStar = ComputeLogPriorZ_Full(Zstar,ParStar);

        alpha = LogLikStar  - LogLik   - 0.5*(Zestar')*Zestar + 0.5*(Ze')*Ze - 0.5*(Vestar')*Vestar +0.5*(Ve')*Ve;
        if isnan(alpha)
            disp('nan')
        end
        alpha
%         if Par.GradCorr
%             for k = 1:length(Names)
%                 alpha = alpha + log(Par.(Names{k}).Prior(Names{k},ParStar)) - log(Par.(Names{k}).Prior(Names{k},Par));
%             end
%         end
        
        if ( alpha >= log(rand) )
           Z = Zstar;
           Ve = Vestar;
           LogLikGain = LogLikStar - LogLik;
           ZpriorGain = - 0.5*(Zestar')*Zestar + 0.5*(Ze')*Ze;
           LogLik = LogLikStar;
%            LogPrior = LogPriorStar;
           Par = ParStar;
           Accepted(iter) = 1;
        else
           Accepted(iter) = 0;
        end


        disp(['Acc= ' num2str(mean(Accepted)) '   h=' num2str(h)])
%         h = exp(log(h) + 0.98^iter*(mean(Accepted)-0.7));
% 

    elseif strcmp(Par.theta_sampler,'JointHMC2')
        % Advanced, uniform on theta
        
        % Update Ze : HMC / MALA
        Ve=(randn(1,2*(N)+length(Par.Names.Estimated)))';
        Vestar = Ve;
        Zestar = Z;
        Names = Par.Names.Estimated;
        for i = 1:length(Names)
            Zestar(length(Z) + Par.(Names{i}).Index,1) = Par.(Names{i}).TransfValue;
        end
        Ze = Zestar;
        ParStar = Par;
        Zestar_h = Ze;
        Vestar_hd2 = Ve;
        Vestar_h = Ve;
        
        for i = 1:nsteps
            
            Zstar = Zestar(1:end-length(Names));
            for k = 1:length(Names)
                ParStar.(Names{k}).TransfValue = Zestar(length(Zstar) + ParStar.(Names{k}).Index);
            end
            try
                ParStar = TransfToNoTransf(ParStar);
            catch
                'stop';
            end        
            Grad = -ComputeScore_Full(Zstar,Y,Vol,VolDer,ParStar);
            Zestar_h(1:2*N,1) = 4/(4+h^2)*(Zestar(1:2*N) + h*Vestar(1:2*N) - h^2/4*Zestar(1:2*N) - h^2/2*Grad(1:2*N));
            Vestar_hd2(1:2*N) =  Vestar(1:2*N) - h/2 * (Zestar(1:2*N) + Zestar_h(1:2*N))/2 - h/2 *Grad(1:2*N);
            Zestar_h(2*N+1:end,1) = 400/(400+h^2)*(Zestar(2*N+1:end) + h*Vestar(2*N+1:end) - h^2/400*Zestar(2*N+1:end) - h^2/2*Grad(2*N+1:end));
            Vestar_hd2(2*N+1:end) =  Vestar(2*N+1:end) - h/200 * (Zestar(2*N+1:end) + Zestar_h(2*N+1:end))/2 - h/2 *Grad(2*N+1:end);
            
            
            ParStar_h = ParStar;
            Zstar_h = Zestar_h(1:end-length(Names));
            for k = 1:length(Names)
                ParStar_h.(Names{k}).TransfValue = Zestar_h(length(Zstar) + ParStar.(Names{k}).Index);
                if isnan(ParStar_h.(Names{k}).TransfValue)
                    'stop';
                end
            end
            try
                ParStar_h = TransfToNoTransf(ParStar_h);
            catch
                'stop';
            end
            
            Grad= -ComputeScore_Full(Zstar_h,Y,Vol,VolDer,ParStar_h); % gradient
            Vestar_h(1:2*N) = Vestar_hd2(1:2*N) - h/2*(Zestar(1:2*N) + Zestar_h(1:2*N))/2 - h/2 * Grad(1:2*N);
            Vestar_h(2*N+1:end) = Vestar_hd2(2*N+1:end) - h/200*(Zestar(2*N+1:end) + Zestar_h(2*N+1:end))/2 - h/2 * Grad(2*N+1:end);
            Zestar = Zestar_h;
            Vestar = Vestar_h;
            Zstar = Zstar_h;
            ParStar = ParStar_h;

        end

        % Accept / reject Z
        LogLikStar = ComputeLogLikZ_Full(Zstar,Y,Vol,ParStar);
%         LogPriorStar = ComputeLogPriorZ_Full(Zstar,ParStar);

        alpha = LogLikStar  - LogLik   - 0.5*(Zstar')*Zstar + 0.5*(Z')*Z - 0.5*(Vestar')*Vestar +0.5*(Ve')*Ve - 0.005*(Zestar(2*N+1:end)')*Zestar(2*N+1:end) + 0.005*(Ze(2*N+1:end)')*Ze(2*N+1:end);
        if isnan(alpha)
            disp('nan')
        end
        alpha
%         if Par.GradCorr
%             for k = 1:length(Names)
%                 alpha = alpha + log(Par.(Names{k}).Prior(Names{k},ParStar)) - log(Par.(Names{k}).Prior(Names{k},Par));
%             end
%         end
        
        if ( alpha >= log(rand) )
           Z = Zstar;
           Ve = Vestar;
           LogLik = LogLikStar;
%            LogPrior = LogPriorStar;
           Par = ParStar;
           Accepted(iter) = 1;
        else
           Accepted(iter) = 0;
        end


        disp(['Acc= ' num2str(mean(Accepted)) '   h=' num2str(h)])
%         h = exp(log(h) + 0.98^iter*(mean(Accepted)-0.7));
% 


    elseif strcmp(Par.theta_sampler,'JointHMC')
        
        % advanced on Z, normal on theta
        
        
        % Update Ze : HMC / MALA
        Vz=(randn(1,2*(N)))';
        Vp=(randn(1,length(Par.Names.Estimated)))';
        Vzstar = Vz;
        Vpstar = Vp;
        Zstar = Z;
        Pstar = zeros(length(Par.Names.Estimated),1);
        Names = Par.Names.Estimated;
        for i = 1:length(Names)
            Pstar(Par.(Names{i}).Index) = Par.(Names{i}).TransfValue;
        end
        P = Pstar;
        ParStar = Par;
        Par.Zfixed = 0;
        Par.thetafixed = 0;

        for i = 1:nsteps
            
            for k = 1:length(Names)
                ParStar.(Names{k}).TransfValue = Pstar(ParStar.(Names{k}).Index);
            end
            try
                ParStar = TransfToNoTransf(ParStar);
            catch
                'stop';
            end        
            Grad = -ComputeScore_Full(Zstar,Y,Vol,VolDer,ParStar);
            Zstar_h = 4/(4+h^2)*(Zstar + h*Vzstar - h^2/4*Zstar - h^2/2*Grad(1:length(Zstar)));
            Vzstar_hd2 =  Vzstar - h/2 * (Zstar + Zstar_h)/2 - h/2 *Grad(1:length(Zstar));
            
            Vpstar_hd2 = Vpstar  - h/2*Grad(length(Zstar)+1:end); % the -h/2 Pstar corresponding to the gradient of the prior is included in Grad
            Pstar_h = Pstar + h*Vpstar_hd2;
            ParStar_h = ParStar;
            for k = 1:length(Names)
                ParStar_h.(Names{k}).TransfValue = Pstar_h(ParStar.(Names{k}).Index);
                if isnan(ParStar_h.(Names{k}).TransfValue)
                    'stop';
                end
            end
            try
                ParStar_h = TransfToNoTransf(ParStar_h);
            catch
                'stop';
            end
            
            Grad= -ComputeScore_Full(Zstar_h,Y,Vol,VolDer,ParStar_h); % gradient
            Vpstar_h = Vpstar_hd2 - h/2*Grad(length(Zstar)+1:end); % again prior gradient in Grad;
            Vzstar_h = Vzstar_hd2 - h/2*(Zstar + Zstar_h)/2 - h/2 * Grad(1:length(Zstar));
            Zstar = Zstar_h;
            Vzstar = Vzstar_h;
            Pstar = Pstar_h;
            Vpstar = Vpstar_h;
            ParStar = ParStar_h;

        end

        % Accept / reject Z
        LogLikStar = ComputeLogLikZ_Full(Zstar,Y,Vol,ParStar);
%         LogPriorStar = ComputeLogPriorZ_Full(Zstar,ParStar);

        alpha = LogLikStar  - LogLik   - 0.5*(Zstar')*Zstar + 0.5*(Z')*Z - 0.5*(Vzstar')*Vzstar +0.5*(Vz')*Vz - 0.5*(Vpstar')*Vpstar +0.5*(Vp')*Vp ;
%         alpha =  - 0.5*(Vpstar')*Vpstar +0.5*(Vp')*Vp ;
        if (not(isreal(alpha)))
            disp('nan')
        end
        LogPriorStar = ComputeLogPriorZ_Full(ParStar);
        LogPrior     = ComputeLogPriorZ_Full(Par);
                
        alpha = alpha + LogPriorStar - LogPrior;
        
           
        if ( alpha >= log(rand) )
           Z = Zstar;
           Vz = Vzstar;
           Vp = Vpstar;
           LogLik = LogLikStar;
           LogPriorTheta = LogPriorStar;
           LogPriorZ = - 0.5*(Z')*Z;
%            LogPrior = LogPriorStar;
           Par = ParStar;
           Accepted(iter) = 1;
        else
           Accepted(iter) = 0;
        end
        LogPost = LogLik  + LogPrior - 0.5*(Z')*Z;


        disp(['Acc= ' num2str(mean(Accepted)) '   h=' num2str(h)])
%         h = exp(log(h) + 0.98^iter*(mean(Accepted)-0.7));


    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        % Update Z : HMC / MALA
        Vz=(randn(1,2*(N)))';
        Vzstar = Vz;
        Zstar = Z;
        Par.Zfixed = 0;
        Par.thetafixed = 1;

        for i = 1:nsteps
            
            Grad = -ComputeScore_Full(Zstar,Y,Vol,VolDer,Par);
            Zstar_h = 4/(4+hZ^2)*(Zstar + hZ*Vzstar - hZ^2/4*Zstar - hZ^2/2*Grad(1:length(Zstar)));
            Vzstar_hd2 =  Vzstar - hZ/2 * (Zstar + Zstar_h)/2 - hZ/2 *Grad(1:length(Zstar));
            Grad= -ComputeScore_Full(Zstar_h,Y,Vol,VolDer,Par); % gradient
            Vzstar_h = Vzstar_hd2 - hZ/2*(Zstar + Zstar_h)/2 - hZ/2 * Grad(1:length(Zstar));
            Zstar = Zstar_h;
            Vzstar = Vzstar_h;
          
        end

        % Accept / reject Z
        LogLikStar = ComputeLogLikZ_Full(Zstar,Y,Vol,Par);
%         LogPriorStar = ComputeLogPriorZ_Full(Zstar,ParStar);

        alpha = LogLikStar  - LogLik   - 0.5*(Zstar')*Zstar + 0.5*(Z')*Z - 0.5*(Vzstar')*Vzstar +0.5*(Vz')*Vz ;
        if (not(isreal(alpha)))
            disp('nan')
        end

        
        if ( alpha >= log(rand) )
           Z = Zstar;
           Vz = Vzstar;
           LogLik = LogLikStar;
           LogPriorZ = - 0.5*(Z')*Z;
%            LogPrior = LogPriorStar;
           AcceptedZ(iter) = 1;
        else
           AcceptedZ(iter) = 0;
        end
        
        
        
        
        
        % Update Par : HMC / MALA
        Vp=(randn(1,length(Par.Names.Estimated)))';
        Vpstar = Vp;
        Pstar = zeros(length(Par.Names.Estimated),1);
        Names = Par.Names.Estimated;
        for i = 1:length(Names)
            Pstar(Par.(Names{i}).Index) = Par.(Names{i}).TransfValue;
        end
        P = Pstar;
        ParStar = Par;
        ParStar.Zfixed = 1;
        ParStar.thetafixed = 0;

       
        
        for i = 1:nsteps
            
            for k = 1:length(Names)
                ParStar.(Names{k}).TransfValue = Pstar(ParStar.(Names{k}).Index);
            end
            try
                ParStar = TransfToNoTransf(ParStar);
            catch
                'stop';
            end        
            Grad = -ComputeScore_Full(Z,Y,Vol,VolDer,ParStar);
            Vpstar_hd2 = Vpstar  - hP/2 *Grad;
            Pstar_h = Pstar  + hP*Vpstar_hd2;
            ParStar_h = ParStar;
            for k = 1:length(Names)
                ParStar_h.(Names{k}).TransfValue = Pstar_h(ParStar.(Names{k}).Index);
                if isnan(ParStar_h.(Names{k}).TransfValue)
                    'stop';
                end
            end
            try
                ParStar_h = TransfToNoTransf(ParStar_h);
            catch
                'stop';
            end
            
            Grad= -ComputeScore_Full(Z,Y,Vol,VolDer,ParStar_h); % gradient
            Vpstar_h = Vpstar_hd2 - hP/2*Grad; % again prior gradient in Grad;
            Pstar = Pstar_h;
            Vpstar = Vpstar_h;
            ParStar = ParStar_h;
        end

        % Accept / reject Z
        LogLikStar = ComputeLogLikZ_Full(Z,Y,Vol,ParStar);
%         LogPriorStar = ComputeLogPriorZ_Full(Zstar,ParStar);

        alpha = LogLikStar  - LogLik  - 0.5*(Vpstar')*Vpstar +0.5*(Vp')*Vp ;%- 0.5*(Pstar')*Pstar +0.5*(P')*P;
        if (not(isreal(alpha)))
            disp('nan')
        end
        LogPriorStar = ComputeLogPriorZ_Full(ParStar);
        LogPrior     = ComputeLogPriorZ_Full(Par);

        alpha = alpha + LogPriorStar - LogPrior;
%            
        if ( alpha >= log(rand) )
           Vp = Vpstar;
           LogLik = LogLikStar;
           LogPriorTheta = LogPriorStar;
           LogPriorZ = - 0.5*(Z')*Z;
%            LogPrior = LogPriorStar;
           Par = ParStar;
           AcceptedPar(iter) = 1;
        else
           AcceptedPar(iter) = 0;
        end
        
        
        
        

        disp(['AccZ= ' num2str(mean(AcceptedZ)) '   hZ=' num2str(hZ)])
        disp(['AccPar= ' num2str(mean(AcceptedPar)) '   hPar=' num2str(hP)  ])
%         h = exp(log(h) + 0.98^iter*(mean(Accepted)-0.7));


    
    elseif strcmp(Par.theta_sampler,'GibbsRW')
        % Update Z : HMC / MALA
%         V=(randn(1,2*(N)))';
%         Vstar = V;
%         Zstar = Z;
% 
%         for i = 1:nsteps
%             Grad = -ComputeScore_Full(Zstar,Y,Vol,VolDer,Par);
%             Zstar_h = 4/(4+hZ^2)*(Zstar + hZ*Vstar - hZ^2/4*Zstar - hZ^2/2*Grad);
%             Vstar_hd2 =  Vstar - hZ/2 * (Zstar + Zstar_h)/2 - hZ/2 *Grad;
%             Grad= -ComputeScore_Full(Zstar_h,Y,Vol,VolDer,Par); % gradient
%             Vstar_h = Vstar_hd2 - hZ/2*(Zstar + Zstar_h)/2 - hZ/2 * Grad;
%             Zstar = Zstar_h;
%             Vstar = Vstar_h;
% 
%         end
% 
%         % Accept / reject Z
%         Lstar = ComputeLogLikZ_Full(Zstar,Y,Vol,Par);
%         logPropRatio= -0.5*Zstar'*Zstar - 0.5*Vstar'*Vstar + 0.5*Z'*Z+0.5*V'*V;
% 
%         alpha = Lstar - L + logPropRatio;
%         if ( alpha >= log(rand) )
%            Z = Zstar;
%            V = Vstar;
%            L = Lstar;
%            AcceptedZ(iter) = 1;
%         else
%            AcceptedZ(iter) = 0;
%         end
% 
%     %         disp(['AccZ= ' num2str(mean(AcceptedZ)) '   h=' num2str(hZ)])
%     %         hZ = exp(log(hZ) + 0.98^iter*(mean(AcceptedZ)-0.6));
% 
% 
%         if Par.H.Estimated
%             % Update H Parameters : RW MH
%             StarPar = Par;
%             StarPar.H.TransfValue =  StarPar.H.TransfValue + randn*hH;
%             StarPar = TransfToNoTransf(StarPar);
% 
%             Lstar = ComputeLogLikZ_Full(Z,Y,Vol,StarPar);
%             % unif prior -> nothing
%             % RW -> just derivatives
%             alpha = Lstar - L + log(Par.H.Corr('H',Par)) - log(Par.H.Corr('H',StarPar));
% 
%             if ( alpha >= log(rand) )
%                Par = StarPar;
%                L = Lstar;
%                AcceptedH(iter) = 1;
%             else
%                AcceptedH(iter) = 0;
%             end
%         end
% 
%     %         hH = exp(log(hH) + 0.98^iter*(mean(AcceptedH)-0.23));
% 
%         if Par.sigma_X.Estimated
%          % Update Sig Parameters : RW MH
%             StarPar = Par;
%             StarPar.sigma_X.TransfValue =  StarPar.sigma_X.TransfValue + randn(1,1)*hsig;
%             StarPar = TransfToNoTransf(StarPar);
% 
%             Lstar = ComputeLogLikZ_Full(Z,Y,Vol,StarPar);
%             % unif prior -> nothing
%             % RW -> just derivatives
%             alpha = Lstar - L + log(Par.sigma_X.Corr('sigma_X',Par)) - log(Par.sigma_X.Corr('sigma_X',StarPar));
% 
%             if ( alpha >= log(rand) )
%                Par = StarPar;
%                L = Lstar;
%                Acceptedsig(iter) = 1;
%             else
%                Acceptedsig(iter) = 0;
%             end
%         end
% 
%         disp(['AccH= ' num2str(mean(AcceptedH)) '   h=' num2str(hH) '  ' 'Accsig= ' num2str(mean(Acceptedsig)) '   h=' num2str(hsig) '  ' 'AccZ= ' num2str(mean(AcceptedZ)) '   h=' num2str(hZ)])
%     %         hH = exp(log(hH) + 0.98^iter*(mean(AcceptedH)-0.23));

    elseif  strcmp(Par.theta_sampler,'GibbsHMC')
        % Update Z : HMC / MALA
        Par.Zfixed = 0;
        Par.thetafixed = 1;
        V=(randn(1,2*(N)))';
        Vstar = V;
        Zstar = Z;

        for i = 1:nsteps
            Grad = -ComputeScore_Full(Zstar,Y,Vol,VolDer,Par);
            Zstar_h = 4/(4+hZ^2)*(Zstar + hZ*Vstar - hZ^2/4*Zstar - hZ^2/2*Grad);
            Vstar_hd2 =  Vstar - hZ/2 * (Zstar + Zstar_h)/2 - hZ/2 *Grad;
            Grad= -ComputeScore_Full(Zstar_h,Y,Vol,VolDer,Par); % gradient
            Vstar_h = Vstar_hd2 - hZ/2*(Zstar + Zstar_h)/2 - hZ/2 * Grad;
            Zstar = Zstar_h;
            Vstar = Vstar_h;

        end

        % Accept / reject Z
        LogLikstar = ComputeLogLikZ_Full(Zstar,Y,Vol,Par);
        logPropRatio= -0.5*Zstar'*Zstar - 0.5*Vstar'*Vstar + 0.5*Z'*Z+0.5*V'*V;

        alpha = LogLikstar - LogLik + logPropRatio;
        if 0%( alpha >= log(rand) )
           LogLikGain = LogLikstar - LogLik;
           ZpriorGain = - 0.5*(Zstar')*Zstar + 0.5*(Z')*Z;
           Z = Zstar;
           V = Vstar;
           LogLik = LogLikstar;
           AcceptedZ(iter) = 1;
        else
           AcceptedZ(iter) = 0;
        end

    %         disp(['AccZ= ' num2str(mean(AcceptedZ)) '   h=' num2str(hZ)])
    %         hZ = exp(log(hZ) + 0.98^iter*(mean(AcceptedZ)-0.6));

    
    
    
        % theta given Z
        Par.Zfixed = 1;
        Par.thetafixed = 0;
        Vthetastar=(randn(1,length(Par.Names.Estimated)))';
        thetastar = [];
        Names = Par.Names.Estimated;
        for i = 1:length(Names)
            thetastar(Par.(Names{i}).Index,1) = Par.(Names{i}).TransfValue;
        end
        theta = thetastar;
        Vtheta = Vthetastar;
        ParStar = Par;
        htheta = hP;

        Names = Par.Names.Estimated;
        for i = 1:nsteps
            
            for k = 1:length(Names)
                ParStar.(Names{k}).TransfValue = thetastar(ParStar.(Names{k}).Index);
            end
            try
                ParStar = TransfToNoTransf(ParStar);
            end        
            Grad = -ComputeScore_Full(Z,Y,Vol,VolDer,ParStar)';
            thetastar_h = 4/(4+htheta^2)*(thetastar + htheta*Vthetastar - htheta^2/4*thetastar - htheta^2/2*Grad);
            Vthetastar_hd2 =  Vthetastar - htheta/2 * (thetastar + thetastar_h)/2 - htheta/2 *Grad;
            
            ParStar_h = ParStar;
            for k = 1:length(Names)
                ParStar_h.(Names{k}).TransfValue = thetastar_h(ParStar.(Names{k}).Index);
                if isnan(ParStar_h.(Names{k}).TransfValue)
                    'stop';
                end
            end
            try
                ParStar_h = TransfToNoTransf(ParStar_h);
            catch
                'stop';
            end
            
            Grad= -ComputeScore_Full(Z,Y,Vol,VolDer,ParStar_h)'; % gradient
            Vthetastar_h = Vthetastar_hd2 - htheta/2*(thetastar + thetastar_h)/2 - htheta/2 * Grad;
            thetastar = thetastar_h;
            Vthetastar = Vthetastar_h;
            ParStar = ParStar_h;

        end

        % Accept / reject theta
        LogLikstar = ComputeLogLikZ_Full(Z,Y,Vol,ParStar);
        logPropRatio=  - 0.5*Vthetastar'*Vthetastar + 0.5*Vtheta'*Vtheta;
        
        

        alpha = LogLikstar*0 - LogLik*0 + logPropRatio;
        LogPriorStar = ComputeLogPriorZ_Full(ParStar);
        LogPrior     = ComputeLogPriorZ_Full(Par);

        alpha = alpha + LogPriorStar - LogPrior;
%         if Par.GradCorr
%             for k = 1:length(Names)
%                 alpha = alpha - log(Par.(Names{k}).Corr(Names{k},Par)) + log(Par.(Names{k}).Corr(Names{k},ParStar));
%             end
%         end
        alpha
        if ( alpha >= log(rand) )
           Vtheta = Vthetastar;
           theta = thetastar;
           LogLik = LogLikstar;
           LogPriorTheta = LogPriorStar;
           Par = ParStar;
           Acceptedtheta(iter) = 1;
        else
           Acceptedtheta(iter) = 0;
        end
    
    
    

        disp(['Acctheta= ' num2str(mean(Acceptedtheta)) '   htheta=' num2str(htheta) '  ' 'AccZ= ' num2str(mean(AcceptedZ)) '   hZ=' num2str(hZ) ])
    %         hH = exp(log(hH) + 0.98^iter*(mean(AcceptedH)-0.23));

    
        
        
    elseif  strcmp(Par.theta_sampler,'Blocks')
        % Update Z : HMC / MALA
%         
%         BlockSize = floor(2*(N)/Par.d);
%         
%         for i  = 1:Par.d
%             Zstar = Z;
%             inds = (i-1)*BlockSize+1:min(i*BlockSize,2*(N-1));
%             Zstar(inds) = randn(1,length(inds));
%             
%             Lstar = ComputeLogLikZ(Zstar,Y,Vol,Par);
% 
%             alpha = Lstar - L;
%             if ( alpha >= log(rand) )
%                Z = Zstar;
%                L = Lstar;
%                AcceptedZ(i,iter) = 1;
%             else
%                AcceptedZ(i,iter) = 0;
%             end
%         end
% 
%     
%     
% 
%         disp(['Acc = ' num2str(mean(AcceptedZ'),3) '  d =' num2str(Par.d) ])
%     %         hH = exp(log(hH) + 0.98^iter*(mean(AcceptedH)-0.23));
% 
%     
        
        
    end

    Bh = Z_to_Bh(Z,N,step,Par);
    X = Bh_to_X_Full(Bh,step,Par);
    
    Names = Par.Names.Estimated;
    for k = 1:length(Names)
        Thetas(Par.(Names{k}).Index,iter) = Par.(Names{k}).Value;
        TransfThetas(Par.(Names{k}).Index,iter) = Par.(Names{k}).TransfValue;
    end
    
    out_Lpriorthetas(iter) = LogPriorTheta;
    out_LpriorZ(iter) = LogPriorZ;
    out_Ls(iter) = LogLik ;%+ LogPrior;
    out_Lposts(iter) = LogPost ;
    out_Zs(iter,:) = Z(obsstep:obsstep:2*(N-1));
    out_Bhs(iter,:) = Bh(obsstep:obsstep:N-1);
    out_Xs(iter,:) = X(obsstep:obsstep:N-1);
    try
        out_LogLikGain(iter,:) = LogLikGain ;
        out_ZpriorGain(iter,:) = ZpriorGain ;
    end
end %MCMC loop

TellParValues(Par)


N = (nobs-1)/step;
obsstep = N/(nobs-1);

% Res.toc = toc;
Res.TransfThetas = TransfThetas;
Res.Thetas = Thetas;
Res.Ztrue = Data.Z;
try
    Res.ParTrue = Data.ParTrue;
end
Res.out_Ls = out_Ls;
Res.out_Lposts = out_Lposts;
Res.out_Lpriorthetas = out_Lpriorthetas;
Res.out_LpriorZ = out_LpriorZ;
Res.out_Zs = out_Zs;
Res.out_Xs = out_Xs;
Res.out_LogLikGain = out_LogLikGain;
Res.out_ZpriorGain = out_ZpriorGain;
Res.Z = Z;

if strcmp(Par.theta_sampler,'JointHMC')
    Res.h = h;
elseif strcmp(Par.theta_sampler,'GibbsHMC')
    Res.hZ = hZ;
    Res.hP = hP;    
else
    Res.hZ = hZ;
    Res.hH = hH;
end
Res.Data = Data;
Res.Par = Par; 