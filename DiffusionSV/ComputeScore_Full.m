function Score = ComputeScore_Full(Z,Y,Vol,VolDer,Par)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


%% dL/dZ

H = Par.H.Value;
sigma_X = Par.sigma_X.Value;
mu = Par.mu.Value;
rho = Par.rho.Value;
kappa = Par.kappa.Value;


N = length(Z)/2;
nobs = length(Y);
step = (nobs-1)/N;
npoints = N/(nobs-1);

Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X_Full(Bh,H,step,sigma_X,kappa);


Vols_s = Vol(X);
Vols_s_prime = VolDer(X);

% %  This was from first trials with simple model. After that we had to
% compute directly dlog(Y|B)/dB
% % % a = dlog(Y|X)/dX
% % a = zeros(N-1,1);
% % denoms = zeros(1,nobs);
% % currentk=1;
% % denoms(1) = sum(Vols_s(1:currentk*npoints-1).^2)*step;
% % for j = 1:N-1
% %     k = floor(j/npoints)+1;
% %     if k>currentk
% %         denoms(k) = sum(Vols_s((k-1)*npoints:k*npoints-1).^2)*step;
% %         currentk=k;
% %     end
% %     a(j) = -Vols_s_prime(j)*Vols_s(j)*step/denoms(k);
% %     a(j) = a(j) + (Y(k+1)-Y(k-1+1))^2*Vols_s_prime(j)*Vols_s(j)*step/(denoms(k)^2);
% % end
% % 
% % 
% % 
% % % a * d X / d Delta B
% % % a*matrix with down-left half filled with ones 
% % tmp = zeros(1,length(a));
% % tmp(end) = sigma_X *a(end);
% % for i = length(a)-1:-1:1
% %    tmp(i) = sigma_X *a(i) + tmp(i+1);
% % end


%%% tmp = dlog(Y|B)/dB
b = zeros(N,1);
tmp = zeros(N,1);
denoms = zeros(1,nobs);
noms = zeros(1,nobs);
currentk=nobs-1;

for k = 1:nobs-1
    denoms(k) = (1-rho^2) *sum(Vols_s((k-1)*npoints+1:k*npoints).^2)*step;
    noms(k)    = (Y(k+1) - Y(k-1+1) - sum(mu - Vols_s((k-1)*npoints+1:k*npoints).^2/2)*step - rho* sum(Vols_s((k-1)*npoints+1:k*npoints)  .* Bh((k-1)*npoints+1:k*npoints)));
end

% denoms(currentk) = (1-rho^2) *sum(Vols_s((currentk-1)*npoints+1:currentk*npoints).^2)*step;
% noms(currentk)    = (Y(currentk+1) - Y(currentk-1+1) - sum(mu - Vols_s((currentk-1)*npoints+1:currentk*npoints).^2/2)*step - rho* sum(Vols_s((currentk-1)*npoints+1:currentk*npoints)  .* Bh((currentk-1)*npoints+1:currentk*npoints)));
% 
% denoms(currentk-1) = (1-rho^2) * sum(Vols_s((currentk-1-1)*npoints+1:(currentk-1)*npoints).^2)*step;
% noms(currentk-1)    = (Y(currentk-1+1) - Y(currentk-1-1+1) - sum(mu - Vols_s((currentk-1-1)*npoints+1:(currentk-1)*npoints).^2/2)*step - rho* sum(Vols_s((currentk-1-1)*npoints+1:(currentk-1)*npoints)  .* Bh((currentk-1-1)*npoints+1:(currentk-1)*npoints)));

b(N) = 0;
tmp(N) = rho* Vols_s(N) * noms(currentk)/(denoms(currentk));
for j = N-1:-1:1
    k = ceil((j+1)/npoints);
%     if k<currentk
%         currentk=k;
%         if k>1
%             denoms(k-1) = (1-rho^2) * sum(Vols_s((k-1-1)*npoints+1:(k-1)*npoints).^2)*step;
%             noms(k-1)   = (Y(k-1+1) - Y(k-1-1+1) - sum(mu - Vols_s((k-1-1)*npoints+1:(k-1)*npoints).^2/2)*step - rho* sum(Vols_s((k-1-1)*npoints+1:(k-1)*npoints).* Bh((k-1-1)*npoints+1:(k-1)*npoints)));
%         end
%     end
    
    c = b(j+1);
    
    b(j) = (1+kappa*step) * c - (1-rho^2) * sigma_X  * Vols_s_prime(j+1)*Vols_s(j+1)*step/denoms(k);
    
    b(j) = b(j) + 1/denoms(k) * ( - sigma_X  * Vols_s_prime(j+1) * Vols_s(j+1) * step + rho* sigma_X * Vols_s_prime(j+1) * Bh(j+1))*noms(k);
    
    b(j) = b(j) + noms(k)^2 * (1-rho^2) * sigma_X  * Vols_s_prime(j+1) * Vols_s(j+1) * step/(denoms(k)^2);
    
    if ceil((j)/npoints)<k
        tmp(j) = b(j) + rho* Vols_s(j) * noms(k-1)/(denoms(k-1));
    else
        tmp(j) = b(j) + rho* Vols_s(j) * noms(k)/(denoms(k));
    end
end

tmp = tmp';
a = tmp;

% P*Delta*M
DiagOfLambda = ComputeDiagOfLambda(N,step,H);


if or(strcmp(Par.theta_sampler,'GibbsRW'),or(strcmp(Par.theta_sampler,'JointHMC'),and(strcmp(Par.theta_sampler,'GibbsHMC'),Par.thetafixed)))
    Score = tmp;
    Score = [Score zeros(1,length(Score))];
    Score = ifft(Score);%,'symmetric');
    Score = sqrt(2*N)*sqrt(real(DiagOfLambda)).*Score;
    % Score = MultBySqrtDiag(Lambda,Score);
    Score = MultByMfromRight(N,Score);

    Score = real(Score);
end
%% dL/dsigma

if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.sigma_X.Estimated)
   
    %%% tmp = dlog(Y|B)/dsigma_X
%     denoms = zeros(1,nobs);
%     noms = zeros(1,nobs);
    currentk=1;

%     for k = 1:nobs-1
%         denoms(k) = (1-rho^2) *sum(Vols_s((k-1)*npoints+1:k*npoints).^2)*step;
%         noms(k)    = (Y(k+1) - Y(k-1+1) - sum(mu - Vols_s((k-1)*npoints+1:k*npoints).^2/2)*step - rho* sum(Vols_s((k-1)*npoints+1:k*npoints)  .* Bh((k-1)*npoints+1:k*npoints)));
%     end


    b = 0; % b = dX_i/d\sigma_X
    tmp = 0;
    for j = 2:N
        k = ceil((j)/npoints);
%         if k>currentk
%             currentk=k;
%         end

        b = (1+kappa*step) * b + Bh(j-1);

        tmp = tmp - (1-rho^2) * b * Vols_s_prime(j)*Vols_s(j)*step/denoms(k);

        tmp = tmp + 1/denoms(k) * ( - b  * Vols_s_prime(j) * Vols_s(j) * step + rho* b * Vols_s_prime(j) * Bh(j))*noms(k);

        tmp = tmp + noms(k)^2 * (1-rho^2) * b  * Vols_s_prime(j) * Vols_s(j) * step/(denoms(k)^2);

    end

    
    if Par.GradCorr
        tmp = tmp/(Par.sigma_X.Corr('sigma_X',Par));
    end
    if isnan(tmp)
        'stop';
    end
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(Z)+Par.sigma_X.Index) = tmp;
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.sigma_X.Index,1) = tmp;
    end
end


%% d L / d k

if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.kappa.Estimated)
   
    %%% tmp = dlog(Y|B)/dk
%     denoms = zeros(1,nobs);
%     noms = zeros(1,nobs);
    currentk=1;

%     for k = 1:nobs-1
%         denoms(k) = (1-rho^2) *sum(Vols_s((k-1)*npoints+1:k*npoints).^2)*step;
%         noms(k)    = (Y(k+1) - Y(k-1+1) - sum(mu - Vols_s((k-1)*npoints+1:k*npoints).^2/2)*step - rho* sum(Vols_s((k-1)*npoints+1:k*npoints)  .* Bh((k-1)*npoints+1:k*npoints)));
%     end


    b = 0; % b = dX_i/d\sigma_X
    tmp = 0;
    for j = 2:N
        k = ceil((j)/npoints);
%         if k>currentk
%             currentk=k;
%         end

        b = (1+kappa*step) * b + X(j-1)*step;

        tmp = tmp - (1-rho^2) * b * Vols_s_prime(j)*Vols_s(j)*step/denoms(k);

        tmp = tmp + 1/denoms(k) * ( - b  * Vols_s_prime(j) * Vols_s(j) * step + rho* b * Vols_s_prime(j) * Bh(j))*noms(k);

        tmp = tmp + noms(k)^2 * (1-rho^2) * b  * Vols_s_prime(j) * Vols_s(j) * step/(denoms(k)^2);

    end

    
    
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(Z)+Par.kappa.Index) = tmp;
        if Par.GradCorr
            Score(length(Z)+Par.kappa.Index) = Score(length(Z)+Par.kappa.Index)/(Par.kappa.Corr('kappa',Par));
        end
        if isnan(Score(length(Z)+Par.kappa.Index))
            'stop';
        end
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.kappa.Index) = tmp;
        if Par.GradCorr
            Score(Par.kappa.Index,1) = Score(Par.kappa.Index)/(Par.kappa.Corr('kappa',Par));
        end
        if isnan(Score(Par.kappa.Index))
            'stop';
        end
    end
   
end




%% dL/dH

if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.H.Estimated)

    % dL/dB: a
    
    tmp = a;
    
    % d Delta Bj / d h
    DerOfSqrtDiagLambda = ComputeDerOfSqrtDiagLambda(N,step,H);
    
%     tmp = [tmp zeros(1,length(tmp))];
%     tmp = ifft(tmp);
%     tmp = real(DerOfSqrtDiagLambda.*tmp);
    
    tmp2 = MultByM(N,Z);
    tmp2 = sqrt(2*N)*(DerOfSqrtDiagLambda'.*tmp2);
    tmp2 = ifft(tmp2);
    tmp2 = tmp2(1:(N));
    
    
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(Z)+Par.H.Index) = real(tmp*tmp2);
        if Par.GradCorr
            Score(length(Z)+Par.H.Index) = Score(length(Z)+Par.H.Index)/(Par.H.Corr('H',Par));
        end
        if isnan(Score(length(Z)+Par.H.Index))
            'stop';
        end
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.H.Index) = real(tmp*tmp2);
        if Par.GradCorr
            Score(Par.H.Index,1) = Score(Par.H.Index)/(Par.H.Corr('H',Par));
        end
        if isnan(Score(Par.H.Index))
            'stop';
        end
    end
end




%% dL/dmu

if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.mu.Estimated)
    tmp = sum(noms(1:nobs-1)./denoms(1:nobs-1));
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(Z)+Par.mu.Index) = tmp;
        if Par.GradCorr
            Score(length(Z)+Par.mu.Index) = Score(length(Z)+Par.mu.Index)/(Par.mu.Corr('mu',Par));
        end
        if isnan(Score(length(Z)+Par.mu.Index))
            'stop';
        end
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.mu.Index) = tmp;
        if Par.GradCorr
            Score(Par.mu.Index,1) = Score(Par.mu.Index)/(Par.mu.Corr('mu',Par));
        end
        if isnan(Score(Par.mu.Index))
            'stop';
        end
    end
end


%% dL/d rho

if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.rho.Estimated)
    
    tmp = (nobs-1)*rho/(1-rho^2);
    
    for k = 1:nobs-1
        tmp = tmp + sum(Vols_s((k-1)*npoints+1:(k)*npoints).* Bh((k-1)*npoints+1:(k)*npoints))*noms(k)/denoms(k);
        tmp = tmp - rho*noms(k)^2/((1-rho^2)*denoms(k));
    end
    
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(Z)+Par.rho.Index) = tmp;
        if Par.GradCorr
            Score(length(Z)+Par.rho.Index) = Score(length(Z)+Par.rho.Index)/(Par.rho.Corr('mu',Par));
        end
        if isnan(Score(length(Z)+Par.rho.Index))
            'stop';
        end
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.rho.Index) = tmp;
        if Par.GradCorr
            Score(Par.rho.Index,1) = Score(Par.rho.Index)/(Par.rho.Corr('mu',Par));
        end
        if isnan(Score(Par.rho.Index))
            'stop';
        end
    end
end



    
    
