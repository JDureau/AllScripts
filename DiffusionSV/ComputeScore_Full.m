function Score = ComputeScore_Full(Z,Y,Vol,VolDer,Par)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


%% dL/dZ

H = Par.H.Value;
sigma_X = Par.sigma_X.Value;
mu_Y = Par.mu_Y.Value;
rho = Par.rho.Value;
kappa = Par.kappa.Value;
mu_X = Par.mu_X.Value;

N = length(Z)/2;
nobs = length(Y);
step = (nobs-1)/N;
npoints = N/(nobs-1);



ks = zeros(1,N);
for j = 1:N
    ks(j) = ceil((j)/npoints);
end

Bh = Z_to_Bh(Z,N,step,Par);
X = Bh_to_X_Full(Bh,step,Par);

Vols_s = Vol(X);
Vols_s_prime = VolDer(X);


try
    Par.thetafixed;
catch
    Par.thetafixed = 0;
end
try
    Par.Zfixed;
catch
    Par.Zfixed = 0;
end


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
    noms(k)    = (Y(k+1) - Y(k-1+1) - sum(mu_Y - Vols_s((k-1)*npoints+1:k*npoints).^2/2)*step - rho* sum(Vols_s((k-1)*npoints+1:k*npoints)  .* Bh((k-1)*npoints+1:k*npoints)));
end

% denoms(currentk) = (1-rho^2) *sum(Vols_s((currentk-1)*npoints+1:currentk*npoints).^2)*step;
% noms(currentk)    = (Y(currentk+1) - Y(currentk-1+1) - sum(mu - Vols_s((currentk-1)*npoints+1:currentk*npoints).^2/2)*step - rho* sum(Vols_s((currentk-1)*npoints+1:currentk*npoints)  .* Bh((currentk-1)*npoints+1:currentk*npoints)));
% 
% denoms(currentk-1) = (1-rho^2) * sum(Vols_s((currentk-1-1)*npoints+1:(currentk-1)*npoints).^2)*step;
% noms(currentk-1)    = (Y(currentk-1+1) - Y(currentk-1-1+1) - sum(mu - Vols_s((currentk-1-1)*npoints+1:(currentk-1)*npoints).^2/2)*step - rho* sum(Vols_s((currentk-1-1)*npoints+1:(currentk-1)*npoints)  .* Bh((currentk-1-1)*npoints+1:(currentk-1)*npoints)));

Grads = zeros(N,1);
bsigma_X = zeros(1,N);
bkappa = zeros(1,N);
bmu_X = zeros(1,N);
bX0 = ones(1,N);
bZ = zeros(1,N);

bZ(N) = 1;
tmp(N) = rho* Vols_s(N) * noms(currentk)/(denoms(currentk));
nomsfull = zeros(N,1);
denomsfull = zeros(N,1);
nomsfull(1) = noms(1);
denomsfull(1) = denoms(1);
for j = N-1:-1:1
    k = ceil((j+1)/npoints);
    if k<currentk
        currentk=k;
%         if k>1
%             denoms(k-1) = (1-rho^2) * sum(Vols_s((k-1-1)*npoints+1:(k-1)*npoints).^2)*step;
%             noms(k-1)   = (Y(k-1+1) - Y(k-1-1+1) - sum(mu - Vols_s((k-1-1)*npoints+1:(k-1)*npoints).^2/2)*step - rho* sum(Vols_s((k-1-1)*npoints+1:(k-1)*npoints).* Bh((k-1-1)*npoints+1:(k-1)*npoints)));
%         end
    end
    ks(j) = k;
    nomsfull(j+1) = noms(ks(j));
    denomsfull(j+1) = denoms(ks(j));
    c = b(j+1);
    
%     bZ(j ) = (1-kappa*step) * bZ(j+1);
 
%     b(j) = - (1-rho^2) * sigma_X  * Vols_s_prime(j+1)*Vols_s(j+1)*step/denoms(k);

    b(j) = (1-kappa*step) * c  - (1-rho^2) * sigma_X  * Vols_s_prime(j+1)*Vols_s(j+1)*step/denoms(ks(j));
    
    b(j) = b(j) + 1/denoms(ks(j)) * ( - sigma_X  * Vols_s_prime(j+1) * Vols_s(j+1) * step + rho* sigma_X * Vols_s_prime(j+1) * Bh(j+1))*noms(ks(j));
    
    b(j) = b(j) + noms(ks(j))^2 * (1-rho^2) * sigma_X  * Vols_s_prime(j+1) * Vols_s(j+1) * step/(denoms(ks(j))^2);
    
    if ceil((j)/npoints)<ks(j)
        tmp(j) = b(j) + rho* Vols_s(j) * noms(ks(j)-1)/(denoms(ks(j)-1));
    else
        tmp(j) = b(j) + rho* Vols_s(j) * noms(ks(j))/(denoms(ks(j)));
    end
end

Grads = Grads - (1-rho^2) * Vols_s_prime .* Vols_s *step./denomsfull ;
Grads = Grads + (-Vols_s_prime.* Vols_s *step + rho*Vols_s_prime.*Bh)./denomsfull.*nomsfull;
Grads = Grads + (nomsfull.^2).*(1-rho^2) .* Vols_s_prime.* Vols_s*step./(denomsfull.^2);





tmp = tmp';
a = tmp;

% P*Delta*M
DiagOfLambda = ComputeDiagOfLambda(N,step,H);


if 1;%or(strcmp(Par.theta_sampler,'GibbsRW'),or(strcmp(Par.theta_sampler,'JointHMC'),and(strcmp(Par.theta_sampler,'GibbsHMC'),Par.thetafixed)))
    Score = tmp;
    Score = [Score zeros(1,length(Score))];
    Score = ifft(Score);%,'symmetric');
    Score = sqrt(2*N)*sqrt(real(DiagOfLambda)).*Score;
    % Score = MultBySqrtDiag(Lambda,Score);
    Score = MultByMfromRight(N,Score);

    Score = real(Score);
    
    % accunting for the prior
    if Par.Prior
        Score = Score ;%- Z;
    end
end

%%%% Computing b = dX/dpar for par = {sigma_X,kappa,mu_X,X0}
bsigma_X = zeros(1,N);
bkappa = zeros(1,N);
bmu_X = zeros(1,N);
bX0 = ones(1,N);

for j = 1:N
    
    
    
    if j >1
        bsigma_X(j) = (1-kappa*step) * bsigma_X(j-1) + Bh(j-1);
        bkappa(j) = (1-kappa*step) * bkappa(j-1) - X(j-1)*step + mu_X*step;
        bmu_X(j) = (1-kappa*step) * bmu_X(j-1) + kappa*step;
        bX0(j) = (1-kappa*step) * bX0(j-1) ;
    end
    
    
%     Grads(j) = - (1-rho^2) * Vols_s_prime(j)*Vols_s(j)*step/denoms(k);
%     Grads(j) = Grads(j) + 1/denoms(k) * ( - Vols_s_prime(j) * Vols_s(j) * step + rho * Vols_s_prime(j) * Bh(j))*noms(k);
%     Grads(j) = Grads(j) + noms(k)^2 * (1-rho^2) * Vols_s_prime(j) * Vols_s(j) * step/(denoms(k)^2);
    
end

Gradsigma_X = Grads' * bsigma_X';
Gradkappa = Grads' * bkappa';
Gradmu_X = Grads' * bmu_X';
GradX0 = Grads' * bX0';


if Par.sigma_X.Estimated
%     if Par.GradCorr
%         Gradsigma_X = -Gradsigma_X*Par.sigma_X.CorrDer('sigma_X',Par)^2/(Par.sigma_X.Corr('sigma_X',Par)^2);
%     end
%     if isnan(tmp)
%         'stop';
%     end
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(Z)+Par.sigma_X.Index) = Gradsigma_X;
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.sigma_X.Index,1) = Gradsigma_X;
    end
end
if Par.kappa.Estimated
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(Z)+Par.kappa.Index) = Gradkappa;
%         if Par.GradCorr
%             Score(length(Z)+Par.kappa.Index) = Score(length(Z)+Par.kappa.Index)/(Par.kappa.Corr('kappa',Par));
%         end
%         if isnan(Score(length(Z)+Par.kappa.Index))
% %             'stop'
%         end
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.kappa.Index) = Gradkappa;
        if Par.GradCorr
            Score(Par.kappa.Index,1) = Score(Par.kappa.Index)/(Par.kappa.Corr('kappa',Par));
        end
        if isnan(Score(Par.kappa.Index))
%             'stop'
        end
    end
end
if Par.mu_X.Estimated
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(Z)+Par.mu_X.Index) = Gradmu_X;
%         if Par.GradCorr
%             Score(length(Z)+Par.mu_X.Index) = Score(length(Z)+Par.mu_X.Index)/(Par.mu_X.Corr('mu_X',Par));
%         end
%         if isnan(Score(length(Z)+Par.mu_X.Index))
% %             'stop'
%         end
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.mu_X.Index) = Gradmu_X;
        if Par.GradCorr
            Score(Par.mu_X.Index,1) = Score(Par.mu_X.Index)/(Par.mu_X.Corr('mu_X',Par));
        end
        if isnan(Score(Par.mu_X.Index))
%             'stop'
        end
    end
end
if Par.X0.Estimated
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(Z)+Par.X0.Index) = GradX0;
%         if Par.GradCorr
%             Score(length(Z)+Par.X0.Index) = Score(length(Z)+Par.X0.Index)/(Par.X0.Corr('X0',Par));
%         end
%         if isnan(Score(length(Z)+Par.X0.Index))
%             'stop'
%         end
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.X0.Index) = GradX0;
        if Par.GradCorr
            Score(Par.X0.Index,1) = Score(Par.X0.Index)/(Par.X0.Corr('X0',Par));
        end
        if isnan(Score(Par.X0.Index))
            'stop'
        end
    end 
end
    




%% dL/dsigma
% 
% if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.sigma_X.Estimated)
%    
%     %%% tmp = dlog(Y|B)/dsigma_X
% %     denoms = zeros(1,nobs);
% %     noms = zeros(1,nobs);
%     currentk=1;
% 
% %     for k = 1:nobs-1
% %         denoms(k) = (1-rho^2) *sum(Vols_s((k-1)*npoints+1:k*npoints).^2)*step;
% %         noms(k)    = (Y(k+1) - Y(k-1+1) - sum(mu - Vols_s((k-1)*npoints+1:k*npoints).^2/2)*step - rho* sum(Vols_s((k-1)*npoints+1:k*npoints)  .* Bh((k-1)*npoints+1:k*npoints)));
% %     end
% 
% 
%     b = 0; % b = dX_i/d\sigma_X
%     tmp = 0;
%     for j = 1:N
%         k = ceil((j)/npoints);
%         if k>currentk
%             currentk=k;
%         end
% 
%         if j ==1
%             b = 0;
%         else
%             b = (1-kappa*step) * b + Bh(j-1);
%         end
%             
%         tmp = tmp - (1-rho^2) * b * Vols_s_prime(j)*Vols_s(j)*step/denoms(k);
% 
%         tmp = tmp + 1/denoms(k) * ( - b  * Vols_s_prime(j) * Vols_s(j) * step + rho* b * Vols_s_prime(j) * Bh(j))*noms(k);
% 
%         tmp = tmp + noms(k)^2 * (1-rho^2) * b  * Vols_s_prime(j) * Vols_s(j) * step/(denoms(k)^2);
% 
%     end
% 
%     
%     if Par.GradCorr
%         tmp = tmp/(Par.sigma_X.Corr('sigma_X',Par));
%     end
%     if isnan(tmp)
%         'stop';
%     end
%     if strcmp(Par.theta_sampler,'JointHMC')
%         Score(length(Z)+Par.sigma_X.Index) = tmp;
%     elseif strcmp(Par.theta_sampler,'GibbsHMC')
%         Score(Par.sigma_X.Index,1) = tmp;
%     end
% end


%% d L / d k

% if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.kappa.Estimated)
%    
%     %%% tmp = dlog(Y|B)/dk
% %     denoms = zeros(1,nobs);
% %     noms = zeros(1,nobs);
%     currentk=1;
% 
% %     for k = 1:nobs-1
% %         denoms(k) = (1-rho^2) *sum(Vols_s((k-1)*npoints+1:k*npoints).^2)*step;
% %         noms(k)    = (Y(k+1) - Y(k-1+1) - sum(mu - Vols_s((k-1)*npoints+1:k*npoints).^2/2)*step - rho* sum(Vols_s((k-1)*npoints+1:k*npoints)  .* Bh((k-1)*npoints+1:k*npoints)));
% %     end
% 
% 
%     b = 0; % b = dX_i/d\sigma_X
%     tmp = 0;
%     for j = 1:N
%         k = ceil((j)/npoints);
% %         if k>currentk
% %             currentk=k;
% %         end
% 
% %         if j == 1
% %             b = 0;
% %         else
% %             b = (1-kappa*step) * bkappa(j) - X(j-1)*step + mu_X*step;
% %         end
%         
%         tmp = tmp - (1-rho^2) * bkappa(j) * Vols_s_prime(j)*Vols_s(j)*step/denoms(k);
% 
%         tmp = tmp + 1/denoms(k) * ( - bkappa(j)  * Vols_s_prime(j) * Vols_s(j) * step + rho* bkappa(j) * Vols_s_prime(j) * Bh(j))*noms(k);
% 
%         tmp = tmp + noms(k)^2 * (1-rho^2) * bkappa(j)  * Vols_s_prime(j) * Vols_s(j) * step/(denoms(k)^2);
% 
%     end
% 
%     
%     
%     if strcmp(Par.theta_sampler,'JointHMC')
%         Score(length(Z)+Par.kappa.Index) = tmp;
%         if Par.GradCorr
%             Score(length(Z)+Par.kappa.Index) = Score(length(Z)+Par.kappa.Index)/(Par.kappa.Corr('kappa',Par));
%         end
%         if isnan(Score(length(Z)+Par.kappa.Index))
%             'stop'
%         end
%     elseif strcmp(Par.theta_sampler,'GibbsHMC')
%         Score(Par.kappa.Index) = tmp;
%         if Par.GradCorr
%             Score(Par.kappa.Index,1) = Score(Par.kappa.Index)/(Par.kappa.Corr('kappa',Par));
%         end
%         if isnan(Score(Par.kappa.Index))
%             'stop'
%         end
%     end
%    
% end


%% d L / d mu_X

% if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.mu_X.Estimated)
%    
%     %%% tmp = dlog(Y|B)/dk
% %     denoms = zeros(1,nobs);
% %     noms = zeros(1,nobs);
%     currentk=1;
% 
% %     for k = 1:nobs-1
% %         denoms(k) = (1-rho^2) *sum(Vols_s((k-1)*npoints+1:k*npoints).^2)*step;
% %         noms(k)    = (Y(k+1) - Y(k-1+1) - sum(mu - Vols_s((k-1)*npoints+1:k*npoints).^2/2)*step - rho* sum(Vols_s((k-1)*npoints+1:k*npoints)  .* Bh((k-1)*npoints+1:k*npoints)));
% %     end
% 
% 
%     b = 0; % b = dX_i/d\sigma_X
%     tmp = 0;
%     for j = 1:N
%         k = ceil((j)/npoints);
% %         if k>currentk
% %             currentk=k;
% %         end
% 
% %         if j == 1
% %             b = 0;
% %         else
% %             b = (1-kappa*step) * bmu_X(j) + kappa*step;
% %         end
%         
%         tmp = tmp - (1-rho^2) * bmu_X(j) * Vols_s_prime(j)*Vols_s(j)*step/denoms(k);
% 
%         tmp = tmp + 1/denoms(k) * ( - bmu_X(j)  * Vols_s_prime(j) * Vols_s(j) * step + rho* bmu_X(j) * Vols_s_prime(j) * Bh(j))*noms(k);
% 
%         tmp = tmp + noms(k)^2 * (1-rho^2) * bmu_X(j)  * Vols_s_prime(j) * Vols_s(j) * step/(denoms(k)^2);
% 
%     end
% 
%     
%     
%     if strcmp(Par.theta_sampler,'JointHMC')
%         Score(length(Z)+Par.mu_X.Index) = tmp;
%         if Par.GradCorr
%             Score(length(Z)+Par.mu_X.Index) = Score(length(Z)+Par.mu_X.Index)/(Par.mu_X.Corr('mu_X',Par));
%         end
%         if isnan(Score(length(Z)+Par.mu_X.Index))
%             'stop'
%         end
%     elseif strcmp(Par.theta_sampler,'GibbsHMC')
%         Score(Par.mu_X.Index) = tmp;
%         if Par.GradCorr
%             Score(Par.mu_X.Index,1) = Score(Par.mu_X.Index)/(Par.mu_X.Corr('mu_X',Par));
%         end
%         if isnan(Score(Par.mu_X.Index))
%             'stop'
%         end
%     end
%    
% end


%% d L / d X0

% if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.X0.Estimated)
%    
%     %%% tmp = dlog(Y|B)/dk
% %     denoms = zeros(1,nobs);
% %     noms = zeros(1,nobs);
%     currentk=1;
% 
% %     for k = 1:nobs-1
% %         denoms(k) = (1-rho^2) *sum(Vols_s((k-1)*npoints+1:k*npoints).^2)*step;
% %         noms(k)    = (Y(k+1) - Y(k-1+1) - sum(mu - Vols_s((k-1)*npoints+1:k*npoints).^2/2)*step - rho* sum(Vols_s((k-1)*npoints+1:k*npoints)  .* Bh((k-1)*npoints+1:k*npoints)));
% %     end
% 
% 
%     b = 1; % b = dX_i/d\sigma_X
%     tmp = 0;
%     for j = 1:N
%         k = ceil((j)/npoints);
% %         if k>currentk
% %             currentk=k;
% %         end
% %         if j == 1
% %             b = 1;
% %         else
% %             b = (1-kappa*step) * b ;
% %         end
%         
%         tmp = tmp - (1-rho^2) * bX0(j) * Vols_s_prime(j)*Vols_s(j)*step/denoms(k);
% 
%         tmp = tmp + 1/denoms(k) * ( - bX0(j)  * Vols_s_prime(j) * Vols_s(j) * step + rho* bX0(j) * Vols_s_prime(j) * Bh(j))*noms(k);
% 
%         tmp = tmp + noms(k)^2 * (1-rho^2) * bX0(j)  * Vols_s_prime(j) * Vols_s(j) * step/(denoms(k)^2);
% 
%     end
% 
%     
%     
%     if strcmp(Par.theta_sampler,'JointHMC')
%         Score(length(Z)+Par.X0.Index) = tmp;
%         if Par.GradCorr
%             Score(length(Z)+Par.X0.Index) = Score(length(Z)+Par.X0.Index)/(Par.X0.Corr('X0',Par));
%         end
%         if isnan(Score(length(Z)+Par.X0.Index))
%             'stop'
%         end
%     elseif strcmp(Par.theta_sampler,'GibbsHMC')
%         Score(Par.X0.Index) = tmp;
%         if Par.GradCorr
%             Score(Par.X0.Index,1) = Score(Par.X0.Index)/(Par.X0.Corr('X0',Par));
%         end
%         if isnan(Score(Par.X0.Index))
%             'stop'
%         end
%     end
%    
% end





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
            Score(length(Z)+Par.H.Index) = Score(length(Z)+Par.H.Index);%/(Par.H.Corr('H',Par));
        end
%         if isnan(Score(length(Z)+Par.H.Index))
%             'stop';
%         end
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

if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.mu_Y.Estimated)
    tmp = sum(noms(1:nobs-1)./denoms(1:nobs-1));
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(Z)+Par.mu_Y.Index) = tmp;
%         if Par.GradCorr
%             Score(length(Z)+Par.mu_Y.Index) = Score(length(Z)+Par.mu_Y.Index)/(Par.mu_Y.Corr('mu_Y',Par));
%         end
%         if isnan(Score(length(Z)+Par.mu_Y.Index))
%             'stop';
%         end
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.mu.Index) = tmp;
        if Par.GradCorr
            Score(Par.mu_Y.Index,1) = Score(Par.mu_Y.Index)/(Par.mu_Y.Corr('mu_Y',Par));
        end
        if isnan(Score(Par.mu_Y.Index))
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
%         if Par.GradCorr
%             Score(length(Z)+Par.rho.Index) = Score(length(Z)+Par.rho.Index)/(Par.rho.Corr('rho',Par));
%         end
%         if isnan(Score(length(Z)+Par.rho.Index))
%             'stop';
%         end
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.rho.Index) = tmp;
        if Par.GradCorr
            Score(Par.rho.Index,1) = Score(Par.rho.Index)/(Par.rho.Corr('rho',Par));
        end
        if isnan(Score(Par.rho.Index))
            'stop';
        end
    end
end

if Par.Prior
    Names = Par.Names.Estimated; % ATTENTION: this is only for full update. Needs to be adapted for gibbs
    for i = 1:length(Names)
        ind = length(Z)+Par.(Names{i}).Index;
    %     Score(1:ind-1)   = Score(1:ind-1)+log(Par.(Names{i}).Corr(Names{i},Par));
    %     Score(ind+1:end) = Score(ind+1:end)+log(Par.(Names{i}).Corr(Names{i},Par));
        temp = Par.(Names{i}).Prior(Names{i},Par);
        tempDer = Par.(Names{i}).DerPrior(Names{i},Par);
        Score(ind) = Score(ind)*(Par.(Names{i}).Corr(Names{i},Par)) + (temp*Par.(Names{i}).CorrDer(Names{i},Par)+tempDer*Par.(Names{i}).Corr(Names{i},Par))/(temp*Par.(Names{i}).Corr(Names{i},Par));
    end
end


    
    
