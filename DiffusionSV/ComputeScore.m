function Score = ComputeScore(Z,Y,Vol,VolDer,Par)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


%% dL/dZ

H = Par.H.Value;
sigma_X = Par.sigma_X.Value;


N = length(Z)/2+1;
nobs = length(Y);
step = (nobs-1)/N;
npoints = N/(nobs-1);

Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh,sigma_X);


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


%%% b = dlog(Y|B)/dB
b = zeros(N-1,1);
denoms = zeros(1,nobs);
currentk=nobs-1;
denoms(currentk) = sum(Vols_s((currentk-1)*npoints:currentk*npoints-1).^2)*step;
for j = N-1:-1:1
    k = floor(j/npoints)+1;
    if k<currentk
        denoms(k) = sum(Vols_s(max(1,(k-1)*npoints):k*npoints-1).^2)*step;
        currentk=k;
    end
    if j == N-1
        c = 0;
    else
        c = b(j+1);
    end
    b(j) = c - sigma_X * Vols_s_prime(j)*Vols_s(j)*step/denoms(k);
    b(j) = b(j) + (Y(k+1)-Y(k-1+1))^2 * sigma_X * Vols_s_prime(j) * Vols_s(j) * step/(denoms(k)^2);
end

tmp = b';



% P*Delta*M
DiagOfLambda = ComputeDiagOfLambda(N,step,H);


if or(strcmp(Par.theta_sampler,'GibbsRW'),or(strcmp(Par.theta_sampler,'JointHMC'),and(strcmp(Par.theta_sampler,'GibbsHMC'),Par.thetafixed)))
    Score = tmp;
    Score = [Score zeros(1,length(Score))];
    Score = ifft(Score);%,'symmetric');
    Score = sqrt(real(DiagOfLambda)).*Score;
    % Score = MultBySqrtDiag(Lambda,Score);
    Score = MultByMfromRight(N,Score);

    Score = real(Score);
end
%% dL/dsigma

if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.sigma_X.Estimated)
   
    % dL/dX = a
    
    tmp = 0;
    tmp2 = 0;
    for i = 1:length(Bh)
        tmp = tmp +  real(Bh(i));
        tmp2 = tmp2 + a(i)*tmp;
    end
    if Par.GradCorr
        tmp2 = tmp2/(Par.sigma_X.Corr('sigma_X',Par));
    end
    if isnan(tmp2)
        'stop';
    end
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(a)*2+Par.sigma_X.Index) = tmp2;
    elseif strcmp(Par.theta_sampler,'GibbsHMC')
        Score(Par.sigma_X.Index,1) = tmp2;
    end
end

%% dL/dH

if and(or(strcmp(Par.theta_sampler,'JointHMC'),and(Par.Zfixed,strcmp(Par.theta_sampler,'GibbsHMC'))),Par.H.Estimated)

    % dL/dX: a
    
    % (dL/dXi) x (dXi / d\Delta Bj)
    % a * d X / d Delta B
    % a*matrix with down-left half filled with ones 
    tmp = zeros(1,length(a));
    tmp(end) = sigma_X *a(end);
    for i = length(a)-1:-1:1
       tmp(i) = sigma_X *a(i) + tmp(i+1);
    end
    
    % d Delta Bj / d h
    DerOfSqrtDiagLambda = ComputeDerOfSqrtDiagLambda(N,step,H);
    
%     tmp = [tmp zeros(1,length(tmp))];
%     tmp = ifft(tmp);
%     tmp = real(DerOfSqrtDiagLambda.*tmp);
    
    tmp2 = MultByM(N,Z);
    tmp2 = (DerOfSqrtDiagLambda'.*tmp2);
    tmp2 = ifft(tmp2);
    tmp2 = tmp2(1:(N-1));
    
    
    if strcmp(Par.theta_sampler,'JointHMC')
        Score(length(a)*2+Par.H.Index) = real(tmp*tmp2);
        if Par.GradCorr
            Score(length(a)*2+Par.H.Index) = Score(length(a)*2+Par.H.Index)/(Par.H.Corr('H',Par));
        end
        if isnan(Score(length(a)*2+Par.H.Index))
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










    
    
