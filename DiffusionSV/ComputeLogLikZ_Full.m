function LogLik = ComputeLogLikZ_Full(Z,Y,Vol,Par)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


H = Par.H.Value;
sigma_X = Par.sigma_X.Value;
mu_Y = Par.mu_Y.Value;
rho = Par.rho.Value;
kappa = Par.kappa.Value;


N = length(Z)/2;
nobs = length(Y);
step = (nobs-1)/N;
npoints = N/(nobs-1);

Bh = Z_to_Bh(Z,N,step,Par);
X = Bh_to_X_Full(Bh,step,Par);


subplot(2,1,1)
plot(X)

subplot(2,1,2)
plot(Vol(X))
pause(0.01)

LogLik = 0;
% mean = sum((mu-Vol(X(1:npoints)).^2/2)*step);
% mean = mean + rho*sum(Vol(X(1:npoints-1)).*Bh(1:npoints-1));
% LogLik = log(normpdf(Y(2), mean ,sqrt(1-rho^2)*sqrt(sum(Vol(X(1:npoints-1)).^2)*step)));   
% ests = zeros(nobs,1);
% ests(2) = sqrt(sum(Vol(X(1:npoints-1)).^2)*step);
for i = 2:nobs
    mean = Y(i-1);
    mean = mean + sum((mu_Y-Vol(X((i-2)*npoints+1:(i-1)*npoints)).^2/2)*step);
    mean = mean + rho*sum(Vol(X((i-2)*npoints+1:(i-1)*npoints)).*Bh((i-2)*npoints+1:(i-1)*npoints));
    LogLik = LogLik + log(normpdf(Y(i),mean,sqrt(1-rho^2)*sqrt(sum(Vol(X((i-2)*npoints+1:(i-1)*npoints)).^2)*step)));   
%     ests(i) = Y(i-1) + sqrt(sum(Vol(X((i-2)*npoints:(i-1)*npoints-1)).^2)*step);
end

Names = Par.Names.Estimated;
if Par.GradCorr
    for j = 1:length(Names)
        LogLik = LogLik ;%- log(Par.(Names{j}).Corr(Names{j},Par));
    end
end




    