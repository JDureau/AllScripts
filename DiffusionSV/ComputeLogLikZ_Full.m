function [LogLik LogTerm1 LogTerm2 LogTerm3] = ComputeLogLikZ_Full(Z,Obss,Vol,Par)



    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


H = Par.H.Value;
sigma_X = Par.sigma_X.Value;
mu_Y = Par.mu_Y.Value;
rho = Par.rho.Value;
kappa = Par.kappa.Value;
Y = Obss.Y;
Yx = Obss.Yx;

N = length(Z)/2;
nobs = length(Y);
step = (nobs-1)/N/253;
npoints = N/(nobs-1);

Bh = Z_to_Bh(Z,N,step,Par);
X = Bh_to_X_Full(Bh,step,Par);


subplot(2,1,1)
plot(X)

subplot(2,1,2)
plot(Vol(X))
pause(0.01)

vals = [0];
LogLik = 0;
% mean = sum((mu-Vol(X(1:npoints)).^2/2)*step);
% mean = mean + rho*sum(Vol(X(1:npoints-1)).*Bh(1:npoints-1));
% LogLik = log(normpdf(Y(2), mean ,sqrt(1-rho^2)*sqrt(sum(Vol(X(1:npoints-1)).^2)*step)));   
% ests = zeros(nobs,1);
% ests(2) = sqrt(sum(Vol(X(1:npoints-1)).^2)*step);
LogTerm1 = 0;
LogTerm2 = 0;
LogTerm3 = 0;
killit = 0;
for i = 2:nobs
    mean = Y(i-1);
    mean = mean + sum((mu_Y-Vol(X((i-2)*npoints+1:(i-1)*npoints)).^2/2)*step);
    mean = mean + rho*sum(Vol(X((i-2)*npoints+1:(i-1)*npoints)).*Bh((i-2)*npoints+1:(i-1)*npoints));
    LogLik = LogLik + log(normpdf(Y(i),mean,sqrt(1-rho^2)*sqrt(sum(Vol(X((i-2)*npoints+1:(i-1)*npoints)).^2)*step)));   
    %     ests(i) = Y(i-1) + sqrt(sum(Vol(X((i-2)*npoints:(i-1)*npoints-1)).^2)*step);
    vals(i) = LogLik;
    if sqrt(1-rho^2)*sqrt(sum(Vol(X((i-2)*npoints+1:(i-1)*npoints)).^2)*step) < 0.001
        killit = 1;
    end
%     disp(sqrt(1-rho^2)*sqrt(sum(Vol(X((i-2)*npoints+1:(i-1)*npoints)).^2)*step))
    LogTerm1 = LogTerm1 - log(max(0.0000001,sqrt(2*pi)*sqrt(1-rho^2)*sqrt(sum(Vol(X((i-2)*npoints+1:(i-1)*npoints)).^2)*step)));
    LogTerm2 = LogTerm2 - (Y(i)-mean)^2/(2*(1-rho^2)*(sum(Vol(X((i-2)*npoints+1:(i-1)*npoints)).^2)*step));
end
for i = 1:nobs-1
    LogTerm3 = LogTerm3 + log(normpdf(X(i*npoints),Yx(i),Par.tau.Value));
end
if killit
    LogLik = -Inf;
else
    LogLik = LogTerm1 + LogTerm2;
    if Par.obsx 
        LogLik = LogLik + LogTerm3;
    end
end
% Names = Par.Names.Estimated;
% if Par.GradCorr
%     for j = 1:length(Names)
%         LogLik = LogLik + log(Par.(Names{j}).Corr(Names{j},Par));
%     end
% end



    