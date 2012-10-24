function LogLik = ComputeLogLikX(X,Y,H,sigma)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


N = length(X)+1;
nobs = length(Y);
step = (nobs-1)/N;
npoints = N/(nobs-1);

LogLik = 0;
LogLik = log(normpdf(Y(2),0,sqrt(sum(sigma(X(1:npoints-1)).^2)*step)));   
for i = 3:nobs
    LogLik = LogLik + log(normpdf(Y(i),Y(i-1),sqrt(sum(sigma(X((i-2)*npoints:(i-1)*npoints-1)).^2)*step)));   
end