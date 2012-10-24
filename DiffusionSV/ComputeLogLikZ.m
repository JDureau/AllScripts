function LogLik = ComputeLogLikZ(Z,Y,H,sigma)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


N = length(Z)/2+1;
nobs = length(Y);
step = (nobs-1)/N;
npoints = N/(nobs-1);

Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh);

LogLik = 0;
LogLik = log(normpdf(Y(2),0,sqrt(sum(sigma(X(1:npoints-1)).^2)*step)));   
for i = 3:nobs
    LogLik = LogLik + log(normpdf(Y(i),Y(i-1),sqrt(sum(sigma(X((i-2)*npoints:(i-1)*npoints-1)).^2)*step)));   
end
