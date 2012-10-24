function LogLik = ComputeLogLikZ(Z,Y,H,sigma)

N = length(Z)/2;
nobs = length(Y);
step = (nobs-1)/N;
npoints = N/(nobs-1);

Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh);

LogLik = 0;
for i = 1:nobs-1
    LogLik = LogLik + log(normpdf(Y(i+1),0,sqrt(sum(sigma(X((i-1)*npoints+1:i*npoints)).^2))));   
end