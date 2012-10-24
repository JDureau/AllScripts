function LogLik = ComputeLogLik(X,Y)

LogLik = 0;
for i = 1:npoints
    LogLik = LogLik + log(normpdf(Y(i),0,sqrt(sum(X((i-1)*npoints:i*npoints-1).^2))));   
end