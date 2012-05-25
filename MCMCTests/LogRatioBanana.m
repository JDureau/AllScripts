function LogRatio = LogRatioBanana(x,xstar,Parameters)



% numerator (xstar)
LogLik = log(fBanana(xstar,Parameters));
LogNum = LogLik;

% denominator (x)
LogLik = log(fBanana(x,Parameters));
LogDenom = LogLik;

LogRatio = LogNum - LogDenom;

