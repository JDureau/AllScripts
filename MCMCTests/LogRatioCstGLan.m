function LogRatio = LogRatioCstGLan(x,xstar,Parameters)

sigma = Parameters.RealStd;
mu = Parameters.RealMean;
Epsil = Parameters.Epsil;
G = Parameters.G;

% numerator (xstar)
der1 = exp(xstar)/((1+exp(xstar))^2);
der2 = exp(xstar)/((1+exp(xstar)))-mu;
GradLogLik = -1/(sigma^2)*der1*der2;

Correc = log(abs(exp(xstar)/((1+exp(xstar))^2)));
LogLik = log(normpdf(invlogit(xstar),Parameters.RealMean,Parameters.RealStd));
LogTran = log(normpdf(x, xstar + Epsil^2/2*G^-1*GradLogLik, Epsil*sqrt(G^-1)));
LogNum = Correc + LogLik + LogTran;

% denominator (x)
der1 = exp(x)/((1+exp(x))^2);
der2 = exp(x)/((1+exp(x)))-mu;
GradLogLik = -1/(sigma^2)*der1*der2;

Correc = log(abs(exp(x)/((1+exp(x))^2)));
LogLik = log(normpdf(invlogit(x),Parameters.RealMean,Parameters.RealStd));
LogTran = log(normpdf(xstar, x + Epsil^2/2*G^-1*GradLogLik, Epsil*sqrt(G^-1)));
LogDenom = Correc + LogLik + LogTran;

LogRatio = LogNum - LogDenom;