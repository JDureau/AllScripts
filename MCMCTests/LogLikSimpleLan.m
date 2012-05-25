function LogProb = LogLikSimpleLan(x,xstar,Parameters)

sigma = Parameters.RealStd;
mu = Parameters.RealMean;
Epsil = Parameters.Epsil;

der1 = exp(x)/((1+exp(x))^2);
der2 = exp(x)/((1+exp(x)))-mu;
GradLogLik = -1/(sigma^2)*der1*der2;

Correc = log(abs(exp(xstar)/((1+exp(xstar))^2)));
LogLik = log(normpdf(invlogit(xstar),Parameters.RealMean,Parameters.RealStd)/(normcdf(1,Parameters.RealMean,Parameters.RealStd)-normcdf(0,Parameters.RealMean,Parameters.RealStd)));
LogTran = log(normpdf(xstar, x + Epsil^2/2*GradLogLik, Epsil));
LogProb = Correc + LogLik + LogTran;