function xstar = SampleSimpleLan(x,Parameters)

sigma = Parameters.RealStd;
mu = Parameters.RealMean;
Epsil = Parameters.Epsil;

der1 = exp(x)/((1+exp(x))^2);
der2 = exp(x)/((1+exp(x)))-mu;
GradLogLik = -1/(sigma^2)*der1*der2;

xstar = x + Epsil^2/2*GradLogLik + Epsil*randn(1,1);