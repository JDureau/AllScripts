function xstar = SampleCstGLan(x,Parameters)

PotValue = x;
LogLik = eval(Parameters.LogLikFun);
G = Parameters.G;
sigma = Parameters.RealStd;
mu = Parameters.RealMean;
Epsil = Parameters.Epsil;

der1 = exp(x)/((1+exp(x))^2);
der2 = exp(x)/((1+exp(x)))-mu;
GradLogLik = -1/(sigma^2)*der1*der2;

xstar = x + Epsil^2/2*G^-1*GradLogLik + Epsil*sqrt(G^-1)*randn(1,1);