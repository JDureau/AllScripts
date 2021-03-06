function xstar = SampleGDyn(x,Parameters)

PotValue = x;
LogLik = eval(Parameters.LogLikFun);

sigma = Parameters.RealStd;
mu = Parameters.RealMean;
der1 = (2*exp(x)*(1+exp(x))-3*exp(3*x))/((1+exp(x))^4);
der2 = -mu*exp(x)/((1+exp(x))^3);
Der2 = -1/sigma^2*(der1+der2);
G = -Der2;


sigma = Parameters.RealStd;
mu = Parameters.RealMean;
Epsil = Parameters.Epsil;

der1 = exp(x)/((1+exp(x))^2);
der2 = exp(x)/((1+exp(x)))-mu;
GradLogLik = -1/(sigma^2)*der1*der2;

xstar = x + Epsil^2/2*GradLogLik + Epsil*sqrt(G^-1)*randn(1,1);