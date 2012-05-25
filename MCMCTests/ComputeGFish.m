function G = ComputeGFish(x,Parameters)

sigma = Parameters.RealStd;
mu = Parameters.RealMean;

der1 = (2*exp(x)*(1+exp(x))-3*exp(3*x))/((1+exp(x))^4);
der2 = -mu*exp(x)/((1+exp(x))^3);

Der2 = -1/sigma^2*(der1+der2);

G = -Der2;

if G<0
    disp('pb with G')
end