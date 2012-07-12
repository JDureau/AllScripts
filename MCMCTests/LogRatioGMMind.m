function LogRatio = LogRatioGMMind(x,xstar,Parameters)

Sigmas = Parameters.Sigmas;
Nc = Parameters.Nc;
Dens = Parameters.Dens;

% numerator (xstar)
LogqTempStar = log(pdf(Mixture,xstar));
LogLikxstar = log(Parameters.f(xstar,Parameters));


% Star to Temp 
LogqStarTemp = log(pdf(Mixture,x));
LogLikx = log(Parameters.f(xs,Parameters));

LogNum = LogLikxstar + LogqStarTemp;
LogDenom = LogLikx + LogqTempStar;
LogRatio = LogNum - LogDenom;

if isnan(LogRatio)
    disp('Ratio GMM rand problem')
end