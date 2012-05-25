function LogRatio = LogRatioGMM(x,xstar,Parameters)

Epsil = Parameters.Epsil;

% numerator (xstar)
temp = 0;
post = posterior(Parameters.SamplDens,x);
for IndComp = 1:Parameters.NbComps
    temp = temp + (mvnpdf(xstar,x + Parameters.Epsil^2/2*(Parameters.Mus(IndComp,:)-x),Parameters.Epsil^2*(squeeze(Parameters.Sigmas(:,:,IndComp))))*post(IndComp));
end
LogqTempStar = log(temp);
LogLikxstar = log(Parameters.f(xstar,Parameters));


% Star to Temp 
temp = 0;
post = posterior(Parameters.SamplDens,xstar);
for IndComp = 1:Parameters.NbComps
    temp = temp + (mvnpdf(x,xstar + Parameters.Epsil^2/2*(Parameters.Mus(IndComp,:)-xstar),Parameters.Epsil^2*(squeeze(Parameters.Sigmas(:,:,IndComp))))*post(IndComp));
end
LogqStarTemp = log(temp);
LogLikx = log(Parameters.f(x,Parameters));

LogNum = LogLikxstar + LogqStarTemp;
LogDenom = LogLikx + LogqTempStar;
LogRatio = LogNum - LogDenom;