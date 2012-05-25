function LogRatio = LogRatioMALA(x,xstar,Parameters)

Epsil = Parameters.Epsil;
epsilon = 10^(-10);

% numerator (xstar)
fx = log(max(eps^3,Parameters.f(x,Parameters)));
Grad = zeros(length(x),1);
for i = 1:Parameters.Dim
    xpdx = x;
    xpdx(i) = xpdx(i)+epsilon;
    fxpdx = log(max(eps,Parameters.f(xpdx,Parameters)));
    Grad(i,1) = (fxpdx-fx)/epsilon;
end
temp =  mvnpdf(xstar',x'+Epsil^2/2*Parameters.ScalingCov*Grad,squeeze(Epsil^2*Parameters.ScalingCov));
LogqTempStar = log(temp);
LogLikxstar = log(max(eps^3,Parameters.f(xstar,Parameters)));


% Star to Temp 
fx = log(max(eps^3,Parameters.f(xstar,Parameters)));
Grad = zeros(length(x),1);
for i = 1:Parameters.Dim
    xpdx = xstar;
    xpdx(i) = xpdx(i)+epsilon;
    fxpdx = log(max(eps,Parameters.f(xpdx,Parameters)));
    Grad(i,1) = (fxpdx-fx)/epsilon;
end
temp =  mvnpdf(x',xstar'+Epsil^2/2*Parameters.ScalingCov*Grad,squeeze(Epsil^2*Parameters.ScalingCov));
LogqStarTemp = log(temp);
LogLikx = log(max(eps^3,Parameters.f(x,Parameters)));

LogNum = LogLikxstar + LogqStarTemp;
LogDenom = LogLikx + LogqTempStar;
LogRatio = LogNum - LogDenom;