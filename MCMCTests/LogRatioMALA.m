function LogRatio = LogRatioMALA(x,xstar,Parameters)

Epsil = Parameters.Epsil;
epsilon = 10^(-10);


% numerator (xstar)
Sigma = Parameters.ScalingCov;
% fx = log(Parameters.f(x,Parameters));
% Grad = zeros(length(x),1);
% for i = 1:Parameters.Dim
%     xpdx = x;
%     xpdx(i) = xpdx(i)+epsilon;
%     fxpdx = log(Parameters.f(xpdx,Parameters));
%     Grad(i,1) = (fxpdx-fx)/epsilon;
% end
Grad = ComputeBananaGrad(x,Parameters);
mu = x'+Epsil^2/2*Sigma*Grad;
temp =  mvnpdf(xstar',mu,squeeze(Epsil^2*Sigma));
LogqTempStar = log(temp);
LogLikxstar = log(Parameters.f(xstar,Parameters));


% Star to Temp 
% fx = log(Parameters.f(xstar,Parameters));
% Grad = zeros(length(x),1);
% for i = 1:Parameters.Dim
%     xpdx = xstar;
%     xpdx(i) = xpdx(i)+epsilon;
%     fxpdx = log(Parameters.f(xpdx,Parameters));
%     Grad(i,1) = (fxpdx-fx)/epsilon;
% end
Grad = ComputeBananaGrad(xstar,Parameters);
mu = xstar'+Epsil^2/2*Sigma*Grad;
temp =  mvnpdf(x',mu,squeeze(Epsil^2*Sigma));
LogqStarTemp = log(temp);
LogLikx = log(Parameters.f(x,Parameters));

LogNum = LogLikxstar + LogqStarTemp;
LogDenom = LogLikx + LogqTempStar;
LogRatio = LogNum - LogDenom;