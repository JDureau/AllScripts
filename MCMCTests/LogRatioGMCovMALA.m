function LogRatio = LogRatioRMALA(x,xstar,Parameters)

Epsil = Parameters.Epsil;
epsilon = 10^(-10);

% numerator (xstar)
Mixture = Parameters.Dens;
liks = max(eps,posterior(Mixture,x));
[b,maxind] = max(liks);
mu = Parameters.Mus(maxind,:)';
Sigma = squeeze(Parameters.Sigmas(:,:,maxind));
fx = log(max(eps^3,Parameters.f(x,Parameters)));
Grad = zeros(length(x),1);
for i = 1:Parameters.Dim
    xpdx = x;
    xpdx(i) = xpdx(i)+epsilon;
    fxpdx = log(max(eps,Parameters.f(xpdx,Parameters)));
    Grad(i,1) = (fxpdx-fx)/epsilon;
end
GDerivs = Parameters.GDerivs(x',Parameters);
mu = x'+Epsil^2/2*Sigma*Grad;
% temps = {};
% tempmu = zeros(Parameters.Dim,1);
% for j = 1:Parameters.Dim
%     temps{j} = (Gm1*GDerivs{j}*Gm1);
%     for i = 1:Parameters.Dim
%         tempmu(i) = tempmu(i) + temps{j}(i,j);
%     end 
% end
% mu = mu - Epsil^2*tempmu;
% temps = {};
% tempmu = zeros(Parameters.Dim,1);
% for j = 1:Parameters.Dim
%     temps{j} = trace(Gm1*GDerivs{j});
%     for i = 1:Parameters.Dim
%         tempmu(i) = tempmu(i) + Gm1(i,j)*temps{j};
%     end 
% end
% mu = mu + Epsil^2/2*tempmu;
temp =  mvnpdf(xstar',mu,squeeze(Epsil^2*Sigma));
LogqTempStar = log(temp);
LogLikxstar = log(max(eps^3,Parameters.f(xstar,Parameters)));


% Star to Temp 
Mixture = Parameters.Dens;
liks = max(eps,posterior(Mixture,x));
[b,maxind] = max(liks);
mu = Parameters.Mus(maxind,:)';
Sigma = squeeze(Parameters.Sigmas(:,:,maxind));
fx = log(max(eps^3,Parameters.f(xstar,Parameters)));
Grad = zeros(length(x),1);
for i = 1:Parameters.Dim
    xpdx = xstar;
    xpdx(i) = xpdx(i)+epsilon;
    fxpdx = log(max(eps,Parameters.f(xpdx,Parameters)));
    Grad(i,1) = (fxpdx-fx)/epsilon;
end
% G = Parameters.G(xstar',Parameters);
% GDerivs = Parameters.GDerivs(xstar',Parameters);
mu = xstar'+Epsil^2/2*Sigma*Grad;
% temps = {};
% tempmu = zeros(Parameters.Dim,1);
% for j = 1:Parameters.Dim
%     temps{j} = (Gm1*GDerivs{j}*Gm1);
%     for i = 1:Parameters.Dim
%         tempmu(i) = tempmu(i) + temps{j}(i,j);
%     end 
% end
% mu = mu - Epsil^2*tempmu;
% temps = {};
% tempmu = zeros(Parameters.Dim,1);
% for j = 1:Parameters.Dim
%     temps{j} = trace(Gm1*GDerivs{j});
%     for i = 1:Parameters.Dim
%         tempmu(i) = tempmu(i) + Gm1(i,j)*temps{j};
%     end 
% end
% mu = mu + Epsil^2/2*tempmu;
temp =  mvnpdf(x',mu,squeeze(Epsil^2*Sigma));
LogqStarTemp = log(temp);
LogLikx = log(max(eps^3,Parameters.f(x,Parameters)));

LogNum = LogLikxstar + LogqStarTemp;
LogDenom = LogLikx + LogqTempStar;
LogRatio = LogNum - LogDenom;