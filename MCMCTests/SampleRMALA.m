function xstar = SampleRMALA(x,Parameters)
% This is like MALA but the covariance matrix is the local one.



Epsil = Parameters.Epsil;

epsilon = 10^(-10);
fx = log(max(eps,Parameters.f(x,Parameters)));
Grad = zeros(length(x),1);
for i = 1:Parameters.Dim
    xpdx = x;
    xpdx(i) = xpdx(i)+epsilon;
    fxpdx = log(max(eps,Parameters.f(xpdx,Parameters)));
    Grad(i,1) = (fxpdx-fx)/epsilon;
end
G = Parameters.G(x,Parameters);
Gm1 = G^-1;
GDerivs = Parameters.GDerivs(x,Parameters);

mu = x'+Epsil^2/2*Gm1*Grad;

temps = {};
tempmu = zeros(Parameters.Dim,1);
for j = 1:Parameters.Dim
    temps{j} = (Gm1*GDerivs{j}*Gm1);
    for i = 1:Parameters.Dim
        tempmu(i) = tempmu(i) + temps{j}(i,j);
    end 
end
mu = mu - Epsil^2*tempmu;

temps = {};
tempmu = zeros(Parameters.Dim,1);
for j = 1:Parameters.Dim
    temps{j} = trace(Gm1*GDerivs{j});
    for i = 1:Parameters.Dim
        tempmu(i) = tempmu(i) + Gm1(i,j)*temps{j};
    end 
end
mu = mu + Epsil^2/2*tempmu;

xstar = mvnrnd(mu,squeeze(Epsil^2*Gm1));

if isnan(xstar)
    disp('Pb Sample MALA')
end