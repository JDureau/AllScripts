function xstar = SampleGMCovMALA(x,Parameters)
% This is like MALA but the covariance matrix is the local one.



Epsil = Parameters.Epsil;

Mixture = Parameters.Dens;

liks = posterior(Mixture,x);
[b,maxind] = max(liks);
Sigma = squeeze(Parameters.Sigmas(:,:,maxind));

epsilon = 10^(-10);
fx = log(max(eps,Parameters.f(x,Parameters)));
Grad = zeros(length(x),1);
for i = 1:Parameters.Dim
    xpdx = x;
    xpdx(i) = xpdx(i)+epsilon;
    fxpdx = log(Parameters.f(xpdx,Parameters));
    Grad(i,1) = (fxpdx-fx)/epsilon;
end

mu = x'+Epsil^2/2*Sigma*Grad;
xstar = mvnrnd(mu,squeeze(Epsil^2*Sigma));

if isnan(xstar)
    disp('Pb Sample MALA')
end