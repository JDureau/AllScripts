function xstar = SampleMALA(x,Parameters)




Epsil = Parameters.Epsil;





Sigma = Parameters.ScalingCov;

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
