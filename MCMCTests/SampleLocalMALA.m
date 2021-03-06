function xstar = SampleLocalMALA(x,Parameters)
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

xstar = mvnrnd(x'+Epsil^2/2*Gm1*Grad,squeeze(Epsil^2*Gm1));

if isnan(xstar)
    disp('Pb Sample MALA')
end