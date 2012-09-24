function bigx = SampleHMC(bigx,Parameters)


Epsil = Parameters.Epsil;
Mixture = Parameters.Dens;
Dim = Parameters.Dim;
x = bigx(1:Dim);

M = Parameters.ScalingCov^(-1);
Mm1 = M^(-1);
% pnp1 = mvnrnd(zeros(Parameters.Dim,1),M)';
pnp1 = mvnrnd(zeros(Parameters.Dim,1),eye(Dim))';
p = pnp1;
% clf
for i = 1:50
    Grad = Parameters.fGrad(x',Parameters);
    tmpp = p + Epsil/2*Grad;
    x = (x' + Epsil*tmpp)';
    Grad = Parameters.fGrad(x',Parameters);
    p = tmpp + Epsil/2*Grad;
end
pause(0.01)
xstar = x;
pstar = p;
bigx = [xstar pstar' pnp1'];