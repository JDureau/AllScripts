function bigx = SampleGMCovHMC(bigx,Parameters)
% This is like MALA but the covariance matrix is the local one.


Epsil = Parameters.Epsil;
Mixture = Parameters.Dens;
Dim = Parameters.Dim;
x = bigx(1:Dim);

liks = posterior(Mixture,x);
[b,maxind] = max(liks);
M = squeeze(Parameters.Sigmas(:,:,maxind));
Mm1 = M^(-1);
pnp1 = mvnrnd(zeros(Parameters.Dim,1),M)';
p = pnp1;
clf
% vals = zeros(1000,2);
scattercloud(Parameters.TrueSamples(:,1),Parameters.TrueSamples(:,2));
hold on
for i = 1:20
%     i
    Grad = ComputeGMMGrad(x',Parameters);
    tmpp = p + Epsil/2*Grad;
    x = (x' + Epsil*Mm1*tmpp)';
    Grad = ComputeGMMGrad(x',Parameters);
    p = tmpp + Epsil/2*Grad;
    
    hold on
    plot(x(1),x(2),'.')
%     disp([x' p Grad])
    xlim([-10 10])
    ylim([-10 10])
%     vals(i,:) = x';
%     if rand(1,1)<0.1
%         pause(0.001)
%     end
end
% hist(vals)
% clf
% scattercloud(vals(:,1),vals(:,2))
pause(0.01)
xstar = x;
pstar = p;
bigx = [xstar pstar' pnp1'];