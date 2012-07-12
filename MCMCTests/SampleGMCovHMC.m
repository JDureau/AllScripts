function bigx = SampleGMCovHMC(bigx,Parameters)
% This is like MALA but the covariance matrix is the local one.


Epsil = Parameters.Epsil;
Mixture = Parameters.Dens;
Dim = Parameters.Dim;
x = bigx(1:Dim);

% liks = posterior(Mixture,x);
% if strcmp(Parameters.Mode,'Log')
%     ind = PickRandInd(log(1+liks));
% elseif strcmp(Parameters.Mode,'L')
%     ind = PickRandInd(liks);
% end
% 
liks = posterior(Mixture,x);
[b,maxind] = max(liks);
M = (squeeze(Parameters.Sigmas(:,:,maxind)))^-(1);
Mm1 = M^(-1);
pnp1 = mvnrnd(zeros(Parameters.Dim,1),M)';
p = pnp1;
clf
% vals = zeros(50,2);
% scattercloud(Parameters.TrueSamples(1:500,1),Parameters.TrueSamples(1:500,2));
% hold on
% scattercloudGM(Parameters.TrueSamples(1:500,1),Parameters.TrueSamples(1:500,2),Parameters.Dens);
% hold on
% tmp = [x' x'+p];
% plot(tmp(1,:),tmp(2,:),'r')
for i = 1:100
%     i
    Grad = Parameters.fGrad(x',Parameters);
    tmpp = p + Epsil/2*Grad;
    x = (x' + Epsil*Mm1*tmpp)';
    tmp = [x' x'+tmpp];
%     plot(tmp(1,:),tmp(2,:),'r')

    Grad = Parameters.fGrad(x',Parameters);
    p = tmpp + Epsil/2*Grad;
    
%     hold on
%     plot(x(1),x(2),'.')
% %     disp([x' p Grad])
%     xlim([-10 10])
%     ylim([-10 10])
%     vals(i,:) = x';
%     if rand(1,1)<0.1
%         pause(0.001)
%     end
end
% xlim([min(vals(:,1)) max(vals(:,1))])
% ylim([min(vals(:,2)) max(vals(:,2))])

% hist(vals)
% clf
% scattercloud(vals(:,1),vals(:,2))
% pause(0.01)
% if rand(1,1)<0.02
%     die
% end
xstar = x;
pstar = p;
bigx = [xstar pstar' pnp1'];