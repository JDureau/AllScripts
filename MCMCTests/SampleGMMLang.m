function xstar = SampleGMMLang(x,Parameters)

% In this one I remove the sigma^-1 term in the drift to check if that's
% not better

Mixture = Parameters.Dens;


liks = max(eps,posterior(Mixture,x));
if strcmp(Parameters.Mode,'Log')
    ind = PickRandInd(log(1+liks));
elseif strcmp(Parameters.Mode,'L')
    ind = PickRandInd(liks);
end


mu = Parameters.Mus(ind,:)';
Epsil = Parameters.Epsil;

% clf
% plot(x(1),x(2),'ob', 'MarkerFaceColor','b')
% hold on
% tmp = x'+Epsil^2*(mu-x')/2;
% plot(tmp(1),tmp(2),'ok', 'MarkerFaceColor','k')
% 
% X = [mvnrnd(tmp,squeeze(Epsil^2*Parameters.Dens.Sigma(:,:,ind)),1000)];
% TempDens = gmdistribution.fit(X,1);
% ezcontour(@(x,y)pdf(TempDens,[x y]),[TempDens.mu(1)-3*sqrt(TempDens.Sigma(1,1)) TempDens.mu(1)+3*sqrt(TempDens.Sigma(1,1))],[TempDens.mu(2)-3*sqrt(TempDens.Sigma(2,2)) TempDens.mu(2)+3*sqrt(TempDens.Sigma(2,2))],300);


Sigma = squeeze(Parameters.Sigmas(:,:,ind));

xstar = mvnrnd(x'+Epsil^2*(mu-x')/2,Epsil^2*Sigma);
% plot(xstar(1),xstar(2),'og', 'MarkerFaceColor','g')
% hold off
% xlim([-20 20])
% ylim([-40 20])

% pause(0.1)
