function [] = scattercloudGM(X1,X2,Density)

% clf
% scattercloud(X1,X2)
cols = rand(Density.NComponents,3);
hold on
for i = 1:Density.NComponents
    thetas = 0:0.01:2*pi;
    [V,D] = eig(squeeze(Density.Sigma(:,:,i)));
    mu = Density.mu(i,:)';
    tmp = repmat(mu,1,length(thetas))  + D(1,1)*repmat(cos(thetas),2,1).*repmat(V(:,1),1,length(thetas)) + D(2,2)*repmat(sin(thetas),2,1).*repmat(V(:,2),1,length(thetas)) ;
    plot(tmp(1,:),tmp(2,:),'col',cols(i,:))
    text(mu(1),mu(2),num2str(Density.PComponents(i)))
end
hold off
