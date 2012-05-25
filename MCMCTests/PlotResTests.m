function [] = PlotResTests(Res,Mu,Sigma)

dim = length(Res.ESS);

figure(1)
for i = 1:dim
    subplot(dim,1,i)
    hist(Res.Vals(i,:));
    title([num2str(Res.ESS(i)/length(Res.Vals(i,:))*100) '%']);
end

figure(2)
for i = 1:dim
    subplot(dim,1,i)
    plot(Res.Vals(i,:))
    if nargin>1
        hold on
        xis = 1:length(Res.Vals(i,:));
        plot(xis,Mu(i)*ones(size(xis)),'g');
        plot(xis,(Mu(i)+sqrt(Sigma(i,i)))*ones(size(xis)),'r');
        plot(xis,(Mu(i)-sqrt(Sigma(i,i)))*ones(size(xis)),'r');
        hold off
    end
end
        
disp(['AccRate: ' num2str(Res.AccRate)])