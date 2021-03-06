function [] = PlotForecast(Res,record,go_back)

if nargin == 2
    go_back = 0;
end

Paths = Res.Paths;
Data = Res.Data;
toplot = [7:9];

agregated = zeros(size(Paths,1),size(Paths,2),size(Paths,3)+size(record,3)-go_back);
agregated(:,:,1:size(Paths,3)-go_back) = Paths(:,:,1:size(Paths,3)-go_back);
agregated(:,:,size(Paths,3)-go_back+1:size(Paths,3)+size(record,3)-go_back) = record;


for i = 1:length(toplot)
    ind = toplot(i);
    subplot(length(toplot),1,i)
    ciplot(quantile(squeeze(agregated(:,ind,:)),0.025),quantile(squeeze(agregated(:,ind,:)),0.975),[172,215,255]/255)
    hold on
    ciplot(quantile(squeeze(agregated(:,ind,:)),0.25),quantile(squeeze(agregated(:,ind,:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(agregated(:,ind,:))),'k','LineWidth',2)
    
    
    for j = 2:size(Data.Observations,2)
        if Data.Observations(ind,j)>0
            plot(Data.Instants(j)*ones(1,3),Data.Observations(ind,j)*[0.9 1 1.1],'col',[55,233,30]/255,'LineWidth',3)
        end
    end
    
    if ind == 7
        title('Female Sex Workers','FontWeight','bold')
    elseif ind == 8
        title('Clients','FontWeight','bold')
    elseif ind == 9
        title('Condom use','FontWeight','bold')
    end        
    hold off
    set(gca,'XTick',[0:(600)/5:(600)])
    tmp = size(agregated,3);
    n = ceil(tmp/120);
    xlim([0 n*120])
    set(gca,'XTick',[0:120:n*120])
    Dates = {1985};
    for j = 2:n+1
        Dates{j} = Dates{j-1} + 5;
    end
    set(gca,'XTickLabel',Dates)
end

    