function [] = PlotResGoogle(Res,Resol)

Data = Res.Data;

% Months = ['jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may';'jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'];
% FirstFound = 0;
% SecondFound = 0;
% for i = 1:size(Months,1)
%     if and(not(FirstFound),Months(i,:) == StartMonth)
%         IndFirstMonth = i;
%         FirstFound = 1;
%     end
%     if and(and(FirstFound,not(SecondFound)), Months(i,:) == EndMonth)
%         IndSecondMonth = i;
%         SecondFound = 1;
%     end
% end
% MonthsVect = Months(IndFirstMonth:IndSecondMonth,:);
  
PathsInstant = 0:Res.Parameters.ComputationTStep:sum(Data.NbComputingSteps);
 
try     
    Paths = Res.Paths;
catch
    Paths = Res.CompletePaths;
end        


dates = {};
delta = floor(length(Res.Data.ObservationsDates)/Resol);
inds = delta:delta:length(Data.ObservationsDates);
for i = 1:length(inds)
    dates{i} = [Res.Data.ObservationsDates{inds(i)}(9:10) '/' Res.Data.ObservationsDates{inds(i)}(6:7)];
end

try
    toplot = [5];
    
    figure(1)
    clf
    
    for i = 1:length(toplot)
        subplot(length(toplot)+2,1,i)
        ciplot(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.025),quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.975),[172,215,255]/255)
        hold on
        ciplot(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.25),quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.75),[100,153,251]/255)
        plot(mean(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps))))),'k','LineWidth',2)
        plot(Data.Observations(5,:),'g','LineWidth',2)
        hold off
        xlim([1 size(Data.Observations,2)])
        set(gca,'XTick',[size(Data.Observations,2)/Resol:size(Data.Observations,2)/Resol:size(Data.Observations,2)])
%         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        set(gca,'XTickLabel',dates)

        title('Inlienza-Like Illnesses Incidence')
    end
    
    subplot(length(toplot)+2,1,length(toplot)+1)
    ciplot(quantile(squeeze(exp(Paths(:,6,:))/Res.Parameters.gammam1.Value^-1.*Paths(:,1,:)/Res.Parameters.TotalPopulation),0.025),quantile(squeeze(exp(Paths(:,6,:))/Res.Parameters.gammam1.Value^-1.*Paths(:,1,:)/Res.Parameters.TotalPopulation),0.975),[172,215,255]/255)
    hold on
    ciplot(quantile(squeeze(exp(Paths(:,6,:))/Res.Parameters.gammam1.Value^-1.*Paths(:,1,:)/Res.Parameters.TotalPopulation),0.25),quantile(squeeze(exp(Paths(:,6,:))/Res.Parameters.gammam1.Value^-1.*Paths(:,1,:)/Res.Parameters.TotalPopulation),0.75),[100,153,251]/255)
    plot(mean(squeeze(exp(Paths(:,6,:))/Res.Parameters.gammam1.Value^-1.*Paths(:,1,:)/Res.Parameters.TotalPopulation)),'k','LineWidth',2)
    t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
    plot(t,0.7*ones(size(Data.Instants)),'--k','LineWidth',2)
    plot(t,1.1*ones(size(Data.Instants)),'--k','LineWidth',2)
    try
        plot(Res.RtPath,'g','LineWidth',2) 
    end
    hold off
    xlim([t(1) t(end)])
    set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
%     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
    set(gca,'XTickLabel',dates)
    title('R_t')
    ylim([0 3])
    
    subplot(3,1,3)
    hist(Res.Thetas(Res.Parameters.RInitProp.Index,:))
    hold on
    plot(0.47,0,'og', 'MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
    hold off
    title('Proportion of initially resistent individuals')
    colormap('default')
end
% try
%     figure(2)
%     Names = Res.Parameters.Names.Estimated;
%     k = size(Res.Thetas,1);
%     for i = 1:k
%         subplot(ceil(sqrt(k)),ceil(sqrt(k)),i)
%     %     plot(ResRWML.Thetas(i,:),ResRWML.Paths(:,8,end),'.')
%         hist(Res.Thetas(i,:))
%         hold on
%         title(Names{i})
%     end
% end
