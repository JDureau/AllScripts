function [] = PlotMarcWithObs(Res,Resol)

Data = Res.Data;
  
PathsInstant = 0:Res.Parameters.ComputationTStep:sum(Data.NbComputingSteps);
 
try     
    Paths = Res.Paths;
catch
    Paths = Res.CompletePaths;
end        
Parameters = Res.Parameters;

dates = {};
delta = floor(length(Res.Data.Dates)/Resol);
Resol = floor(length(Res.Data.Dates)/delta);
inds = delta:delta:delta*Resol;
for i = 1:length(inds)
    dates{i} = [num2str(Res.Data.Dates{inds(i)}.Day) ' ' Res.Data.Dates{inds(i)}.MonthInLetters];
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
        set(gca,'XTick',[delta:delta:length(Data.Dates)])
%         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        set(gca,'XTickLabel',dates)

        title('Estimated Total Influenza Incidence')
    end
    
    % smoothed
    subplot(length(toplot)+2,1,length(toplot)+1)
    temp = diag(Res.Thetas(Parameters.gammam1.Index,:))*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
    ciplot(quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.025),quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.975),[172,215,255]/255)
    hold on
    ciplot(quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.25),quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.75),[100,153,251]/255)
    plot(mean(temp(:,max(1,cumsum(Data.NbComputingSteps)))),'k','LineWidth',2)
    t = 1:length(Data.NbComputingSteps);
    plot(t,1*ones(length(Data.NbComputingSteps),1),'--k','LineWidth',2)
    try
        plot(Res.RtPath,'g','LineWidth',2) 
    end
    hold off
    xlim([1 length(Data.NbComputingSteps)])
    set(gca,'XTick',[delta:delta:length(Data.Dates)])
%     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
%     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
    set(gca,'XTickLabel',dates)
    title('R_t')
    ymax = 3;
    ylim([0 ymax])
    hold on
    yis1 = ymax/100:ymax/100:ymax*2/5;
    yis2 = ymax*3/5:ymax/100:ymax;
    plot(8*ones(size(yis1)),yis1,'r')
    plot(8*ones(size(yis2)),yis2,'r')
    text(8,ymax*1/2,'sch. closure','HorizontalAlignment','center')
    plot(14*ones(size(yis1)),yis1,'r')
    plot(14*ones(size(yis2)),yis2,'r')
    text(14,ymax*1/2,'end of holidays','HorizontalAlignment','center')
    hold off
    
%     subplot(length(toplot)+2,1,length(toplot)+1)
%     temp = diag(Res.Thetas(Parameters.gamma.Index,:).^-1)*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
%     ciplot(quantile(temp,0.025),quantile(temp,0.975),[172,215,255]/255)
%     hold on
%     ciplot(quantile(temp,0.25),quantile(temp,0.75),[100,153,251]/255)
%     plot(mean(temp),'k','LineWidth',2)
%     t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
%     plot(t,1*ones(size(Data.Instants)),'--k','LineWidth',2)
%     try
%         plot(Res.RtPath,'g','LineWidth',2) 
%     end
%     hold off
%     xlim([0 Data.Instants(end)])
%     TicksInds = Data.Instants(inds);
%     set(gca,'XTick',TicksInds)
% %     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
% %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
%     set(gca,'XTickLabel',dates)
%     title('R_t')
%     ymax = 3;
%     ylim([0 ymax])
%     hold on
%     yis1 = ymax/100:ymax/100:ymax*2/5;
%     yis2 = ymax*3/5:ymax/100:ymax;
%     plot(Data.Instants(8)*ones(size(yis1)),yis1,'r')
%     plot(Data.Instants(8)*ones(size(yis2)),yis2,'r')
%     text(Data.Instants(8),ymax*1/2,'sch. closure','HorizontalAlignment','center')
%     plot(Data.Instants(14)*ones(size(yis1)),yis1,'r')
%     plot(Data.Instants(14)*ones(size(yis2)),yis2,'r')
%     text(Data.Instants(14),ymax*1/2,'end of holidays','HorizontalAlignment','center')
%     hold off
%     
%     
% smoothed
    subplot(length(toplot)+2,1,length(toplot)+2)
    temp = squeeze(exp(Paths(:,7,:)));
    ciplot(quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.025),quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.975),[172,215,255]/255)
    hold on
    ciplot(quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.25),quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.75),[100,153,251]/255)
    plot(mean(temp(:,max(1,cumsum(Data.NbComputingSteps)))),'k','LineWidth',2)
    t = 1:length(Data.NbComputingSteps);
    plot(t,1*ones(length(Data.NbComputingSteps),1),'--k','LineWidth',2)

    hold off
    xlim([1 length(Data.NbComputingSteps)])
    set(gca,'XTick',[delta:delta:length(Data.Dates)])
%     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
%     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
    set(gca,'XTickLabel',dates)
    title('\rho_{t}')
    ymax = 0.2;
    ylim([0 ymax])
    hold on
    yis1 = ymax/100:ymax/100:ymax*2/5;
    yis2 = ymax*3/5:ymax/100:ymax;
    plot(8*ones(size(yis1)),yis1,'r')
    plot(8*ones(size(yis2)),yis2,'r')
    text(8,ymax*1/2,'sch. closure','HorizontalAlignment','center')
    plot(14*ones(size(yis1)),yis1,'r')
    plot(14*ones(size(yis2)),yis2,'r')
    text(14,ymax*1/2,'end of holidays','HorizontalAlignment','center')
    hold off

%     subplot(length(toplot)+2,1,length(toplot)+1)
%     ciplot(quantile(squeeze(exp(Paths(:,6,:))),0.025),quantile(squeeze(exp(Paths(:,6,:))),0.975),[172,215,255]/255)
%     hold on
%     ciplot(quantile(squeeze(exp(Paths(:,6,:))),0.25),quantile(squeeze(exp(Paths(:,6,:))),0.75),[100,153,251]/255)
%     plot(mean(squeeze(exp(Paths(:,6,:)))),'k','LineWidth',2)
%     t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
%     try
%         plot(Res.RtPath,'g','LineWidth',2) 
%     end
%     hold off
%     xlim([0 Data.Instants(end)])
%     TicksInds = Data.Instants(inds);
%     set(gca,'XTick',TicksInds)
% %     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
% %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
%     set(gca,'XTickLabel',dates)
%     title('\beta_t')
%     ymax = 3;
%     ylim([0 ymax])
%     hold on
%     yis1 = ymax/100:ymax/100:ymax*2/5;
%     yis2 = ymax*3/5:ymax/100:ymax;
%     plot(Data.Instants(8)*ones(size(yis1)),yis1,'r')
%     plot(Data.Instants(8)*ones(size(yis2)),yis2,'r')
%     text(Data.Instants(8),ymax*1/2,'sch. closure','HorizontalAlignment','center')
%     plot(Data.Instants(14)*ones(size(yis1)),yis1,'r')
%     plot(Data.Instants(14)*ones(size(yis2)),yis2,'r')
%     text(Data.Instants(14),ymax*1/2,'end of holidays','HorizontalAlignment','center')
%     hold off
%     
%     subplot(3,1,3)
%     hist(Res.Thetas(Res.Parameters.RInitProp.Index,:))
%     title('Proportion of initially resistent individuals')
%     colormap('default')
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
