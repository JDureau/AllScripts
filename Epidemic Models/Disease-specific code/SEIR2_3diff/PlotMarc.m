function [] = PlotMarc(Res,Resol,figind)

Data = Res.Data;
  
PathsInstant = 0:Res.Parameters.ComputationTStep:sum(Data.NbComputingSteps);
 
try     
    Paths = Res.Paths;
catch
    Paths = Res.CompletePaths;
end        
Parameters = Res.Parameters;

dates = {};
try
    delta = floor((length(Res.Data.Dates)-1)/Resol);
    Resol = floor((length(Res.Data.Dates)-1)/delta);
    inds = delta:delta:(delta)*Resol;
    for i = 1:length(inds)
        dates{i} = [num2str(Res.Data.Dates{inds(i)}.Day) ' ' Res.Data.Dates{inds(i)}.MonthInLetters];
    end
    
catch
    delta = floor((length(Res.Data.Instants)-1)/Resol);
    Resol = floor((length(Res.Data.Instants)-1)/delta);
end

try
    if strcmp(Parameters.FigureMode,'BW') 
             toplot = [5];

             if Parameters.Color
                 Dark  = [100,153,251]/255;
                 Light = [172,215,255]/255;
                 VeryLight = [225 225 225]/255;
                 DotsColor = [34 139 34]/255;
                 SimplotColor = [34 139 34]/255;
             else
                 Dark  = [150 150 150]/255;
                 Light = [189 189 189]/255;
                 VeryLight = [225 225 225]/255;
                 DotsColor = [0 0 0]/255;
                 SimplotColor = [0 0 0]/255;
             end
               
               
             
               
    
               
    %     figure(1)
        try 
            if figind == 1
%                 clf
            end
            rows = 2;
            cols = 2;
        catch
            clf
            rows = 2;
            cols = 1;
            figind = 1;
        end
        i=1;
        nsmooth = Parameters.nsmooth;
        subplot(rows,cols,figind)
        ciplot((quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.025)),(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.975)),Light)
        hold on
        try Parameters.holidays
            ymax = max(mean(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1)))));
            tmpymax = ymax*1.4;
            h = area([7 13],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
            h = area([20 22],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
            set(gca,'Layer','top')
            ciplot((quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.025)),(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.975)),Light)
            ciplot(smooth(quantile(squeeze(exp(Paths(:,6,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,6,:))),0.975),nsmooth),Light)
        end
        
        ciplot((quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.25)),(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.75)),Dark)
        plot((mean(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))))),':k','LineWidth',1.5)
        plot((Data.Observations(5,:)),'.','Color',DotsColor,'MarkerSize',16)
        hold off
        xlim([1 size(Data.Observations,2)-1])
        try
            set(gca,'XTick',[delta:delta:length(Data.Dates)]-1)
    %         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
            set(gca,'XTickLabel',dates)
        catch
            set(gca,'XTick',[0:delta:length(Res.Data.Instants)])
        end
        title('Total Influenza Incidence','FontSize',14)
        ylabel('Incidence','FontSize',12)
        xlabel('t (weeks)','FontSize',12)
        set(gca,'FontSize',12)

        subplot(rows,cols,min(2+figind,rows*cols))
        esttraj = mean(squeeze(exp(Paths(:,6,:))));
        ymax = max(esttraj);
        ymin = min(esttraj);
        ampl = ymax-ymin;
        ylim([ymin-0.7*ampl ymax+0.8*ampl])
  
        ciplot(smooth(quantile(squeeze(exp(Paths(:,6,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,6,:))),0.975),nsmooth),Light)
        hold on
        try Parameters.holidays
            ciplot(smooth(quantile(squeeze(exp(Paths(:,6,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,6,:))),0.75),nsmooth),Dark)
            plot(mean(squeeze(exp(Paths(:,6,:)))),':k','LineWidth',1)
            plot(1000000,ymax,'.','Color',DotsColor,'MarkerSize',16)
            tmpymax = ymax+0.8*ampl;
            h = area([1+Data.Instants(7) 1+Data.Instants(13)],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
            h = area([1+Data.Instants(20) 1+Data.Instants(22)],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
            set(gca,'Layer','top')
            ciplot(smooth(quantile(squeeze(exp(Paths(:,6,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,6,:))),0.975),nsmooth),Light)
        end
        ylim([ymin-0.7*ampl ymax+0.75*ampl])
        ciplot(smooth(quantile(squeeze(exp(Paths(:,6,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,6,:))),0.75),nsmooth),Dark)
        plot(mean(squeeze(exp(Paths(:,6,:)))),':k','LineWidth',1.5)
        
        t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
        
        
        try
            plot(1000000,ymax,'.','Color',DotsColor,'MarkerSize',16)
            if figind == 2
                width = 4;
            else
                width = 1.5;
            end
            plot(Parameters.ComputationTStep:Parameters.ComputationTStep:Data.Instants(end),Data.RealBetaTraj(1:length(Parameters.ComputationTStep:Parameters.ComputationTStep:Data.Instants(end))),'-','Color',SimplotColor,'LineWidth',width)
        end
%         hold off
        
        try
            xlim([0 (size(Data.Observations,2)-1)*7/Parameters.ComputationTStep])
            set(gca,'XTick',[0:delta*7/Parameters.ComputationTStep:Data.Instants(end)])
            set(gca,'XTickLabel',0:delta:52)
        end
        try Parameters.holidays
            xlim([0 (size(Data.Observations,2)-2)*7/Parameters.ComputationTStep])
            set(gca,'XTick',[delta*7/Parameters.ComputationTStep:delta*7/Parameters.ComputationTStep:length(Data.Dates)*21]-2*7/Parameters.ComputationTStep)
    %         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
            set(gca,'XTickLabel',dates)
        end
        
        title('Effective contact rate','FontSize',14)
        ylabel('\beta_t','FontSize',12)
        xlabel('t (weeks)','FontSize',12)
        set(gca,'FontSize',12)
        
%         hold on
%         yis1 = ymax/100:ymax/100:ymax*2/5;
%         yis2 = ymax*3/5:ymax/100:ymax;
%         set(gca,'FontSize',12)
        hold off
  
    else
        die
    end
catch
        
        if size(Paths,2)<10
            toplot = [5];

        %     figure(1)
            clf


            % logged:
            for i = 1:length(toplot)
                subplot(length(toplot)+1,1,i)
                ciplot((quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.025)),(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.975)),[172,215,255]/255)
                hold on
                ciplot((quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.25)),(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.75)),[100,153,251]/255)
                plot((mean(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))))),'k','LineWidth',2)
                plot((Data.Observations(5,:)),'g','LineWidth',2)
                hold off
                xlim([1 size(Data.Observations,2)])
                try
                    set(gca,'XTick',[delta:delta:length(Data.Dates)])
            %         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
                    set(gca,'XTickLabel',dates)
                end
                title('Total Influenza Incidence')
                ylabel('Incidence')
                xlabel('t (weeks)')
            end



            subplot(length(toplot)+1,1,length(toplot)+1)
            nsmooth = 35;
            ciplot(smooth(quantile(squeeze(exp(Paths(:,6,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,6,:))),0.975),nsmooth),[172,215,255]/255)
            hold on
            ciplot(smooth(quantile(squeeze(exp(Paths(:,6,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,6,:))),0.75),nsmooth),[100,153,251]/255)
            plot(mean(squeeze(exp(Paths(:,6,:)))),'k','LineWidth',2)
            t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
            try
                plot(Res.BetatPath,'g','LineWidth',2) 
            end
            hold off
            xlim([0 Data.Instants(end)])
            try
                TicksInds = Data.Instants(inds);
                set(gca,'XTick',TicksInds)
            %     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
            %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
                set(gca,'XTickLabel',dates)
            end
            title('Effective contact rate')
            ylabel('\beta_t')
            xlabel('t (weeks)')
            ymax = 2;
        %     ylim([0.7 ymax])
            hold on
            yis1 = ymax/100:ymax/100:ymax*2/5;
            yis2 = ymax*3/5:ymax/100:ymax;
            hold off
        else
            toplot = [9 10];

        %     figure(1)
            clf


            % logged:
            for i = 1:length(toplot)
                subplot(length(toplot)+1,1,i)
                ciplot((quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.025)),(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.975)),[172,215,255]/255)
                hold on
                ciplot((quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.25)),(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.75)),[100,153,251]/255)
                plot((mean(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))))),'k','LineWidth',2)
                plot((Data.Observations(toplot(i),:)),'g','LineWidth',2)
                hold off
                xlim([1 size(Data.Observations,2)])
                try
                    set(gca,'XTick',[delta:delta:length(Data.Dates)])
            %         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
                    set(gca,'XTickLabel',dates)
                end
                title('Total Influenza Incidence')
                ylabel('Incidence')
                xlabel('t (weeks)')
            end



            subplot(length(toplot)+1,1,length(toplot)+1)
            nsmooth = 35;
            ciplot(smooth(quantile(squeeze(exp(Paths(:,11,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,11,:))),0.975),nsmooth),[172,215,255]/255)
            hold on
            ciplot(smooth(quantile(squeeze(exp(Paths(:,11,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,11,:))),0.75),nsmooth),[100,153,251]/255)
            plot(mean(squeeze(exp(Paths(:,11,:)))),'k','LineWidth',2)
            t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
            try
                plot(Res.BetatPath,'g','LineWidth',2) 
            end
            hold off
            xlim([0 Data.Instants(end)])
            try
                TicksInds = Data.Instants(inds);
                set(gca,'XTick',TicksInds)
            %     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
            %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
                set(gca,'XTickLabel',dates)
            end
            title('Effective contact rate')
            ylabel('\beta_t')
            xlabel('t (weeks)')
            ymax = 2;
        %     ylim([0.7 ymax])
            hold on
            yis1 = ymax/100:ymax/100:ymax*2/5;
            yis2 = ymax*3/5:ymax/100:ymax;
            hold off
            
            
        end
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
