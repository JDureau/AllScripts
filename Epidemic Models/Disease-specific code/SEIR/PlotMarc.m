function [] = PlotMarc(Res,Resol,figind)

try
    Data = Res.Data;
end
try
    PathsInstant = 0:Res.Parameters.ComputationTStep:sum(Data.NbComputingSteps);
end
try     
    Paths = Res.Paths;
catch
    Paths = Res.CompletePaths;
end        
try
    Parameters = Res.Parameters;
end
dates = {};
try
    delta = floor((length(Res.Data.Dates)-1)/Resol);
    Resol = floor((length(Res.Data.Dates)-1)/delta);
    inds = delta:delta:(delta)*Resol;
    for i = 1:length(inds)
        dates{i} = [num2str(Res.Data.Dates{inds(i)}.Day) ' ' Res.Data.Dates{inds(i)}.MonthInLetters];
    end
    
catch
    try
        delta = floor((length(Res.Data.Instants)-1)/Resol);
        Resol = floor((length(Res.Data.Instants)-1)/delta);
    end
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
        plot((median(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))))),':k','LineWidth',1.5)
        plot((Data.Observations(5,:)),'.','Color',DotsColor,'MarkerSize',26)
        hold off
        xlim([1 size(Data.Observations,2)-1])
        try
            set(gca,'XTick',[delta:delta:length(Data.Dates)]-1)
    %         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
            set(gca,'XTickLabel',dates)
        catch
            set(gca,'XTick',[0:delta:length(Res.Data.Instants)])
        end
        title('Total Influenza Incidence','FontSize',26)
        ylabel('Incidence','FontSize',22)
%         xlabel('t (weeks)','FontSize',22)
        set(gca,'FontSize',22)

        subplot(rows,cols,min(2+figind,rows*cols))
        esttraj = median(squeeze(exp(Paths(:,6,:))));
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
        plot(median(squeeze(exp(Paths(:,6,:)))),':k','LineWidth',1.5)
        
        t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
        
        
        try
%             plot(1000000,ymax,'.','Color',DotsColor,'MarkerSize',16)
            if figind == 2
                width = 4;
            else
                width = 1.1;
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
        
        title('Effective contact rate','FontSize',26)
        ylabel('\beta_t','FontSize',22)
%         xlabel('t (weeks)','FontSize',22)
        set(gca,'FontSize',22)
        
%         hold on
%         yis1 = ymax/100:ymax/100:ymax*2/5;
%         yis2 = ymax*3/5:ymax/100:ymax;
%         set(gca,'FontSize',12)
        hold off
  
    else
        die
    end
catch
       
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
             
        if size(Res.Paths,2) <8
            toplot = [5];
        else
            toplot = [9 10];
        end
    %     figure(1)
        clf



        for i = 1:length(toplot)
            if size(Res.Paths,2) <8
                subplot(length(toplot)+1,1,i)
            else
                subplot(length(toplot)+2,2,2*(i-1)+1:2*(i-1)+2)
            end
            ciplot((quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.025)),(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.975)),[172,215,255]/255)
            hold on
            try Parameters.holidays
                ymax = max(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.975));
                tmpymax = 1000;%ymax*1.4;
                h = area([7 13],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
                h = area([20 22],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
                set(gca,'Layer','top')
                ciplot((quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.025)),(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.975)),[172,215,255]/255)
                ciplot((quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.25)),(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.75)),[100,153,251]/255)
                ylim([0 tmpymax])
            end            
            ciplot((quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.25)),(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.75)),[100,153,251]/255)
            plot((mean(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))))),':k','LineWidth',1)
            plot((Data.Observations(toplot(i),:)),'.','Color',DotsColor,'MarkerSize',10)
            hold off
            xlim([1 size(Data.Observations,2)-1])
            try
                set(gca,'XTick',[delta:delta:length(Data.Dates)]-1)
        %         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
                set(gca,'XTickLabel',dates)
            catch
                set(gca,'XTick',[0:delta:length(Res.Data.Instants)])
            end
%             try
%                 set(gca,'XTick',[delta:delta:length(Data.Dates)])
%         %         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
%                 set(gca,'XTickLabel',dates)
%             end
            if size(Res.Paths,2) <8
                title('Total Influenza Incidence')
            else
                if i == 1
                    title('Influenza incidence among kids')
                else
                    title('Influenza incidence among adults')
                end
            end
            ylabel('Incidence')
%             xlabel('t (weeks)')
        end

    %     % smoothed
    %     subplot(length(toplot)+1,1,length(toplot)+1)
    %     temp = diag(Res.Thetas(Parameters.gamma.Index,:).^-1)*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
    %     ciplot(quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.025),quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.975),[172,215,255]/255)
    %     hold on
    %     ciplot(quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.25),quantile(temp(:,max(1,cumsum(Data.NbComputingSteps))),0.75),[100,153,251]/255)
    %     plot(mean(temp(:,max(1,cumsum(Data.NbComputingSteps)))),'k','LineWidth',2)
    %     t = 1:length(Data.NbComputingSteps);
    %     plot(t,1*ones(length(Data.NbComputingSteps),1),'--k','LineWidth',2)
    %     try
    %         plot(Res.RtPath,'g','LineWidth',2) 
    %     end
    %     hold off
    %     xlim([1 length(Data.NbComputingSteps)])
    %     set(gca,'XTick',[delta:delta:length(Data.Dates)])
    % %     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
    % %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
    %     set(gca,'XTickLabel',dates)
    %     title('R_t')
    %     ymax = 3;
    %     ylim([0 ymax])
    %     hold on
    %     yis1 = ymax/100:ymax/100:ymax*2/5;
    %     yis2 = ymax*3/5:ymax/100:ymax;
    %     plot(8*ones(size(yis1)),yis1,'r')
    %     plot(8*ones(size(yis2)),yis2,'r')
    %     text(8,ymax*1/2,'sch. closure','HorizontalAlignment','center')
    %     plot(14*ones(size(yis1)),yis1,'r')
    %     plot(14*ones(size(yis2)),yis2,'r')
    %     text(14,ymax*1/2,'end of holidays','HorizontalAlignment','center')
    %     hold off

    %     subplot(length(toplot)+2,1,length(toplot)+1)
    %     try
    %         temp = diag(Res.Thetas(Parameters.gamma.Index,:).^-1)*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
    %     catch
    %         try
    %             temp = diag(Res.Thetas(Parameters.gammam1.Index,:))*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
    %         catch
    %             temp = Parameters.gammam1.Value*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
    %         end
    %     end
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

        if size(Res.Paths,2) <8
            toplotdiff = [6];

            subplot(length(toplot)+1,1,length(toplot)+1)
            nsmooth = 35;
            ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
            hold on
            try Parameters.holidays
                ymax = max(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.975));
                tmpymax = 3;%ymax*1.4;
                h = area([7 13],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
                h = area([20 22],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
                set(gca,'Layer','top')
                ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
                ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
                ylim([0 tmpymax])
            end  
            ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
            plot(median(squeeze(exp(Paths(:,toplotdiff ,:)))),':k','LineWidth',1)
            title('Effective contact rate')
            ylabel('\beta_t')
            xlabel('t (weeks)')
        else
            toplotdiff = [11];
            
            subplot(length(toplot)+2,2,length(toplot)*2+1:length(toplot)*2+2)
            nsmooth = 5;
            ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
            hold on
%             try Parameters.holidays
%                 ymax = max(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.975));
%                 tmpymax = 3;%ymax*1.4;
%                 h = area([1+Data.Instants(7) 1+Data.Instants(13)],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
%                 h = area([1+Data.Instants(20) 1+Data.Instants(22)],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
%                 set(gca,'Layer','top')
%                 ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
%                 ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
%                 ylim([0 tmpymax])
%             end  
            ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
            plot(mean(squeeze(exp(Paths(:,toplotdiff ,:)))),':k','LineWidth',1)
            title('Effective contacts from kids to kids')
            ylabel('\beta_{kids \rightarrow kids}')
%             xlabel('t (weeks)')
%             ylim([0 3])

            xlim([0 Data.Instants(end)])
            try Parameters.holidays
                xlim([0 (size(Data.Observations,2)-2)*7/Parameters.ComputationTStep])
                set(gca,'XTick',[delta*7/Parameters.ComputationTStep:delta*7/Parameters.ComputationTStep:length(Data.Dates)*21]-2*7/Parameters.ComputationTStep)
                %         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
                set(gca,'XTickLabel',dates)
            end
 
%             try
%                 TicksInds = Data.Instants(inds);
%                 set(gca,'XTick',TicksInds)
%             %     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
%             %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
%                 set(gca,'XTickLabel',dates)
%             end
            
%             subplot(length(toplot)+2,2,length(toplot)*2+2)
%             nsmooth = 35;
%             if strcmp(Parameters.Problem,'Marc3diff')
%                 toplotdiff = 12;
%                 ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
%                 hold on
%                 ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
%                 plot(mean(squeeze(exp(Paths(:,toplotdiff ,:)))),'k','LineWidth',2)
%             elseif strcmp(Parameters.Problem,'Marc2diff')
%                 mult = Res.Thetas(Parameters.adultsmult.Index,:);
%                 ciplot(repmat(quantile(mult.*squeeze(Res.Thetas(Parameters.beta12init.Index,:)),0.025),1,size(Paths,3)),repmat(quantile(mult.*squeeze(Res.Thetas(Parameters.beta12init.Index,:)),0.975),1,size(Paths,3)),[172,215,255]/255)
%                 hold on
%                 ciplot(repmat(quantile(mult.*squeeze(Res.Thetas(Parameters.beta12init.Index,:)),0.25),1,size(Paths,3)),repmat(quantile(mult.*squeeze(Res.Thetas(Parameters.beta12init.Index,:)),0.75),1,size(Paths,3)),[100,153,251]/255)
%                 plot(repmat(median(mult.*squeeze(Res.Thetas(Parameters.beta12init.Index,:))),1,size(Paths,3)),'k','LineWidth',2)
%             elseif strcmp(Parameters.Problem,'Marc2diffb')
%                 toplotdiff = 11;
%                 add = Res.Thetas(Parameters.kidsadd.Index,:);
%                 mult = Res.Thetas(Parameters.kidsmult.Index,:);
%                 ciplot(smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
%                 hold on
%                 ciplot(smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
%                 plot(mean(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff ,:)))),'k','LineWidth',2)
%             else
%                 add = Res.Thetas(Parameters.kidsadd.Index,:);
%                 mult = Res.Thetas(Parameters.kidsmult.Index,:);
%                 ciplot(smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
%                 hold on
%                 ciplot(smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
%                 plot(mean(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff ,:)))),'k','LineWidth',2)
%             end
%             title('adults to kids')
%             ylabel('\beta_{21}(t)')
%             xlabel('t (weeks)')
            
%             subplot(length(toplot)+1,2,length(toplot)*2+2)
%             nsmooth = 35;
%             if strcmp(Parameters.Problem,'Marc3diff')
%                 toplotdiff = 11;
%                 mult = Res.Thetas(Parameters.adultsmult.Index,:);
%                 ciplot(smooth(quantile(repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
%                 hold on
%                 ciplot(smooth(quantile(repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
%                 plot(mean(repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff ,:)))),'k','LineWidth',2)
%             elseif strcmp(Parameters.Problem,'Marc2diff')
%                 ciplot(repmat(quantile(squeeze(Res.Thetas(Parameters.beta12init.Index,:)),0.025),1,size(Paths,3)),repmat(quantile(squeeze(Res.Thetas(Parameters.beta12init.Index,:)),0.975),1,size(Paths,3)),[172,215,255]/255)
%                 hold on
%                 ciplot(repmat(quantile(squeeze(Res.Thetas(Parameters.beta12init.Index,:)),0.25),1,size(Paths,3)),repmat(quantile(squeeze(Res.Thetas(Parameters.beta12init.Index,:)),0.75),1,size(Paths,3)),[100,153,251]/255)
%                 plot(repmat(median(squeeze(Res.Thetas(Parameters.beta12init.Index,:))),1,size(Paths,3)),'k','LineWidth',2)
%             elseif strcmp(Parameters.Problem,'Marc2diffb')
%                 toplotdiff = 11;
%                 add = Res.Thetas(Parameters.adultsadd.Index,:);
%                 mult = Res.Thetas(Parameters.adultsmult.Index,:);
%                 ciplot(smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
%                 hold on
%                 ciplot(smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
%                 plot(mean(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff ,:)))),'k','LineWidth',2)
%    
%             else
%                 add = Res.Thetas(Parameters.adultsadd.Index,:);
%                 mult = Res.Thetas(Parameters.adultsmult.Index,:);
%                 ciplot(smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
%                 hold on
%                 ciplot(smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
%                 plot(mean(repmat(add',1,size(Paths,3))+repmat(mult',1,size(Paths,3)).*squeeze(exp(Paths(:,toplotdiff ,:)))),'k','LineWidth',2)
%             end
%             title('kids to adults')
%             ylabel('\beta_{12}(t)')
%             xlabel('t (weeks)')
            
            subplot(length(toplot)+2,2,length(toplot)*2+3:length(toplot)*2+4)
            
            try Parameters.holidays
                %for legend first            
                ciplot(1000*smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(1000*quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
                hold on
                ciplot(1000*smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),1000*smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
                plot(1000*mean(squeeze(exp(Paths(:,toplotdiff ,:)))),':k','LineWidth',1)
                plot(1000*(Data.Observations(toplot(i),:)),'.','Color',DotsColor,'MarkerSize',10)
                h = area(1000*[1+Data.Instants(7) 1+Data.Instants(13)],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
%                 legend('95% c.i.','50% c.i.','Post. mean','HPA data','Holidays')
            end
            
            nsmooth = 5;
            if strcmp(Parameters.Problem,'Marc3diff')
                toplotdiff = 13;
                ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
                hold on
                ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
                plot(mean(squeeze(exp(Paths(:,toplotdiff ,:)))),'k','LineWidth',2)
            elseif strcmp(Parameters.Problem,'Marc2diff')
                toplotdiff = 12;
                ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
                hold on
                try Parameters.holidays
                    ymax = max(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.975));
                    tmpymax = 3;%ymax*1.4;
                    h = area([1+Data.Instants(7) 1+Data.Instants(13)],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
                    h = area([1+Data.Instants(20) 1+Data.Instants(22)],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
                    set(gca,'Layer','top')
                    ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
                    ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
                    ylim([0 tmpymax])
                end
                ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
                plot(mean(squeeze(exp(Paths(:,toplotdiff ,:)))),':k','LineWidth',1)
                ylim([0 3])
            elseif strcmp(Parameters.Problem,'Marc2diffb')
                toplotdiff = 12;
                ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
                hold on
                try Parameters.holidays
                    ymax = max(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)+1))),0.975));
                    tmpymax = 3;%ymax*1.4;
                    h = area([1+Data.Instants(7) 1+Data.Instants(13)],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
                    h = area([1+Data.Instants(20) 1+Data.Instants(22)],[tmpymax tmpymax],'FaceColor',VeryLight,'EdgeColor',VeryLight);
                    set(gca,'Layer','top')
                    ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.025),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.975),nsmooth),[172,215,255]/255)
                    ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
                    ylim([0 tmpymax])
                end
                ciplot(smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.25),nsmooth),smooth(quantile(squeeze(exp(Paths(:,toplotdiff,:))),0.75),nsmooth),[100,153,251]/255)
                plot(mean(squeeze(exp(Paths(:,toplotdiff ,:)))),'k','LineWidth',2)
                ylim([0 3])
  
            else
                ciplot(repmat(quantile(squeeze(Res.Thetas(Parameters.beta22init.Index,:)),0.025),1,size(Paths,3)),repmat(quantile(squeeze(Res.Thetas(Parameters.beta22init.Index,:)),0.975),1,size(Paths,3)),[172,215,255]/255)
                hold on
                ciplot(repmat(quantile(squeeze(Res.Thetas(Parameters.beta22init.Index,:)),0.25),1,size(Paths,3)),repmat(quantile(squeeze(Res.Thetas(Parameters.beta22init.Index,:)),0.75),1,size(Paths,3)),[100,153,251]/255)
                plot(repmat(median(squeeze(Res.Thetas(Parameters.beta22init.Index,:))),1,size(Paths,3)),'k','LineWidth',2)
            end
            title('Effective contacts from adults to adults')
            ylabel('\beta_{adults \rightarrow adults}')
%             xlabel('t (weeks)')
             xlim([0 Data.Instants(end)])
            try Parameters.holidays
                xlim([0 (size(Data.Observations,2)-2)*7/Parameters.ComputationTStep])
                set(gca,'XTick',[delta*7/Parameters.ComputationTStep:delta*7/Parameters.ComputationTStep:length(Data.Dates)*21]-2*7/Parameters.ComputationTStep)
                %         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
                set(gca,'XTickLabel',dates)
            end
        end
        
        
        
        t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
        try
            plot(Res.BetatPath,'g','LineWidth',2) 
        end
        hold off
        xlim([0 Data.Instants(end)])
        try
            TicksInds = Data.Instants;
            set(gca,'XTick',TicksInds)
%             set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
        %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
%             set(gca,'XTickLabel',dates)
%             
%             set(gca,'XTick',[delta*7/Parameters.ComputationTStep:delta*7/Parameters.ComputationTStep:length(Data.Dates)*21]-2*7/Parameters.ComputationTStep)
%                 %         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
%                 set(gca,'XTickLabel',dates)
%             set(gca,'XTick',[delta*7/Parameters.ComputationTStep:delta*7/Parameters.ComputationTStep:length(Data.Dates)*21]-2*7/Parameters.ComputationTStep)
%             set(gca,'XTick',[108 159])
%             set(gca,'XTickLabel',{'13 Jul' , '1 Aug'})
%                 xlim([0 168])

            
        end
        
        ymax = 2;
    %     ylim([0.7 ymax])
        hold on
        yis1 = ymax/100:ymax/100:ymax*2/5;
        yis2 = ymax*3/5:ymax/100:ymax;
    %     plot(Data.Instants(8)*ones(size(yis1)),yis1,'r')
    %     plot(Data.Instants(8)*ones(size(yis2)),yis2,'r')
    %     text(Data.Instants(8),ymax*1/2,'sch. closure','HorizontalAlignment','center')
    %     plot(Data.Instants(14)*ones(size(yis1)),yis1,'r')
    %     plot(Data.Instants(14)*ones(size(yis2)),yis2,'r')
    %     text(Data.Instants(14),ymax*1/2,'end of holidays','HorizontalAlignment','center')
        hold off
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
