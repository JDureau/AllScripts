function mPost = SMCToOptimizewithPrior(Pars,Data,Model,Parameters)


Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = Pars(Parameters.(Names{i}).Index);
end

Parameters = UpdateParsTransfToNoTransf(Parameters);
Temp = EstimationSMCsmoothGen(Data, Model, Parameters);

% LogPrior = 0;
% for i = 1:length(Names)
%     tmp = Parameters.(Names{i}).Prior(Names{i},Parameters);
%     LogPrior = LogPrior +log(tmp);
% end
    
if Parameters.Correction
    mPost = -(Temp.LogLik+Temp.LogPrior-Temp.LogCorr);
elseif not(Parameters.Correction)
    mPost = -(Temp.LogLik+Temp.LogPrior);%-Temp.LogCorr);
end

disp(mPost)

% disp(Temp.LogLik)

% HIV
if strcmp(Parameters.Problem,'ImperialHIV')
    subplot(4,1,1)
    plot(mean(squeeze(Temp.CompletePaths(:,1,:))))
    hold on
    plot(quantile(squeeze(Temp.CompletePaths(:,1,:)),0.95),'r')
    plot(quantile(squeeze(Temp.CompletePaths(:,1,:)),0.05),'r')
    plot(Data.Instants,Data.Observations(7,:),'g')
    hold off
    set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    subplot(4,1,2)
    plot(mean(squeeze(Temp.CompletePaths(:,2,:))))
    hold on
    plot(quantile(squeeze(Temp.CompletePaths(:,2,:)),0.95),'r')
    plot(quantile(squeeze(Temp.CompletePaths(:,2,:)),0.05),'r')
    plot(Data.Instants,Data.Observations(8,:),'g')
    hold off
    set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    subplot(4,1,3)
    plot(Data.Instants,Temp.PosteriorMeansRecord(9,:))
    set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    title(Temp.LogLik)
    subplot(4,1,4)
    plot(mean(squeeze(Temp.CompletePaths(:,3,:))))
    hold on
    plot(quantile(squeeze(Temp.CompletePaths(:,3,:)),0.95),'r')
    plot(quantile(squeeze(Temp.CompletePaths(:,3,:)),0.05),'r')
%     try
        plot(Data.Fts,'g')
%     end
    hold off
    title('Ft')
    ylim([0 1])
    pause(0.01)
elseif strcmp(Parameters.Problem,'ImperialHIV2')
    subplot(4,1,1)
    plot(mean(squeeze(Temp.CompletePaths(:,1,:))))
    hold on
    plot(quantile(squeeze(Temp.CompletePaths(:,1,:)),0.95),'r')
    plot(quantile(squeeze(Temp.CompletePaths(:,1,:)),0.05),'r')
    plot(Data.Instants,Data.Observations(7,:),'g')
    hold off
    set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    subplot(4,1,2)
    plot(mean(squeeze(Temp.CompletePaths(:,2,:))))
    hold on
    plot(quantile(squeeze(Temp.CompletePaths(:,2,:)),0.95),'r')
    plot(quantile(squeeze(Temp.CompletePaths(:,2,:)),0.05),'r')
    plot(Data.Instants,Data.Observations(8,:),'g')
    hold off
    set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    subplot(4,1,3)
    if strcmp(Parameters.DiffusionType,'Logistic')
        mu = Parameters.CUinit.Value + Parameters.CUdelta.Value;
        plot(Data.Instants,mu*exp(Temp.PosteriorMeansRecord(9,:))/(1+exp(Temp.PosteriorMeansRecord(9,:))))
    else
        plot(Data.Instants,Temp.PosteriorMeansRecord(9,:))
    end
    
    set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    title(Temp.LogLik)
    subplot(4,1,4)
    if strcmp(Parameters.DiffusionType,'Bertallanfy')
        mu = Parameters.BRmu.Value ;
        m = Parameters.BRmm1.Value +1;
        tmp = mean(squeeze(Temp.CompletePaths(:,3,:)));
        plot(((1-m)*tmp+mu^(1-m)).^(1/(1-m)),'b')
        hold on
        tmp = quantile(squeeze(Temp.CompletePaths(:,3,:)),0.95);
        plot(((1-m)*tmp+mu^(1-m)).^(1/(1-m)),'r')
        tmp = quantile(squeeze(Temp.CompletePaths(:,3,:)),0.05);
        plot(((1-m)*tmp+mu^(1-m)).^(1/(1-m)),'r')
        try
            plot(Data.Fts,'g')
        end
    else
        tmp = mean(squeeze(Temp.CompletePaths(:,3,:)));
        plot(exp(tmp)./(1+exp(tmp)),'b')
        hold on
        tmp = quantile(squeeze(Temp.CompletePaths(:,3,:)),0.95);
        plot(exp(tmp)./(1+exp(tmp)),'r')
        tmp = quantile(squeeze(Temp.CompletePaths(:,3,:)),0.05);
        plot(exp(tmp)./(1+exp(tmp)),'r')
        try
            plot(Data.Fts,'g')
        end
    end
   

    hold off
%     plot(mean(squeeze(Temp.CompletePaths(:,3,:))))
%     hold off
%     hold on
%     plot(quantile(squeeze(Temp.CompletePaths(:,3,:)),0.95),'r')
%     plot(quantile(squeeze(Temp.CompletePaths(:,3,:)),0.05),'r')
% %     try
%         plot(Data.Fts,'g')
% %     end
%     hold off
    title('Ft')
    ylim([0 1])
    pause(0.01)
elseif strcmp(Parameters.Problem,'MarcFlu')   
    try
    % disp(Parameters.SigmaRW.Value)
        subplot(2,1,1)
        plot(Data.Instants,Temp.PosteriorMeansRecord(5,:))
        hold on
        plot(Data.Instants,Data.Observations(5,:),'g')
    %     plot(Data.Instants,Temp.Posterior975Record(3,:),'r')
    %     plot(Data.Instants,Temp.Posterior225Record(3,:),'r')
        hold off
        set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
        set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        subplot(2,1,2)
        plot(Data.Instants,exp(Temp.PosteriorMeansRecord(6,:)))
        hold on
    %     plot(Data.Instants,exp(Temp.Posterior975Record(6,:)),'r')
    %     plot(Data.Instants,exp(Temp.Posterior225Record(6,:)),'r')
        hold off
        set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
        set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        pause(0.01)
        subplot(2,1,2)
        if strcmp(Parameters.DiffusionType,'Add')
            plot(squeeze(mean(exp(Temp.CompletePaths(:,6,:)))))
            hold on
            temp = squeeze(Temp.CompletePaths(:,6,:));
            plot(squeeze(quantile(exp(temp),0.025)),'r')
            plot(squeeze(quantile(exp(temp),0.975)),'r')
            hold off
        else
            plot(squeeze(mean(exp(Temp.CompletePaths(:,6,:)))))
            hold on
            temp = squeeze(Temp.CompletePaths(:,6,:));
            plot(squeeze(quantile(exp(temp),0.025)),'r')
            plot(squeeze(quantile(exp(temp),0.975)),'r')
            hold off
        end
        set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
        set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        title(Temp.LogLik)
        pause(0.01)
    end
elseif strcmp(Parameters.Problem,'MarcFluPlusObs')   
    % disp(Parameters.SigmaRW.Value)
        subplot(3,1,1)
        plot(Data.Instants,Temp.PosteriorMeansRecord(5,:))
        hold on
        plot(Data.Instants,Data.Observations(5,:),'g')
    %     plot(Data.Instants,Temp.Posterior975Record(3,:),'r')
    %     plot(Data.Instants,Temp.Posterior225Record(3,:),'r')
        hold off
        set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
        set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        title(Temp.LogLik)
        subplot(3,1,2)
        plot(Data.Instants,exp(Temp.PosteriorMeansRecord(6,:)))
        hold on
    %     plot(Data.Instants,exp(Temp.Posterior975Record(6,:)),'r')
    %     plot(Data.Instants,exp(Temp.Posterior225Record(6,:)),'r')
        hold off
        set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
        set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        title('\beta_{t}')
        %     subplot(2,1,2)
    %     plot(squeeze(mean(exp(Temp.CompletePaths(:,4,:)))))
    %     hold on
    %     temp = squeeze(Temp.CompletePaths(:,4,:));
    %     plot(squeeze(quantile(exp(temp),0.025)),'r')
    %     plot(squeeze(quantile(exp(temp),0.975)),'r')
    %     hold off
    %     set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
    %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        
        subplot(3,1,3)
        if strcmp(Parameters.DiffusionType,'Add')
            plot(Data.Instants,exp(Temp.PosteriorMeansRecord(6,:)))
        else
            plot(Data.Instants,exp(Temp.PosteriorMeansRecord(7,:)))
        end
        hold on
    %     plot(Data.Instants,exp(Temp.Posterior975Record(6,:)),'r')
    %     plot(Data.Instants,exp(Temp.Posterior225Record(6,:)),'r')
        hold off
        set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
        set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])

        %     subplot(2,1,2)
    %     plot(squeeze(mean(exp(Temp.CompletePaths(:,4,:)))))
    %     hold on
    %     temp = squeeze(Temp.CompletePaths(:,4,:));
    %     plot(squeeze(quantile(exp(temp),0.025)),'r')
    %     plot(squeeze(quantile(exp(temp),0.975)),'r')
    %     hold off
    %     set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
    %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        title('\rho_{t}')
        pause(0.01)
        
elseif strcmp(Parameters.Problem,'Dengue')
    subplot(2,1,1)
    plot(Data.Instants,Temp.PosteriorMeansRecord(11,:))
    hold on 
    plot(Data.Instants,Data.Observations(11,:),'g')
    hold off
    subplot(2,1,2)
    plot(Data.Instants,exp(Temp.PosteriorMeansRecord(10,:)));
    pause(0.01)

elseif strcmp(Parameters.Problem,'Marc2diff')
        subplot(4,1,1)
        plot(Data.Instants,Data.Observations(9,:),'g')
        hold on
        plot(Data.Instants,Temp.PosteriorMeansRecord(9,:))
        hold off
        subplot(4,1,2)
        plot(Data.Instants,Data.Observations(10,:),'g')
        hold on
        plot(Data.Instants,Temp.PosteriorMeansRecord(10,:))
        hold off
        subplot(4,1,3)
        plot(exp(Temp.PosteriorMeansRecord(11,:)))
        subplot(4,1,4)
        plot(exp(Temp.PosteriorMeansRecord(12,:)))
        title(Temp.LogLik)
        pause(0.001)
 elseif strcmp(Parameters.Problem,'Marc3diff')
        subplot(5,1,1)
        plot(Data.Instants,Data.Observations(9,:),'g')
        hold on
        plot(Data.Instants,Temp.PosteriorMeansRecord(9,:))
        hold off
        subplot(5,1,2)
        plot(Data.Instants,Data.Observations(10,:),'g')
        hold on
        plot(Data.Instants,Temp.PosteriorMeansRecord(10,:))
        hold off
        subplot(5,1,3)
        plot(exp(Temp.PosteriorMeansRecord(11,:)))
        subplot(5,1,4)
        plot(exp(Temp.PosteriorMeansRecord(12,:)))
        subplot(5,1,5)
        plot(exp(Temp.PosteriorMeansRecord(13,:)))
        title(Temp.LogLik)
        pause(0.001)

end

