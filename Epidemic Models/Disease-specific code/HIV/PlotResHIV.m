function [] = PlotResHIV(Res,Parameters,FtPath)

Data = Res.Data;

try     
    Paths = Res.Paths;
catch
    Paths = Res.CompletePaths;
end        

if size(Paths,2)>4
    ToPlot = [7,8,9];
else
    ToPlot = [1,2,3];
end

if strcmp(Parameters.TypeWork,'Normal')
    if size(Paths,2)<40
%     figure(1)
clf
    subplot(3,1,1)
    ciplot(quantile(squeeze(Paths(:,ToPlot(1),:)),0.025),quantile(squeeze(Paths(:,ToPlot(1),:)),0.975),[172,215,255]/255)
    xlim([0 620])
    hold on
    ciplot(quantile(squeeze(Paths(:,ToPlot(1),:)),0.25),quantile(squeeze(Paths(:,ToPlot(1),:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(Paths(:,ToPlot(1),:))),'k','LineWidth',2)
    try
        for i = 2:size(Data.Observations,2)
            if Data.Observations(7,i)>0
                plot(Data.Instants(i)*ones(1,2),[Data.Observations(7,i)-2*sqrt(Data.Observations(7,i)*(100-Data.Observations(7,i))/400) Data.Observations(7,i)+2*sqrt(Data.Observations(7,i)*(100-Data.Observations(7,i))/400)],'r','LineWidth',3)
            end
        end
    end
     tmp = quantile(squeeze(Paths(:,ToPlot(1),:)),0.975);
     m = max(tmp);
    ylim([0 m*1.1])
    plot([456 456],[0 m*1.1],'--k','LineWidth',2)
    hold off
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    title('Prevalence among female sex workers','FontWeight','bold')
    xlabel('Time')
    ylabel('Prevalence')
%     legend('95% C.I','50% C.I.','Mean','Data')
    
%     figure(2)
    subplot(3,1,2)
    ciplot(quantile(squeeze(Paths(:,ToPlot(2),:)),0.025),quantile(squeeze(Paths(:,ToPlot(2),:)),0.975),[172,215,255]/255)
    xlim([0 620])
    hold on
    ciplot(quantile(squeeze(Paths(:,ToPlot(2),:)),0.25),quantile(squeeze(Paths(:,ToPlot(2),:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(Paths(:,ToPlot(2),:))),'k','LineWidth',2)
    try 
        for i = 2:size(Data.Observations,2)
            if Data.Observations(8,i)>0
                plot(Data.Instants(i)*ones(1,2),[Data.Observations(8,i)-2*sqrt(Data.Observations(8,i)*(100-Data.Observations(8,i))/400) Data.Observations(8,i)+2*sqrt(Data.Observations(8,i)*(100-Data.Observations(8,i))/400)],'r','LineWidth',3)
            end
        end
    end
     tmp = quantile(squeeze(Paths(:,ToPlot(2),:)),0.975);
     m = max(tmp);
    ylim([0 m*1.1])
    plot([456 456],[0 m*1.1],'--k','LineWidth',2)
    hold off
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])    
    xlabel('Time')
    ylabel('Prevalence')
    title('Prevalence among clients','FontWeight','bold')
%     legend('95% C.I','50% C.I.','Mean','Data') 
    
%     figure(3)
    subplot(3,1,3)
    if strcmp(Parameters.DiffusionType,'Logistic')
        paths = squeeze(Paths(:,ToPlot(3),:));
        indinit = Parameters.CUinit.Index;
        inddelta  = Parameters.CUdelta.Index;
        mu = Res.Thetas(indinit,:) + Res.Thetas(inddelta,:);
        Betas = min(1,max(0,diag(mu)*exp(paths)./(1+exp(paths))));
        ciplot(quantile(Betas,0.025),quantile(Betas,0.975),[172,215,255]/255)
        xlim([0 620])
        hold on
        ciplot(quantile(Betas,0.25),quantile(Betas,0.75),[100,153,251]/255)
        plot(mean(Betas),'k','LineWidth',2)

    else
        try
            if or(strcmp(Parameters.DiffusionType,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
                if Parameters.BRsigma.Estimated
                    tmp = squeeze(Paths(:,ToPlot(3),:));
%                     beta0s = squeeze(Res.Thetas(Parameters.BRbase.Index,:));
                    mus = squeeze(Res.Thetas(Parameters.BRmu.Index,:));
                    ms  = squeeze(Res.Thetas(Parameters.BRmm1.Index,:))+1;
%                     ts = squeeze(Res.Thetas(Parameters.BRtinfl.Index,:));
        %             Bs = 1-(beta0s./mus).^(1-ms);
        %             ks = 1./ts.*log(Bs./(1-ms));
                    FtBer = mean((repmat(1-ms',1,size(Paths,3)).*tmp+repmat(mus.^(1-ms)',1,size(Paths,3))).^(repmat(1./(1-ms)',1,size(Paths,3))));
                    Paths(:,ToPlot(3),:) = ((repmat(1-ms',1,size(Paths,3)).*tmp+repmat(mus.^(1-ms)',1,size(Paths,3))).^(repmat(1./(1-ms)',1,size(Paths,3))));

                end

    %             m = Parameters.BRm.Value;
    %             mu = Parameters.BRmu.Value;
    %             k = Parameters.k;
            elseif strcmp(Parameters.DiffusionType,'Sigmoid')
                if Parameters.Sigmsigma.Estimated
                    tmp = squeeze(Paths(:,3,:));
                    rate = squeeze(Res.Thetas(Parameters.Sigmrate.Index,:));
                    base = squeeze(Res.Thetas(Parameters.Sigmbase.Index,:));
                    mu = squeeze(Res.Thetas(Parameters.Sigmmu.Index,:));
                    tinfl = squeeze(Res.Thetas(Parameters.Sigmtinfl.Index,:));
                    c = 1./(1+exp(tinfl./rate));
                    b = (mu-base).*c./(1-c);
                    a = base - b;
                    indend = size(tmp,2);
                    a = repmat(a',1,indend);
                    b = repmat(b',1,indend);
                    c = repmat(c',1,indend);
    %                 FtSigmSto = mean(a + b./(c*(1+tmp)));
                    Paths(:,ToPlot(3),:) = (a + b./(c.*(1+tmp)));
                end
            

            elseif or(strcmp(Parameters.DiffusionType,'Add'),strcmp(Parameters.DiffusionType,'AddConstr'))
                Paths(:,ToPlot(3),:) = exp(Paths(:,ToPlot(3),:))./(1+exp(Paths(:,ToPlot(3),:)));
            elseif strcmp(Parameters.DiffusionType,'Sigmoid')
                base = Parameters.Sigmbase.Value;
                rate = Parameters.Sigmrate.Value;
                mu = Parameters.Sigmmu.Value;
                Paths = base + (mu-base)./(1+(Paths));
            end
            ciplot(quantile(squeeze(Paths(:,ToPlot(3),:)),0.025),quantile(squeeze(Paths(:,ToPlot(3),:)),0.975),[172,215,255]/255)
            xlim([0 620])
            hold on
            ciplot(quantile(squeeze(Paths(:,ToPlot(3),:)),0.25),quantile(squeeze(Paths(:,ToPlot(3),:)),0.75),[100,153,251]/255)
            plot(mean(squeeze(Paths(:,ToPlot(3),:))),'k','LineWidth',2)

        end
    end
%     for i = 1:9
%         plot(1:620,i*0.1*ones(1,620))
%     end
    if nargin == 3
        plot(FtPath,'--g','LineWidth',2)
    end
    plot([456 456],[0 1],'--k','LineWidth',2)
    ylim([0 1])
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    title('Estimated condom use frequency along time','FontWeight','bold')
    xlabel('Time')
    ylabel('Condom use frequency')
%     legend('95% C.I','50% C.I.','Mean')

    else
        if strcmp(Parameters.DiffusionType,'Logistic')
            paths = Paths;
            indinit = Parameters.CUinit.Index;
            inddelta  = Parameters.CUdelta.Index;
            mu = Res.Thetas(indinit,:) + Res.Thetas(inddelta,:);
            Betas = min(1,max(0,diag(mu)*exp(paths)./(1+exp(paths))));
            ciplot(quantile(Betas,0.025),quantile(Betas,0.975),[172,215,255]/255)
            xlim([0 620])
            hold on
            ciplot(quantile(Betas,0.25),quantile(Betas,0.75),[100,153,251]/255)
            plot(mean(Betas),'k','LineWidth',2)

        else
            try
                if or(strcmp(Parameters.DiffusionType,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
                    if Parameters.BRsigma.Estimated
                        tmp = Paths;
    %                     beta0s = squeeze(Res.Thetas(Parameters.BRbase.Index,:));
                        mus = squeeze(Res.Thetas(Parameters.BRmu.Index,:));
                        ms  = squeeze(Res.Thetas(Parameters.BRmm1.Index,:))+1;
    %                     ts = squeeze(Res.Thetas(Parameters.BRtinfl.Index,:));
            %             Bs = 1-(beta0s./mus).^(1-ms);
            %             ks = 1./ts.*log(Bs./(1-ms));
                        FtBer = mean((repmat(1-ms',1,size(Paths,2)).*tmp+repmat(mus.^(1-ms)',1,size(Paths,2))).^(repmat(1./(1-ms)',1,size(Paths,2))));
                        Paths = ((repmat(1-ms',1,size(Paths,2)).*tmp+repmat(mus.^(1-ms)',1,size(Paths,2))).^(repmat(1./(1-ms)',1,size(Paths,2))));

                    end

        %             m = Parameters.BRm.Value;
        %             mu = Parameters.BRmu.Value;
        %             k = Parameters.k;

                elseif or(strcmp(Parameters.DiffusionType,'Add'),strcmp(Parameters.DiffusionType,'AddConstr'))
                    Paths = exp(Paths)./(1+exp(Paths));
                
                end
                clf
                ciplot(quantile(squeeze(Paths),0.025),quantile(squeeze(Paths),0.975),[172,215,255]/255)
                xlim([0 620])
                hold on
                ciplot(quantile(squeeze(Paths),0.25),quantile(squeeze(Paths),0.75),[100,153,251]/255)
                plot(mean(squeeze(Paths)),'k','LineWidth',2)

            end
        end
    %     for i = 1:9
    %         plot(1:620,i*0.1*ones(1,620))
    %     end
        if nargin == 3
            plot(FtPath,'--g','LineWidth',2)
        end
        plot([456 456],[0 1],'--k','LineWidth',2)
        ylim([0 1])
        set(gca,'XTick',[0:(600)/5:(600)])
        set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
        title('Estimated condom use frequency along time','FontWeight','bold')
        xlabel('Time')
        ylabel('Condom use frequency')
    end
        

elseif strcmp(Parameters.TypeWork,'ISSTDR Poster')
%     figure(1)
    clf
    subplot(3,1,1)
    ciplot(quantile(squeeze(Paths(:,ToPlot(1),:)),0.025),quantile(squeeze(Paths(:,ToPlot(1),:)),0.975),[172,215,255]/255)
    xlim([0 620])
    hold on
    ciplot(quantile(squeeze(Paths(:,ToPlot(1),:)),0.25),quantile(squeeze(Paths(:,ToPlot(1),:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(Paths(:,ToPlot(1),:))),'k','LineWidth',2)
    for i = 2:size(Data.Observations,2)
        if Data.Observations(7,i)>0
            plot(Data.Instants(i)*ones(1,2),[Parameters.ObsMin(i-1) Parameters.ObsMax(i-1)]*100,'r','LineWidth',3)
        end
    end
    
    tmp = quantile(squeeze(Paths(:,ToPlot(1),:)),0.975);
    m = max(tmp);
    ylim([0,m*1.1])
    plot([456 456],[0 m*1.1],'--k','LineWidth',2)
    hold off
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'],'FontSize',14)
    title('Prevalence among Female Sex Workers','FontWeight','bold','FontSize',16)
%     gtext({'\leftarrow Start of','     Avahan'},'FontSize',14)

%     xlabel('Time')
    ylabel('Prevalence','FontSize',15)
%     legend('95% C.I','50% C.I.','Mean','Data')
    

    subplot(3,1,2)
    ciplot(quantile(squeeze(Paths(:,ToPlot(2),:)),0.025),quantile(squeeze(Paths(:,ToPlot(2),:)),0.975),[172,215,255]/255)
    xlim([0 620])
    hold on
    ciplot(quantile(squeeze(Paths(:,ToPlot(2),:)),0.25),quantile(squeeze(Paths(:,ToPlot(2),:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(Paths(:,ToPlot(2),:))),'k','LineWidth',2)
    for i = 2:size(Data.Observations,2)
        if Data.Observations(8,i)>0
            plot(Data.Instants(i)*ones(1,2),[Parameters.ObsMin(i-1) Parameters.ObsMax(i-1)]*100,'r','LineWidth',3)
        end
    end
    
    tmp = quantile(squeeze(Paths(:,ToPlot(2),:)),0.975);
    m = max(tmp);
    ylim([0,m*1.1])
    plot([456 456],[0 m*1.1],'--k','LineWidth',2)
    hold off
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'],'FontSize',14)
    title('Prevalence among Clients','FontWeight','bold','FontSize',16)
%     gtext({'\leftarrow Start of','     Avahan'},'FontSize',14)

%     xlabel('Time')
    ylabel('Prevalence','FontSize',15)

%   
    
%     figure(3)
    subplot(3,1,3)
    tmp = squeeze(Paths(:,ToPlot(3),:));
    ciplot(quantile(exp(tmp)./(1+exp(tmp)),0.025),quantile(exp(tmp)./(1+exp(tmp)),0.975),[172,215,255]/255)
    xlim([0 620])
    hold on
    ciplot(quantile(exp(tmp)./(1+exp(tmp)),0.25),quantile(exp(tmp)./(1+exp(tmp)),0.75),[100,153,251]/255)
    plot(mean(exp(tmp)./(1+exp(tmp))),'k','LineWidth',2)
%     for i = 1:9
%         plot(1:620,i*0.1*ones(1,620))
%     end
    if nargin == 3
        plot(FtPath,'--g','LineWidth',2)
    end
    tmp = quantile(squeeze(Paths(:,ToPlot(3),:)),0.975);
    m = max(tmp);
    ylim([0,m*1.1])
    plot([456 456],[0 m*1.1],'--k','LineWidth',2)
    hold off
    ylim([0 1])
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'],'FontSize',14)
    title('Estimated Condom Use','FontWeight','bold','FontSize',16)
%     xlabel('Time')
    ylabel('Condom Use','FontSize',15)
%     legend('95% C.I','50% C.I.','Mean')
    
elseif strcmp(Parameters.TypeWork,'ISSTDR Poster Val')
   
%   
    clf
%     figure(3)
    ciplot(quantile(squeeze(Paths(:,ToPlot(3),:)),0.025),quantile(squeeze(Paths(:,ToPlot(3),:)),0.975),[172,215,255]/255)
    xlim([0 620])
    hold on
    ciplot(quantile(squeeze(Paths(:,ToPlot(3),:)),0.25),quantile(squeeze(Paths(:,ToPlot(3),:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(Paths(:,ToPlot(3),:))),'k','LineWidth',2)
    m = size(Paths,3);
    xis = 0:m/(308-1):m;
    plot(xis,Res.SimTraj(1:308),'g','LineWidth',4)
%     for i = 1:9
%         plot(1:620,i*0.1*ones(1,620))
%     end
    if nargin == 3
        plot(FtPath,'--g','LineWidth',2)
    end
    hold off
    ylim([0 1])
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'],'FontSize',16)
    title('Estimated Condom Use','FontWeight','bold','FontSize',18)
%     xlabel('Time')
    ylabel('Condom Use','FontSize',16)
%     legend('95% C.I','50% C.I.','Mean')


  
elseif strcmp(Parameters.TypeWork,'Boston Examples')
   
    SmLgth = 15;
%     clf
%     figure(3)


%  clf
if Parameters.PlotIndex == 5
    subplot(4,2,1)
    ciplot(quantile(squeeze(Paths(:,ToPlot(1),:)),0.025),quantile(squeeze(Paths(:,ToPlot(1),:)),0.975),[172,215,255]/255)
    xlim([0 620])
    hold on
    ciplot(quantile(squeeze(Paths(:,ToPlot(1),:)),0.25),quantile(squeeze(Paths(:,ToPlot(1),:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(Paths(:,ToPlot(1),:))),'k','LineWidth',2)
    for i = 2:size(Data.Observations,2)
        if Data.Observations(7,i)>0
            plot(Data.Instants(i)*ones(1,2),[Parameters.ObsMin(i-1) Parameters.ObsMax(i-1)]*100,'r','LineWidth',3)
        end
    end
    
    tmp = quantile(squeeze(Paths(:,ToPlot(1),:)),0.975);
    m = max(tmp);
    ylim([0,m*1.1])
    plot([456 456],[0 m*1.1],'--k','LineWidth',2)
    hold off
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'],'FontSize',14)
    title('Prevalence among Female Sex Workers','FontWeight','bold','FontSize',16)
%     gtext({'\leftarrow Start of','     Avahan'},'FontSize',14)

%     xlabel('Time')
    ylabel('Prevalence','FontSize',15)
%     legend('95% C.I','50% C.I.','Mean','Data')
    

    subplot(4,2,2)
    ciplot(quantile(squeeze(Paths(:,ToPlot(2),:)),0.025),quantile(squeeze(Paths(:,ToPlot(2),:)),0.975),[172,215,255]/255)
    xlim([0 620])
    hold on
    ciplot(quantile(squeeze(Paths(:,ToPlot(2),:)),0.25),quantile(squeeze(Paths(:,ToPlot(2),:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(Paths(:,ToPlot(2),:))),'k','LineWidth',2)
    for i = 2:size(Data.Observations,2)
        if Data.Observations(8,i)>0
            plot(Data.Instants(i)*ones(1,2),[Parameters.ObsMin(i-1) Parameters.ObsMax(i-1)]*100,'r','LineWidth',3)
        end
    end
    
    tmp = quantile(squeeze(Paths(:,ToPlot(2),:)),0.975);
    m = max(tmp);
    ylim([0,m*1.1])
    plot([456 456],[0 m*1.1],'--k','LineWidth',2)
    hold off
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'],'FontSize',14)
    title('Prevalence among Clients','FontWeight','bold','FontSize',16)
%     gtext({'\leftarrow Start of','     Avahan'},'FontSize',14)

%     xlabel('Time')
    ylabel('Prevalence','FontSize',15)
end





    try
        subplot(4,2,2+Parameters.PlotIndex)
    catch
        clf
    end
    
    
    n = min(5000,size(Res.Thetas,2));
    Res.Thetas = Res.Thetas(:,1:n);
    Paths = Paths(1:n,:,:);

    try
        if or(strcmp(Parameters.DiffusionType,'Bertallanfy'),and(strcmp(Parameters.DiffusionType,'Sigmoid'),not(Parameters.Sigmsigma.Estimated)))
            if Parameters.BRsigma.Estimated
                tmp = squeeze(Paths(:,ToPlot(3),:));
                beta0s = squeeze(Res.Thetas(Parameters.BRbase.Index,:));
                mus = squeeze(Res.Thetas(Parameters.BRmu.Index,:));
                ms = squeeze(Res.Thetas(Parameters.BRmm1.Index,:))+1;
                ts = squeeze(Res.Thetas(Parameters.BRtinfl.Index,:));
    %             Bs = 1-(beta0s./mus).^(1-ms);
    %             ks = 1./ts.*log(Bs./(1-ms));
                FtBer = mean((repmat(1-ms',1,size(Paths,3)).*tmp+repmat(mus.^(1-ms)',1,size(Paths,3))).^(repmat(1./(1-ms)',1,size(Paths,3))));
                Paths(:,ToPlot(3),:) = ((repmat(1-ms',1,size(Paths,3)).*tmp+repmat(mus.^(1-ms)',1,size(Paths,3))).^(repmat(1./(1-ms)',1,size(Paths,3))));
                
            else
                Paths = (squeeze(Paths));
                
            end

%             m = Parameters.BRm.Value;
%             mu = Parameters.BRmu.Value;
%             k = Parameters.k;
        elseif strcmp(Parameters.DiffusionType,'Sigmoid')
            if Parameters.Sigmsigma.Estimated
                tmp = squeeze(Paths(:,3,:));
                rate = squeeze(Res.Thetas(Parameters.Sigmrate.Index,:));
                base = squeeze(Res.Thetas(Parameters.Sigmbase.Index,:));
                mu = squeeze(Res.Thetas(Parameters.Sigmmu.Index,:));
                tinfl = squeeze(Res.Thetas(Parameters.Sigmtinfl.Index,:));
                c = 1./(1+exp(tinfl./rate));
                b = (mu-base).*c./(1-c);
                a = base - b;
                a = repmat(a',1,size(Paths,3));
                b = repmat(b',1,size(Paths,3));
                c = repmat(c',1,size(Paths,3));
%                 FtSigmSto = mean(a + b./(c*(1+tmp)));
                Paths(:,ToPlot(3),:) = (a + b./(c.*(1+tmp)));
            else
                Paths = (squeeze(Paths));
            end

        elseif or(strcmp(Parameters.DiffusionType,'Add'),strcmp(Parameters.DiffusionType,'AddConstr'))
            Paths(:,ToPlot(3),:) = exp(Paths(:,ToPlot(3),:))./(1+exp(Paths(:,ToPlot(3),:)));
        end
    end

    if size(Paths,2)<20
        ciplot(smooth(quantile(squeeze(Paths(:,ToPlot(3),:)),0.025),SmLgth,'median'),smooth(quantile(squeeze(Paths(:,ToPlot(3),:)),0.975),SmLgth,'median'),[172,215,255]/255)
        xlim([0 620])
        hold on
        ciplot(smooth(quantile(squeeze(Paths(:,ToPlot(3),:)),0.25),SmLgth,'median'),smooth(quantile(squeeze(Paths(:,ToPlot(3),:)),0.75),SmLgth,'median'),[100,153,251]/255)
        plot(smooth(median(squeeze(Paths(:,ToPlot(3),:))),SmLgth,'median'),'k','LineWidth',2)
        m = size(Paths,3);
        xis = 0:m/(308-1):m;
        try
            plot(Res.Data.BuiltTraj(:,9),'g','LineWidth',4)
        end
        tmp = quantile(squeeze(Paths(:,ToPlot(3),:)),0.975);    

    else
        ciplot(smooth(quantile(squeeze(Paths),0.025),SmLgth,'median'),smooth(quantile(squeeze(Paths),0.975),SmLgth,'median'),[172,215,255]/255)
        xlim([0 620])
        hold on
        ciplot(smooth(quantile(squeeze(Paths),0.25),SmLgth,'median'),smooth(quantile(squeeze(Paths),0.75),SmLgth,'median'),[100,153,251]/255)
        plot(smooth(median(squeeze(Paths)),SmLgth,'median'),'k','LineWidth',2)
        tmp = quantile(squeeze(Paths),0.975);    

    end
    
    m = max(tmp);
    ylim([0 1])
    plot([456 456],[0 1],'--k','LineWidth',2)
%     try
%         plot(Parameters.Sigm,'y','LineWidth',4)
%     end
%     for i = 1:9
%         plot(1:620,i*0.1*ones(1,620))
%     end
    if nargin == 3
        plot(FtPath,'--g','LineWidth',2)
    end
    hold off
    ylim([0 1])
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'],'FontSize',16)
    if Parameters.PlotIndex == 1
        title('Estimated Condom Use','FontWeight','bold','FontSize',18)
    end
%     xlabel('Time')
%     ylabel('Condom Use','FontSize',16)
%     legend('95% C.I','50% C.I.','Mean')
    
elseif strcmp(Parameters.TypeWork,'Boston Examples HIV2')
   
    SmLgth = 15;
    clf
%     figure(3)
    

    subplot(2,1,1)
    ciplot(quantile(squeeze(Paths(:,ToPlot(1),:)),0.025),quantile(squeeze(Paths(:,ToPlot(1),:)),0.975),[172,215,255]/255)
    xlim([0 620])
    hold on
    ciplot(quantile(squeeze(Paths(:,ToPlot(1),:)),0.25),quantile(squeeze(Paths(:,ToPlot(1),:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(Paths(:,ToPlot(1),:))),'k','LineWidth',2)
    for i = 2:size(Data.Observations,2)
        if Data.Observations(7,i)>0
            plot(Data.Instants(i)*ones(1,2),[Data.Observations(7,i)*(0.95) Data.Observations(7,i)*(1.05)],'r','LineWidth',3)
        end
    end
    
    tmp = quantile(squeeze(Paths(:,ToPlot(1),:)),0.975);
    m = max(tmp);
    ylim([0,m*1.1])
    plot([456 456],[0 m*1.1],'--k','LineWidth',2)
    hold off
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'],'FontSize',14)
    title('Prevalence among Female Sex Workers','FontWeight','bold','FontSize',16)
%     gtext({'\leftarrow Start of','     Avahan'},'FontSize',14)

%     xlabel('Time')
    ylabel('Prevalence','FontSize',15)

    subplot(2,1,2)
%     subplot(1,1,1)

    tmp = squeeze(Paths(:,ToPlot(3),:));
    ciplot(smooth(quantile(exp(tmp)./(1+exp(tmp)),0.025),SmLgth,'median'),smooth(quantile(exp(tmp)./(1+exp(tmp)),0.975),SmLgth,'median'),[172,215,255]/255)
    xlim([0 620])
    hold on
    ciplot(smooth(quantile(exp(tmp)./(1+exp(tmp)),0.25),SmLgth,'median'),smooth(quantile(exp(tmp)./(1+exp(tmp)),0.75),SmLgth,'median'),[100,153,251]/255)
    plot(smooth(mean(exp(tmp)./(1+exp(tmp))),SmLgth,'median'),'k','LineWidth',2)
    m = size(Paths,3);
    xis = 0:m/(308-1):m;
    try
%         plot(xis,0.65*ones(size(xis)),'g','LineWidth',4)
        plot(Res.Data.BuiltTraj(:,9),'g','LineWidth',4)
        plot(Parameters.Sigm,'r','LineWidth',4)
    end
%     try
%         plot(Parameters.Sigm,'y','LineWidth',4)
%     end
%     for i = 1:9
%         plot(1:620,i*0.1*ones(1,620))
%     end
    if nargin == 3
        plot(FtPath,'--g','LineWidth',2)
    end
    hold off
    ylim([0 1])
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'],'FontSize',16)
    title('Estimated Condom Use','FontWeight','bold','FontSize',38)
%     xlabel('Time')
    ylabel('Condom Use','FontSize',36)
    set(gca,'Fontsize',36)

%     legend('95% C.I','50% C.I.','Mean')


elseif strcmp(Parameters.TypeWork,'PlayObs')
    figure(1)
    subplot(3,1,1)
    ciplot(quantile(squeeze(Paths(:,ToPlot(1),:)),0.025),quantile(squeeze(Paths(:,ToPlot(1),:)),0.975),[172,215,255]/255)
%     xlim([0 620])
    hold on
    ciplot(quantile(squeeze(Paths(:,ToPlot(1),:)),0.25),quantile(squeeze(Paths(:,ToPlot(1),:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(Paths(:,ToPlot(1),:))),'k','LineWidth',2)
    for i = 2:length(Data.Instants)
        if Data.ObservedVariables(i) == 7
            plot(Data.Instants(i)*ones(1,3),Data.Observations(7,i)*[0.9 1 1.1],'col',[55,233,30]/255,'LineWidth',3)
        end
    end
    hold off
%     set(gca,'XTick',[0:(600)/5:(600)])
%     set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    title('Female Sex Workers')

    subplot(3,1,2)
    ciplot(quantile(squeeze(Paths(:,ToPlot(2),:)),0.025),quantile(squeeze(Paths(:,ToPlot(2),:)),0.975),[172,215,255]/255)
%     xlim([0 620])
    hold on
    ciplot(quantile(squeeze(Paths(:,ToPlot(2),:)),0.25),quantile(squeeze(Paths(:,ToPlot(2),:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(Paths(:,ToPlot(2),:))),'k','LineWidth',2)
    for i = 2:length(Data.Instants)
        if Data.ObservedVariables(i) == 8
            plot(Data.Instants(i)*ones(1,3),Data.Observations(8,i)*[0.9 1 1.1],'col',[55,233,30]/255,'LineWidth',3)
        end
    end
    hold off
%     set(gca,'XTick',[0:(600)/5:(600)])
%     set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])    

    title('Clients')

    subplot(3,1,3)
    ciplot(quantile(squeeze(Paths(:,ToPlot(3),:)),0.025),quantile(squeeze(Paths(:,ToPlot(3),:)),0.975),[172,215,255]/255)
    xlim([0 620])
    hold on
    ciplot(quantile(squeeze(Paths(:,ToPlot(3),:)),0.25),quantile(squeeze(Paths(:,ToPlot(3),:)),0.75),[100,153,251]/255)
    plot(mean(squeeze(Paths(:,ToPlot(3),:))),'k','LineWidth',2)
    set(gca,'XTick',[0:(600)/5:(600)])
    set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'],'FontSize',16)
    title('Estimated Condom Use','FontWeight','bold','FontSize',18)

%     for i = 1:9
%         plot(1:620,i*0.1*ones(1,620))
%     end
    if nargin == 3
        plot(FtPath,'--g','LineWidth',4)
    end
    hold off
    ylim([0 1])
%     set(gca,'XTick',[0:(600)/5:(600)])
%     set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
    title('Condom use frequency')
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
