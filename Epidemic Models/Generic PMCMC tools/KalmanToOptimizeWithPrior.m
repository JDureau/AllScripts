function mPost= KalmanToOptimizeWithPrior(Pars,Data,Model,Parameters)

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = Pars(Parameters.(Names{i}).Index);
end

Parameters = UpdateParsTransfToNoTransf(Parameters);
CheckParametersGen(Parameters)



try
    Temp = EstimationEKFGen(Data, Model, Parameters);
    logs = [];
    LogPrior = 0;
    for i = 1:length(Names)
        temp = Parameters.(Names{i}).Prior(Names{i},Parameters);
       
        logs(i) = log(temp);
        LogPrior = LogPrior +log(temp);
%         Names{i}
%         LogPrior
    end
    
    mPost = -(Temp.LogLik+LogPrior-Parameters.Correction*Temp.LogCorr);
    if Parameters.Correction
        mPost = -(Temp.LogLik+LogPrior-Temp.LogCorr);
    elseif not(Parameters.Correction)
        mPost = -(Temp.LogLik+LogPrior);%-Temp.LogCorr);
    end
%     mPost
  
catch
    mPost= -Inf;
end

try
clf
if strcmp(Parameters.Problem,'ImperialHIV2')
    rd = rand(1);
    if rd >0
        subplot(3,1,1)
        plot(Data.Instants,Temp.PosteriorMeans(7,:))
        hold on
        plot(Data.Instants,Data.Observations(7,:),'g')
        plot(Data.Instants,Temp.Posterior975(7,:),'r')
        plot(Data.Instants,Temp.Posterior025(7,:),'r')
        hold off
        set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
        set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
        subplot(3,1,2)
        plot(Data.Instants,Temp.PosteriorMeans(8,:))
        hold on
        plot(Data.Instants,Data.Observations(8,:),'g')
        plot(Data.Instants,Temp.Posterior975(8,:),'r')
        plot(Data.Instants,Temp.Posterior025(8,:),'r')
        hold off
        set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
        set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
        subplot(3,1,3)
        if or(strcmp(Parameters.DiffusionType,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
            mu = Parameters.BRmu.Value ;
            m = Parameters.BRmm1.Value +1;
            tmp = (min(0,Temp.PosteriorMeans(9,:)));
            plot(Data.Instants,((1-m)*tmp+mu^(1-m)).^(1/(1-m)),'g')
            hold on
            tmp = min(0,Temp.Posterior975(9,:));
            plot(Data.Instants,((1-m)*tmp+mu^(1-m)).^(1/(1-m)),'r')
            tmp = min(0,Temp.Posterior025(9,:));
            plot(Data.Instants,((1-m)*tmp+mu^(1-m)).^(1/(1-m)),'r')
            hold off
        elseif strcmp(Parameters.DiffusionType,'Sigmoid')
            base = Parameters.Sigmbase.Value;
            mu   = Parameters.Sigmmu.Value;
            rate = Parameters.Sigmrate.Value;
            tinfl = Parameters.Sigmbase.Value;
            tmp = max(0,Temp.PosteriorMeans(9,:));
            plot(Data.Instants,base+(mu-base)./(1+tmp),'g')
            hold on
            tmp =  max(0,Temp.Posterior975(9,:));
            plot(Data.Instants,base+(mu-base)./(1+tmp),'r')
            tmp =  max(0,Temp.Posterior025(9,:));
            plot(Data.Instants,base+(mu-base)./(1+tmp),'r')
            hold off
        else
            tmp = Temp.PosteriorMeans(9,:);
            plot(Data.Instants,exp(tmp)./(1+exp(tmp)),'g')
            hold on
            tmp = Temp.Posterior975(9,:);
            plot(Data.Instants,exp(tmp)./(1+exp(tmp)),'r')
            tmp = Temp.Posterior025(9,:);
            plot(Data.Instants,exp(tmp)./(1+exp(tmp)),'r')
            hold off
        end
        ylim([0 1])
        title(['LogPost: ' num2str(-mPost) ]);%' LogPrior: ' num2str(LogPrior)])
        pause(0.001)
    end
elseif strcmp(Parameters.Problem,'GoogleFlu')
% for SEIR
    try
        subplot(2,1,1)
        plot(Data.Instants,(Temp.PosteriorMeans(5,:)))
        hold on
        plot(Data.Instants,Data.Observations(5,:),'g')
        hold off
        set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
        set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        subplot(2,1,2)
        plot(Data.Instants,exp(Temp.PosteriorMeans(6,:)))
        set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
        set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        title(Temp.LogLik)
        pause(0.001)
    end
    elseif strcmp(Parameters.Problem,'GoogleFlu')
% for SEIR
    try
        subplot(2,1,1)
        plot(Data.Instants,(Temp.PosteriorMeans(5,:)))
        hold on
        plot(Data.Instants,Data.Observations(5,:),'g')
        hold off
        set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
        set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        subplot(2,1,2)
        plot(Data.Instants,exp(Temp.PosteriorMeans(6,:)))
        set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
        set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        title(Temp.LogLik)
        pause(0.001)
    end
elseif strcmp(Parameters.Problem,'Dengue')
% for SEIR
    try
        rep1 = Parameters.rep1.Value;
        rep2 = Parameters.rep2.Value;
        subplot(2,1,1)
        plot(Data.Instants,max(eps,rep1*rep2*Temp.PosteriorMeans(11,:)))
        hold on
        plot(Data.Instants,Data.Observations(11,:),'g')
%         plot(Data.Instants,rep1*rep2*Temp.Posterior975(11,:),'r')
%         plot(Data.Instants,rep1*rep2*Temp.Posterior025(11,:),'r')
        hold off
        sizeStep = Data.NbComputingSteps(2);
        subplot(2,1,2)
        plot(Data.Instants,exp(Temp.PosteriorMeans(10,:)))
        title(['LogLik: ' num2str(Temp.LogLik)  ' - LogPrior: ' num2str(Temp.LogPrior)])
%         subplot(2,1,1)
%         plot(log(Temp.Liks))
        pause(0.001)
    end
elseif strcmp(Parameters.Problem,'SIRsimple')
% for SEIR
    try
        
        subplot(2,1,1)
        plot(Data.Observations(3,:),'g')
        hold on
        plot(0.24*Temp.PosteriorMeans(3,:))
        hold off
        subplot(2,1,2)
        plot(exp(Temp.PosteriorMeans(4,:)))
        title(Temp.LogLik)
        pause(0.001)
    end
elseif strcmp(Parameters.Problem,'MarcFluPlusObs')
    subplot(2,1,1)
        plot(Data.Observations(5,:),'g')
        hold on
        plot(Temp.PosteriorMeans(5,:))
        hold off
        subplot(2,1,2)
        plot(exp(Temp.PosteriorMeans(6,:)))
        title(Temp.LogLik)
        pause(0.001)
elseif strcmp(Parameters.Problem,'Marc1diff')
        subplot(2,1,1)
        plot(Data.Observations(5,:),'g')
        hold on
        plot(Temp.PosteriorMeans(5,:))
        hold off
        subplot(2,1,2)
        plot(exp(Temp.PosteriorMeans(6,:)))
        title(Temp.LogLik)
        pause(0.001)
elseif strcmp(Parameters.Problem,'Marc2diff')
        subplot(4,1,1)
        plot(Data.Observations(9,:),'g')
        hold on
        plot(Temp.PosteriorMeans(9,:))
        hold off
        subplot(4,1,2)
        plot(Data.Observations(10,:),'g')
        hold on
        plot(Temp.PosteriorMeans(10,:))
        hold off
        subplot(4,1,3)
        plot(exp(Temp.PosteriorMeans(11,:)))
        subplot(4,1,4)
        plot(exp(Temp.PosteriorMeans(12,:)))
        title(Temp.LogLik)
        pause(0.001)
        
 elseif strcmp(Parameters.Problem,'Marc3diff')
        subplot(5,1,1)
        plot(Data.Observations(9,:),'g')
        hold on
        plot(Temp.PosteriorMeans(9,:))
        hold off
        subplot(5,1,2)
        plot(Data.Observations(10,:),'g')
        hold on
        plot(Temp.PosteriorMeans(10,:))
        hold off
        subplot(5,1,3)
        plot(exp(Temp.PosteriorMeans(11,:)))
        subplot(5,1,4)
        plot(exp(Temp.PosteriorMeans(12,:)))
        subplot(5,1,5)
        plot(exp(Temp.PosteriorMeans(13,:)))
        title(Temp.LogLik)
        pause(0.001)
end
end