function mPost = SMCToOptimizewithPriorNoTransf(Pars,Data,Model,Parameters)


Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).Value = Pars(Parameters.(Names{i}).Index);
    disp(Pars(Parameters.(Names{i}).Index))
end

Parameters = UpdateParsNoTransfToTransf(Parameters);
Temp = EstimationSMCsmoothGen(Data, Model, Parameters);

LogPrior = 0;
ParametersForPriors = Parameters;
for i = 1:length(Names)
    LogPrior = LogPrior +log(eval(Parameters.(Names{i}).Prior));
end
    
mPost = -(Temp.LogLik+LogPrior);

% HIV
try
subplot(4,1,1)
plot(Data.Instants,Temp.PosteriorMeansRecord(7,:))
hold on
plot(Data.Instants,Data.Observations(7,:),'g')
hold off
set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
subplot(4,1,2)
plot(Data.Instants,Temp.PosteriorMeansRecord(8,:))
hold on
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
plot(quantile(squeeze(Temp.CompletePaths(:,3,:)),0.25),'r')
plot(quantile(squeeze(Temp.CompletePaths(:,3,:)),0.75),'r')
hold off
title('Ft')
ylim([0 1])
pause(0.01)
catch

    
% disp(Parameters.SigmaRW.Value)
    subplot(2,1,1)
    plot(Data.Instants,Temp.PosteriorMeansRecord(3,:))
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

    %     subplot(2,1,2)
%     plot(squeeze(mean(exp(Temp.CompletePaths(:,4,:)))))
%     hold on
%     temp = squeeze(Temp.CompletePaths(:,4,:));
%     plot(squeeze(quantile(exp(temp),0.025)),'r')
%     plot(squeeze(quantile(exp(temp),0.975)),'r')
%     hold off
%     set(gca,'XTick',[0:(Data.Instants(end))/12:(Data.Instants(end))])
%     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
    title(Temp.LogLik)
    pause(0.01)

end

