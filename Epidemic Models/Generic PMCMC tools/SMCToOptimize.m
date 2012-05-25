function mLogLik = SMCToOptimize(Pars,Data,Model,Parameters)


Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = Pars(Parameters.(Names{i}).Index);
end

Parameters = UpdateParsTransfToNoTransf(Parameters);
Temp = EstimationSMCsmoothGen(Data, Model, Parameters);

mLogLik = -Temp.LogLik;
subplot(3,1,1)
plot(Data.Instants,Temp.PosteriorMeansRecord(6,:))
hold on
plot(Data.Instants,Data.Observations(6,:),'g')
hold off
set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
subplot(3,1,2)
plot(Data.Instants,Temp.PosteriorMeansRecord(7,:))
hold on
plot(Data.Instants,Data.Observations(7,:),'g')
hold off
set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
subplot(3,1,3)
plot(Data.Instants,Temp.PosteriorMeansRecord(8,:))
set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
title(Temp.LogLik)
pause(0.01)