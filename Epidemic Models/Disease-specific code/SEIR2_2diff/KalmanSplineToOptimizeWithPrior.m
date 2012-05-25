function mPost= KalmanSplineToOptimizeWithPrior(Pars,Data,Model,Parameters)

n = Parameters.NbSplines;
ts = Pars(1:n);
yis = Pars(n+1:2*n);
[ts,I,J] = UNIQUE(ts);


yis = yis(I);
Betas = spline(ts,yis,1:sum(Data.NbComputingSteps));

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = Pars(2*n+Parameters.(Names{i}).Index);
end

Parameters = UpdateParsTransfToNoTransf(Parameters);
CheckParametersGen(Parameters)
Parameters.Betas = Betas;

% disp(Parameters.betainit)
Temp = EstimationEKFGen(Data, Model, Parameters);
try
    Temp = EstimationEKFGen(Data, Model, Parameters);
%     disp(Temp.LogLik + Temp.LogPrior)
%     logs = [];
%     for i = 1:length(Names)
%         temp = Parameters.(Names{i}).Prior(Names{i},Parameters);
%         
%         logs(i) = log(temp);
%         LogPrior = LogPrior +log(temp);
%     end
    
    mPost = -(Temp.LogLik+ Temp.LogPrior-Temp.LogCorr);
%     disp(logs)
%     disp(['Prior ' num2str(LogPrior)])
catch
    mPost = Inf;
    disp('EstimationEKF crashed')
end

% % for HIV
clf
if strcmp(Parameters.Problem,'ImperialHIV')
    rd = rand(1);
    if rd >0
        subplot(3,1,1)
        plot(Data.Instants,Temp.PosteriorMeans(7,:))
        hold on
        plot(Data.Instants,Data.Observations(7,:),'g')
        hold off
        set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
        set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
        subplot(3,1,2)
        plot(Data.Instants,Temp.PosteriorMeans(8,:))
        hold on
        plot(Data.Instants,Data.Observations(8,:),'g')
        hold off
        set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
        set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
        subplot(3,1,3)
        plot(Data.Instants,Temp.PosteriorMeans(9,:))
        set(gca,'XTick',[0:(Data.Instants(end))/5:(Data.Instants(end))])
        set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
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
elseif strcmp(Parameters.Problem,'MarcFlu')
% for SEIR
    try
        subplot(2,1,1)
        plot(Data.Instants,(Temp.PosteriorMeans(5,:)))
        hold on
        plot(Data.Instants,Data.Observations(5,:),'g')
        plot(Data.Instants,Temp.Posterior975(5,:),'r')
        plot(Data.Instants,Temp.Posterior025(5,:),'r')
        hold off
        sizeStep = Data.NbComputingSteps(2);
        set(gca,'XTick',[1:6:36]*sizeStep)
        set(gca,'XTickLabel',{'1 Jun';'13 Jul';'25 Aug';'7 Oct';'19 Nov';'31 Dec';'5 Feb'})
        subplot(2,1,2)
        plot(Data.Instants,exp(Temp.PosteriorMeans(6,:)))
        set(gca,'XTick',[1:6:36]*sizeStep)
        set(gca,'XTickLabel',{'1 Jun';'13 Jul';'25 Aug';'7 Oct';'19 Nov';'31 Dec';'5 Feb'})
        title(['LogLik: ' num2str(Temp.LogLik)  ' - LogPrior: ' num2str(Temp.LogPrior)])
        subplot(2,1,1)
        plot(log(Temp.Liks))
        pause(0.001)
    end
elseif strcmp(Parameters.Problem,'MarcFluPlusObs')
% for SEIR
    try
        subplot(3,1,1)
        plot(Data.Instants,(Temp.PosteriorMeans(5,:)))
        hold on
        plot(Data.Instants,Data.Observations(5,:),'g')
        hold off
        sizeStep = Data.NbComputingSteps(2);
        set(gca,'XTick',[1:6:36]*sizeStep)
        set(gca,'XTickLabel',{'1 Jun';'13 Jul';'25 Aug';'7 Oct';'19 Nov';'31 Dec';'5 Feb'})
        subplot(3,1,2)
        plot(Data.Instants,exp(Temp.PosteriorMeans(6,:)))
        set(gca,'XTick',[1:6:36]*sizeStep)
        set(gca,'XTickLabel',{'1 Jun';'13 Jul';'25 Aug';'7 Oct';'19 Nov';'31 Dec';'5 Feb'})
        title(Temp.LogLik)
        subplot(3,1,3)
        plot(Data.Instants,exp(Temp.PosteriorMeans(7,:)))
        set(gca,'XTick',[1:6:36]*sizeStep)
        set(gca,'XTickLabel',{'1 Jun';'13 Jul';'25 Aug';'7 Oct';'19 Nov';'31 Dec';'5 Feb'})
        title(Temp.LogLik)
        pause(0.001)
    end
end
