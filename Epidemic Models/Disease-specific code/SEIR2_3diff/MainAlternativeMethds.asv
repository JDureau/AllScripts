DataPath = 'H:\My Documents\PhD Work\Data\HPA';

A = load([DataPath '\andre_estimates_31_01.txt']);

Data.Dates = {};
Data.NewCases = {};
for i = 1:size(A,1)
    Data.Dates{i} = i;
    for j = 1:7
    Data.NewCases{i}{j} = A(i,j)*10;
    end
end

InitialDate = struct();
InitialDate.Month = 6;
InitialDate.Day = 1;
InitialDate.Year = 2009;
Data.Dates = ApproxWeeklyDates(InitialDate,35);

plot(A)
legend('0-4','5-14','15-24','25-44','45-64','65+')

PopWeigths = [667600,2461800,5904100,6862500,14417400,12847800,7929300];
PopProps = PopWeigths/sum(PopWeigths)*100;
Weigthed = sum((A(:,1:7)*diag(PopWeigths.^-1)*diag(PopProps/100)*100000)');

plot(A*diag(PopWeigths.^-1)*100000)
hold on
plot(Weigthed,'k','LineWidth',2)
hold off
legend('<1','1-4','5-14','15-24','25-44','45-64','65+')




%% Inference
cd('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Filtering')
addpath('H:\My Documents\PhD Work\Matlab Scripts\General Tools')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes\Resampling\pf_resampling')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Joint Sampling')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\MIF')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Model Selection')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Optimization Approach')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\Models')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\SEIR')


SavePath = 'S:\Results\';

Res = load([SavePath 'Marc_FirstEst_WithPaths.mat']);

ResPost = Res.Res3

xt = zeros(5000*713,1);
xtp1 = zeros(5000*713,1);
for i = 1:5000
    xt((i-1)*713+1:(i)*713) = squeeze(exp(Res.Paths(i,6,1:713)));
    xtp1((i-1)*713+1:(i)*713) = squeeze(exp(Res.Paths(i,6,2:714)));
end
plot(xt,xt+1,'.')

xt = squeeze(exp(Res.Paths(:,6,:)));
Aggr = [squeeze(exp(Res.Paths(:,6,:))) Res.Thetas'];
cor = corr(Aggr*diag(std(Aggr).^-1));

k = length(names)
ind = 1;
indplot = 2
indplot2 = 4
for i = 1:k-1
    for j = i+1:k
        subplot(ceil(sqrt(k*(k-1)/2)),ceil(sqrt(k*(k-1)/2)),ind)
        ind = ind+1;
        scattercloud(ResPost.Thetas(i,:),ResPost.Thetas(j,:))
        xlabel(names{i})
        ylabel(names{j})
        hold on
        plot(ResPost.Thetas(i,inds(indplot)),ResPost.Thetas(j,inds(indplot)),'or')
        plot(ResPost.Thetas(i,inds(indplot2)),ResPost.Thetas(j,inds(indplot2)),'og')
        hold off
    end
end

SEIRModel = struct();
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;


Parameters = ResPost.Parameters;
[bof,inds] = sort(unique(ResPost.LogPosts),'descend');
LogPosts = [];
for j = 1:30
    disp(j)
    ind = inds(j);
    names = Parameters.Names.Estimated;
    for i = 1:length(names)
        indpar = Parameters.(names{i}).Index;
        Parameters.(names{i}).Value = Res.Thetas(indpar,ind);
    end
    Parameters = UpdateParsNoTransfToTransf(Parameters);
    for i = 1:15
        Temp = EstimationSMCsmoothGen(Res.Data, SEIRModel, Parameters);
        LogPosts(j,i) = Temp.LogLik*ResPost.LogPosts(ind)/ResPost.LogLiks(ind);
    end
end
[bof,inds2] = sort(median(LogPosts'));
ind = 2;
names = Parameters.Names.Estimated;
for i = 1:length(names)
    indpar = Parameters.(names{i}).Index;
    Parameters.(names{i}).Value = median(Res.Thetas(indpar,:));
%     Parameters.(names{i}).Value = Res.Thetas(indpar,ind);
end
Parameters = UpdateParsNoTransfToTransf(Parameters);




ResFilt = EstimationSMCfiltGen(Res.Data, SEIRModel, Parameters)


PlotMarc(ResFilt,8)

NbIts = 100
Thetas = zeros(7,NbIts);
Paths = zeros(NbIts,6,714);
LogLiks = [];
for i = 1:length(names)
    Thetas(i,:) = repmat(Parameters.(names{i}).Value,1,NbIts);
end
for i = 1:NbIts
    disp(i)
    Temp = EstimationSMCsmoothGen(Res.Data, SEIRModel, Parameters);
    ind = ceil(rand(1,1)*Parameters.NbParticules);
    Paths(i,:,:) = Temp.CompletePaths(ind,:,:);
    LogLiks(i) = Temp.LogLik;
end
ResSmooth = Res;
ResSmooth.Thetas = Thetas;
ResSmooth.Paths = Paths;
PlotMarc(ResSmooth,8) 
   
PlotMarc(ResSmooth,8)


%% plots 
clf
Resol = 8;
Data = Res.Data;
PathsInstant = 0:Res.Parameters.ComputationTStep:sum(Data.NbComputingSteps); 
Parameters = Res.Parameters;

dates = {};
delta = floor(length(Res.Data.Dates)/Resol);
Resol = floor(length(Res.Data.Dates)/delta);
inds = delta:delta:delta*Resol;
for i = 1:length(inds)
    dates{i} = [num2str(Res.Data.Dates{inds(i)}.Day) ' ' Res.Data.Dates{inds(i)}.MonthInLetters];
end

% % 1: filter
% clf
% Res = ResFilt;
% Paths = Res.Paths;
% subplot(3,1,1)
% temp = diag(Res.Thetas(Parameters.gamma.Index,:).^-1)*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
% ciplot(quantile(temp,0.025),quantile(temp,0.975),[172,215,255]/255)
% hold on
% ciplot(quantile(temp,0.25),quantile(temp,0.75),[100,153,251]/255)
% plot(mean(temp),'k','LineWidth',2)
% t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
% plot(t,1*ones(size(Data.Instants)),'--k','LineWidth',2)
% try
%     plot(Res.RtPath,'g','LineWidth',2) 
% end
% hold off
% xlim([0 Data.Instants(end)])
% TicksInds = Data.Instants(inds);
% set(gca,'XTick',TicksInds)
% %     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
% %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
% set(gca,'XTickLabel',dates)
% title('filter conditional distribution: p(R_{t}|\theta^{*},y_{1:floor(t)})')%|\theta^{*},y_{1:\left\floor t \right\floor})')
% ymax = 2.5;
% ylim([0.5 ymax])
% hold on
% yis1 = ymax/100:ymax/100:ymax*4/6;
% yis2 = ymax*5/6:ymax/100:ymax;
% plot(Data.Instants(8)*ones(size(yis1)),yis1,'r')
% plot(Data.Instants(8)*ones(size(yis2)),yis2,'r')
% text(Data.Instants(8),ymax*9/12,'sch. closure','HorizontalAlignment','center')
% plot(Data.Instants(14)*ones(size(yis1)),yis1,'r')
% plot(Data.Instants(14)*ones(size(yis2)),yis2,'r')
% text(Data.Instants(14),ymax*9/12,'end of holidays','HorizontalAlignment','center')
% hold off


% 2: smoothed filter
Res = ResFilt;
Paths = Res.Paths;
subplot(3,1,1)
temp = diag(Res.Thetas(Parameters.gamma.Index,:).^-1)*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
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
title('a) smoothed filter conditional distribution p(R_{t}|\theta^{*},y_{1:floor(t)})')%|\theta^{*},y_{1:\left\floor t \right\floor})')
ymax = 2.5;
ylim([0.5 ymax])
hold on
yis1 = ymax/100:ymax/100:ymax*4/6;
yis2 = ymax*5/6:ymax/100:ymax;
plot(8*ones(size(yis1)),yis1,'r')
plot(8*ones(size(yis2)),yis2,'r')
text(8,ymax*9/12,'sch. closure','HorizontalAlignment','center')
plot(14*ones(size(yis1)),yis1,'r')
plot(14*ones(size(yis2)),yis2,'r')
text(14,ymax*9/12,'end of holidays','HorizontalAlignment','center')
hold off

% 3: Smoothing distr
Res = ResSmooth;
Paths = Res.Paths;
subplot(3,1,2)
temp = diag(Res.Thetas(Parameters.gamma.Index,:).^-1)*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
ciplot(quantile(temp,0.025),quantile(temp,0.975),[172,215,255]/255)
hold on
ciplot(quantile(temp,0.25),quantile(temp,0.75),[100,153,251]/255)
plot(mean(temp),'k','LineWidth',2)
t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
plot(t,1*ones(size(Data.Instants)),'--k','LineWidth',2)
try
    plot(Res.RtPath,'g','LineWidth',2) 
end
hold off
xlim([0 Data.Instants(end)])
TicksInds = Data.Instants(inds);
set(gca,'XTick',TicksInds)
%     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
%     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
set(gca,'XTickLabel',dates)
title('b) smoother conditional distribution: p(R_{t}|\theta^{*},y_{1:n})')%|\theta^{*},y_{1:\left\floor t \right\floor})')
ymax = 2.5;
ylim([0.5 ymax])
hold on
yis1 = ymax/100:ymax/100:ymax*4/6;
yis2 = ymax*5/6:ymax/100:ymax;
plot(Data.Instants(8)*ones(size(yis1)),yis1,'r')
plot(Data.Instants(8)*ones(size(yis2)),yis2,'r')
text(Data.Instants(8),ymax*9/12,'sch. closure','HorizontalAlignment','center')
plot(Data.Instants(14)*ones(size(yis1)),yis1,'r')
plot(Data.Instants(14)*ones(size(yis2)),yis2,'r')
text(Data.Instants(14),ymax*9/12,'end of holidays','HorizontalAlignment','center')
hold off

% 4: Posterior marginal
Res = ResPost;
Paths = Res.Paths;
subplot(3,1,3)
temp = diag(Res.Thetas(Parameters.gamma.Index,:).^-1)*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
ciplot(quantile(temp,0.025),quantile(temp,0.975),[172,215,255]/255)
hold on
ciplot(quantile(temp,0.25),quantile(temp,0.75),[100,153,251]/255)
plot(mean(temp),'k','LineWidth',2)
t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
plot(t,1*ones(size(Data.Instants)),'--k','LineWidth',2)
try
    plot(Res.RtPath,'g','LineWidth',2) 
end
hold off
xlim([0 Data.Instants(end)])
TicksInds = Data.Instants(inds);
set(gca,'XTick',TicksInds)
%     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
%     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
set(gca,'XTickLabel',dates)
title('c) marginal posterior distribution: p(R_{t}|y_{1:n})')%|\theta^{*},y_{1:\left\floor t \right\floor})')
ymax = 2.5;
ylim([0.5 ymax])
hold on
yis1 = ymax/100:ymax/100:ymax*4/6;
yis2 = ymax*5/6:ymax/100:ymax;
plot(Data.Instants(8)*ones(size(yis1)),yis1,'r')
plot(Data.Instants(8)*ones(size(yis2)),yis2,'r')
text(Data.Instants(8),ymax*9/12,'sch. closure','HorizontalAlignment','center')
plot(Data.Instants(14)*ones(size(yis1)),yis1,'r')
plot(Data.Instants(14)*ones(size(yis2)),yis2,'r')
text(Data.Instants(14),ymax*9/12,'end of holidays','HorizontalAlignment','center')
hold off
