function [] =  PlotfBMpmcmc(Res) 


Data = Res.Data;
Par = Res.Parameters;

step = Data.step;
nobs = Data.nobs;
obsstep = Data.obsstep;
N = Data.N;
k = length(Par.Names.Estimated);

figure(1)
clf
if k>1
    subplot(k+1,k,1:k)
end    
xis = step:step:nobs-1-step;
try
    plot(Data.X(obsstep:obsstep:N-1),'g')
catch
    q50 = quantile(squeeze(Res.Paths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.5);
end    
hold on
try
    q2p5 = quantile(squeeze(Res.Paths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.025);
    q25 = quantile(squeeze(Res.Paths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.25);
    q50 = quantile(squeeze(Res.Paths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.5);
    q75 = quantile(squeeze(Res.Paths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.75);
    q97p5 = quantile(squeeze(Res.Paths(:,1,cumsum(Res.Data.NbComputingSteps)+1)),0.975);
    plot(q2p5,':')
    plot(q25,'--')
    plot(q50,'-')
    plot(q75,'--')
    plot(q97p5,':')
end
legend('True volatility','95% c.i.','50% c.i.','Posterior median')
hold off
ylabel('Volatility','FontSize',12)
xlabel('Time','FontSize',12)


Names = Par.Names.Estimated;
if k>1
    for i = 1:k
        for j = 1:k
            h = subplot(k+1,k,k+(i-1)*k+j);
            if i==j
                ind = Par.(Names{i}).Index;
                [fi,xi]=ksdensity(Res.Thetas(ind,:));
                plot(xi,fi)
                hold on
                try
                    plot(Data.ParTrue.(Names{i}).Value,0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
                end
                hold off
                ylabel(Names{i},'FontSize',12)
            elseif i>j
                indi = Par.(Names{i}).Index;
                indj = Par.(Names{j}).Index;
                plot(Res.Thetas(indj,:),Res.Thetas(indi,:),'.')
                xlabel(Names{j},'FontSize',12)
                ylabel(Names{i},'FontSize',12)
            else
                indi = Par.(Names{i}).Index;
                indj = Par.(Names{j}).Index;
                text(0.5,0.5,num2str(corr(Res.Thetas(indi,:)',Res.Thetas(indj,:)'),2),'FontSize',20);
                set ( h, 'visible', 'off')
            end
        end
    end
end

figure(2)
clf
for i = 1:k
    subplot(ceil(sqrt(k+3)),ceil(sqrt(k+3)),i)
    ind = Par.(Names{i}).Index;
    plot(Res.Thetas(ind,:));
    title(Names{i})
end

subplot(ceil(sqrt(k+3)),ceil(sqrt(k+3)),k+1)
plot(Res.LogLiks)
% ylim([-500 -200])
title('LogLiks')

subplot(ceil(sqrt(k+3)),ceil(sqrt(k+3)),k+2)
plot(Res.LogPosts)
% ylim([-500 -200])
title('LogPosts')

