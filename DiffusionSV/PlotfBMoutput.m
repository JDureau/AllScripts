function [] = PlotfBMoutput(Res)

Data = Res.Data;

step = Data.step;
nobs = Data.nobs;
obsstep = Data.obsstep;
N = Data.N;
k = length(Res.Par.Names.Estimated);

figure(1)
clf
if k>1
    subplot(k+1,k,1:k)
end    
xis = step:step:nobs-1-step;
try
    plot(Data.X(obsstep:obsstep:N-1),'g')
catch
    q50 = quantile(Res.out_Xs,0.5);
end    
hold on
q2p5 = quantile(Res.out_Xs,0.025);
q25 = quantile(Res.out_Xs,0.25);
q50 = quantile(Res.out_Xs,0.5);
q75 = quantile(Res.out_Xs,0.75);
q97p5 = quantile(Res.out_Xs,0.975);
plot(q2p5,':')
plot(q25,'--')
plot(q50,'-')
plot(q75,'--')
plot(q97p5,':')
legend('True volatility','95% c.i.','50% c.i.','Posterior median')
hold off
ylabel('Volatility','FontSize',12)
xlabel('Time','FontSize',12)


Names = Res.Par.Names.Estimated;
if k>1
    for i = 1:k
        for j = 1:k
            h = subplot(k+1,k,k+(i-1)*k+j);
            if i==j
                ind = Res.Par.(Names{i}).Index;
                [fi,xi]=ksdensity(Res.Thetas(ind,:));
                plot(xi,fi)
                hold on
                try
                    plot(Data.ParTrue.(Names{i}).Value,0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
                end
                hold off
                ylabel(Names{i},'FontSize',12)
            elseif i>j
                indi = Res.Par.(Names{i}).Index;
                indj = Res.Par.(Names{j}).Index;
                plot(Res.Thetas(indj,:),Res.Thetas(indi,:),'.')
                xlabel(Names{j},'FontSize',12)
                ylabel(Names{i},'FontSize',12)
            else
                indi = Res.Par.(Names{i}).Index;
                indj = Res.Par.(Names{j}).Index;
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
    ind = Res.Par.(Names{i}).Index;
    plot(Res.Thetas(ind,:));
    title(Names{i})
end

subplot(ceil(sqrt(k+3)),ceil(sqrt(k+3)),k+1)
plot(Res.out_Ls)
% ylim([-500 -200])
title('LogLiks')
tmp = 1;
ind = 0;
while tmp
    subplot(ceil(sqrt(k+3)),ceil(sqrt(k+3)),k+2)
    ind = ind+1;
    plot(Res.out_Lpriorthetas(ind:end))
    tmp = 0;
end
% ylim([-500 -200])
title('LogPriorthetas')

tmp = 1;
ind = 0;
while tmp
    subplot(ceil(sqrt(k+3)),ceil(sqrt(k+3)),k+3)
    ind = ind+1;
    plot(Res.out_LpriorZ(ind:end))
    tmp = 0;
end
% ylim([-500 -200])
title('LogPriorZ')
    
%
% 
% subplot(4,2,3)
% plot(Res.out_Hs);
% ylabel('H traceplot','FontSize',20)
% 
% subplot(4,2,4)
% [fi,xi]=ksdensity(Res.out_Hs);
% plot(xi,fi)
% hold on
% plot(Data.Htrue,0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
% hold off
% ylabel('H posterior estimate','FontSize',20)
% 
% 
% subplot(4,2,5)
% plot(Res.out_sigs);
% ylabel('sigma_X traceplot','FontSize',20)
% 
% subplot(4,2,6)
% [fi,xi]=ksdensity(Res.out_sigs);
% plot(xi,fi)
% hold on
% plot(Data.sigma_Xtrue,0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
% hold off
% ylabel('sigma_X posterior estimate','FontSize',20)
% 
% 
% n = size(Res.out_Xs,2);
% subplot(4,2,7)
% plot(Res.out_Xs(:,floor(n/3)));
% ylabel('X(N/3) traceplot','FontSize',20)
% 
% subplot(4,2,8)
% plot(Res.out_Xs(:,floor(2*n/3)));
% ylabel('X(2*N/3) traceplot','FontSize',20)
% 
% 
% clf
% subplot(3,2,[1 2])
% xis = step:step:nobs-1-step;
% plot(Data.X(obsstep:obsstep:N-1),'g')
% hold on
% q2p5 = quantile(Res.out_Xs,0.025);
% q25 = quantile(Res.out_Xs,0.25);
% q50 = quantile(Res.out_Xs,0.5);
% q75 = quantile(Res.out_Xs,0.75);
% q97p5 = quantile(Res.out_Xs,0.975);
% plot(q2p5,':')
% plot(q25,'--')
% plot(q50,'-')
% plot(q75,'--')
% plot(q97p5,':')
% legend('True volatility','95% c.i.','50% c.i.','Posterior median')
% hold off
% ylabel('Volatility','FontSize',20)
% xlabel('Time (1 obs per time unit)','FontSize',20)
% 
% subplot(3,2,[3 5])
% plot(Res.out_Hs,Res.out_sigs,'.')
% xlabel('H','FontSize',20)
% ylabel('sigma_X','FontSize',20)
% 
% subplot(3,2,4)
% [fi,xi]=ksdensity(Res.out_sigs);
% plot(xi,fi)
% hold on
% plot(Data.sigma_Xtrue,0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
% hold off
% ylabel('sigma_X posterior estimate','FontSize',20)
% 
% subplot(3,2,6)
% [fi,xi]=ksdensity(Res.out_Hs);
% plot(xi,fi)
% hold on
% plot(Data.Htrue,0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
% hold off
% ylabel('H posterior estimate','FontSize',20)
