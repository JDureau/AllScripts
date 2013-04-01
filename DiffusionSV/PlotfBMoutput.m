function [] = PlotfBMoutput(Res)

Data = Res.Data;

step = Data.step;
nobs = Data.nobs;
obsstep = Data.obsstep;
N = Data.N;
k = length(Res.Par.Names.Estimated);
k2 = length(Res.Par.Names.Estimated);% + Res.Par.NbZpar;


try
   Res.Data.X;
   sim = 1;
catch
   sim = 0;
end
    

nbrow = ceil(sqrt(k2+5));
nbcol = ceil(sqrt(k2+5));
if k2+5 <=  (nbrow-1)*nbcol
    nbrow = nbrow-1;
end


figure(1)
clf
if k2>1
    subplot(k2+1,k2,1:k2)
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
% legend('True volatility','95% c.i.','50% c.i.','Posterior median')
hold off
ylabel('Volatility','FontSize',12)
xlabel('Time','FontSize',12)


Names = Res.Par.Names.Estimated;
if k2>=1
    for i = 1:k2
        for j = 1:k2
            h = subplot(k2+1,k2,k2+(i-1)*k2+j);
            if i==j
                if i>k
                    [fi,xi]=ksdensity(Res.TransfThetas(i,:));
                    plot(xi,fi)
                    hold on
                    try
%                         plot(Data.Z(i-k),0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
                    end
                    hold off
                    ylabel(['Z' num2str(i-k)],'FontSize',12)
                else
                    ind = Res.Par.(Names{i}).Index;
                    [fi,xi]=ksdensity(Res.Thetas(ind,:));
                    plot(xi,fi)
                    hold on
                    if sim
                        plot(Data.ParTrue.(Names{i}).Value,0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
                    end
                    hold off
                    ylabel(Names{i},'FontSize',12)
                end
            elseif i>j
                if and(i>k,j>k)
                    plot(Res.TransfThetas(j-k,:),Res.TransfThetas(i-k,:),'.')
                    ylabel(['Z' num2str(i-k)],'FontSize',12)
                    xlabel(['Z' num2str(j-k)],'FontSize',12)
                elseif i>k
                    indj = Res.Par.(Names{j}).Index;
                    plot(Res.Thetas(indj,:),Res.TransfThetas(i-k,:),'.')
                    xlabel(Names{j},'FontSize',12)
                    ylabel(['Z' num2str(i-k)],'FontSize',12)
                elseif j>k 
                    indi = Res.Par.(Names{i}).Index;
                    plot(Res.TransfThetas(j-k,:),Res.Thetas(indi,:),'.')
                    ylabel(Names{i},'FontSize',12)
                    xlabel(['Z' num2str(j-k)],'FontSize',12)
                else
                    indi = Res.Par.(Names{i}).Index;
                    indj = Res.Par.(Names{j}).Index;
                    plot(Res.Thetas(indj,:),Res.Thetas(indi,:),'.')
                    xlabel(Names{j},'FontSize',12)
                    ylabel(Names{i},'FontSize',12)
                end
            else
                if and(i>k,j>k)
                    indi = i;
                    indj = j;
                elseif i>k
                    indi = i;
                    indj = Res.Par.(Names{j}).Index;
                elseif j>k 
                    indi = Res.Par.(Names{i}).Index;
                    indj = j;
                else
                    indi = Res.Par.(Names{i}).Index;
                    indj = Res.Par.(Names{j}).Index;
                end
                try
                    text(0.5,0.5,num2str(corr(Res.TransfThetas(indi,:)',Res.TransfThetas(indj,:)'),2),'FontSize',20);
                catch
                    'y'
                end
                set ( h, 'visible', 'off')
            end
        end
    end
end

figure(2)
clf
for i = 1:k2
    subplot(nbrow,nbcol,i)
    if i >k
        plot(Res.TransfThetas(i,:),'k');
        title(['Z' num2str(i-k)])
    else
        ind = Res.Par.(Names{i}).Index;
        plot(Res.Thetas(ind,:),'k');
        title(Names{i})
    end
end

subplot(nbrow,nbcol,k2+1)
plot(Res.out_Ls,'k')
% ylim([-500 -200])
title('LogLiks')
tmp = 1;
ind = 0;
while tmp
    subplot(nbrow,nbcol,k2+2)
    ind = ind+1;
    plot(Res.out_Lpriorthetas(ind:end),'k')
    tmp = 0;
end
% ylim([-500 -200])
title('LogPrior')

tmp = 1;
ind = 0;
while tmp
    subplot(nbrow,nbcol,k+3)
    ind = ind+1;
    plot(Res.out_LpriorZ,'k')
    tmp = 0;
end
% ylim([-500 -200])
title('LpriorZ')
    
try
    
    max1 = max(Res.out_LT1s-mean(Res.out_LT1s));
    max2 = max(Res.out_LT2s-mean(Res.out_LT2s(2:end)));
    min1 = min(Res.out_LT1s-mean(Res.out_LT1s));
    min2 = min(Res.out_LT2s-mean(Res.out_LT2s(2:end)));
    
    subplot(nbrow,nbcol,k2+4)
    plot(Res.out_LT1s,'k')
    title('variance term (logscale)')
    ylim([min(min1,min2) max(max1,max2)] + mean(Res.out_LT1s))
    
    subplot(nbrow,nbcol,k2+5)
%     plot((Res.out_LT1s-mean(Res.out_LT1s)))
%     hold on
    plot(Res.out_LT2s,'k')
    title('exponential term (logscale)')
    ylim([min(min1,min2) max(max1,max2)] + mean(Res.out_LT2s(2:end)))

%     plot(Res.out_Lpriorthetas-mean(Res.out_Lpriorthetas(2:end)),'g')
%     plot((2*Res.Data.LogTerm1-mean(Res.out_LT1s(2:end)))*ones(size(Res.out_LT2s)),'r')
%     hold off
%     legend('variance term','exp term','prior term','true variance term')
end
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
