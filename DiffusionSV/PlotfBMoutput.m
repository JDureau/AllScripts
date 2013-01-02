function [] = PlotfBMoutput(Res)

Data = Res.Data;

step = Data.step;
nobs = Data.nobs;
obsstep = Data.obsstep;
N = Data.N;

subplot(4,2,1:2)
xis = step:step:nobs-1-step;
plot(Data.X(obsstep:obsstep:N-1),'g')
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
ylabel('Volatility','FontSize',20)
xlabel('Time (1 obs per time unit)','FontSize',20)

subplot(4,2,3)
plot(Res.out_Hs);
ylabel('H traceplot','FontSize',20)

subplot(4,2,4)
[fi,xi]=ksdensity(Res.out_Hs);
plot(xi,fi)
hold on
plot(Data.Htrue,0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
hold off
ylabel('H posterior estimate','FontSize',20)


subplot(4,2,5)
plot(Res.out_sigs);
ylabel('sigma_X traceplot','FontSize',20)

subplot(4,2,6)
[fi,xi]=ksdensity(Res.out_sigs);
plot(xi,fi)
hold on
plot(Data.sigma_Xtrue,0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
hold off
ylabel('sigma_X posterior estimate','FontSize',20)


n = size(Res.out_Xs,2);
subplot(4,2,7)
plot(Res.out_Xs(:,floor(n/3)));
ylabel('X(N/3) traceplot','FontSize',20)

subplot(4,2,8)
plot(Res.out_Xs(:,floor(2*n/3)));
ylabel('X(2*N/3) traceplot','FontSize',20)


clf
subplot(3,2,[1 2])
xis = step:step:nobs-1-step;
plot(Data.X(obsstep:obsstep:N-1),'g')
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
ylabel('Volatility','FontSize',20)
xlabel('Time (1 obs per time unit)','FontSize',20)

subplot(3,2,[3 5])
plot(Res.out_Hs,Res.out_sigs,'.')
xlabel('H','FontSize',20)
ylabel('sigma_X','FontSize',20)

subplot(3,2,4)
[fi,xi]=ksdensity(Res.out_sigs);
plot(xi,fi)
hold on
plot(Data.sigma_Xtrue,0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
hold off
ylabel('sigma_X posterior estimate','FontSize',20)

subplot(3,2,6)
[fi,xi]=ksdensity(Res.out_Hs);
plot(xi,fi)
hold on
plot(Data.Htrue,0,'og','MarkerEdgeColor','g','MarkerFaceColor','g')
hold off
ylabel('H posterior estimate','FontSize',20)

% subplot(3,2,5)
% plot(Res.out_sigs);
% subplot(3,2,6)
% [fi,xi]=ksdensity(Res.out_sigs);
% ylabel('H traceplot','FontSize',20)
% plot(xi,fi)
% hold on
% plot(Data.sigma_Xtrue,0,'og')
% hold off