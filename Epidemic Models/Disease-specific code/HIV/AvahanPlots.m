


SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan';


load([SavePath '/TreatedRessRightPriors.mat'])


clf


Green = [0 205 0]/255


load([SavePath '/TreatedRessRightPriors_Biased.mat'])

indsbig = find(Res.amplin>0.5)


corr(Res.inflptin',Res.inflptout')

corr(Res.inflptin(indsbig)',Res.inflptout(indsbig)')

corr(Res.amplin',Res.amplout')


clf

subplot(1,2,1)
load([SavePath '/TreatedRessRightPriors.mat'])
indsbig = find(Res.amplin>0.5);
indssmall = find(Res.amplin<0.5);
plot(Res.amplin(indsbig),Res.amplout(indsbig),'o','Color',Green,'MarkerFaceColor','g','MarkerSize',10)
hold on
plot(Res.amplin(indssmall),Res.amplout(indssmall),'o','Color',Green,'MarkerFaceColor','g','MarkerSize',4)
xis = 0:0.01:1;
plot(xis,xis,'k')
hold off
xlim([0 1])
ylim([0 1])
title('Present scenario','FontSize',14)
xlabel('Simulated CU amplitude','FontSize',12)
ylabel('Estimated CU amplitude','FontSize',12)

subplot(1,2,2)
load([SavePath '/TreatedRessRightPriors_Biased.mat'])
indsbig = find(Res.amplin>0.5);
indssmall = find(Res.amplin<0.5);
plot(Res.amplin(indsbig),Res.amplout(indsbig),'o','Color',Green,'MarkerFaceColor','g','MarkerSize',10)
hold on
plot(Res.amplin(indssmall),Res.amplout(indssmall),'o','Color',Green,'MarkerFaceColor','g','MarkerSize',4)
xis = 0:0.01:1;
plot(xis,xis,'k')
hold off
xlim([0 1])
ylim([0 1])
title('Negative scenario','FontSize',14)
xlabel('Simulated CU amplitude','FontSize',12)
ylabel('Estimated CU amplitude','FontSize',12)








subplot(1,2,1)
load([SavePath '/TreatedRessRightPriors.mat'])
indsbig = find(Res.amplin>0.5);
indssmall = find(Res.amplin<0.5);
plot(Res.inflptin(indsbig),Res.inflptout(indsbig),'o','Color',Green,'MarkerFaceColor','g','MarkerSize',10)
hold on
plot(Res.inflptin(indssmall),Res.inflptout(indssmall),'o','Color',Green,'MarkerFaceColor','g','MarkerSize',4)
xis = 0:0.01:1;
plot(xis,xis,'k')
hold off
set(gca,'YTick',[120:120:600])
set(gca,'YTickLabel',['1990';'1995';'2000';'2005';'2010'],'FontSize',14)
set(gca,'XTick',[120:120:600])
set(gca,'XTickLabel',['1990';'1995';'2000';'2005';'2010'],'FontSize',14)
xlim([120 600])
ylim([120 600])
title('Present scenario','FontSize',14)
xlabel('Simulated CU shift timing','FontSize',12)
ylabel('Estimated CU shift timing','FontSize',12)


subplot(1,2,2)
load([SavePath '/TreatedRessRightPriors_Biased.mat'])
indsbig = find(Res.amplin>0.5);
indssmall = find(Res.amplin<0.5);
plot(Res.inflptin(indsbig),Res.inflptout(indsbig),'o','Color',Green,'MarkerFaceColor','g','MarkerSize',10)
hold on
plot(Res.inflptin(indssmall),Res.inflptout(indssmall),'o','Color',Green,'MarkerFaceColor','g','MarkerSize',4)
xis = 0:0.01:1;
plot(xis,xis,'k')
hold off
set(gca,'YTick',[120:120:600])
set(gca,'YTickLabel',['1990';'1995';'2000';'2005';'2010'],'FontSize',14)
set(gca,'XTick',[120:120:600])
set(gca,'XTickLabel',['1990';'1995';'2000';'2005';'2010'],'FontSize',14)
xlim([120 600])
ylim([120 600])
title('Negative scenario','FontSize',14)
xlabel('Simulated CU shift timing','FontSize',12)
ylabel('Estimated CU shift timing','FontSize',12)



clf
plot(Res.inflptin,Res.inflptout,'.')



files = {'/TreatedRessRightPriorsLogisticWithCorrs.mat','/TreatedRessRightPriorsLogisticWithCorrs_Biased.mat','/TreatedRessRightPriorsWithCorrs_Restricted.mat','/TreatedRessRightPriorsWithCorrs_AddingObs.mat'};



subplot(1,2,1)
Positions = 1:4:9;
load([SavePath files{1}])
slopes = Res.slopes;
qs = [];
qs(1) = -10;
qs(2) = quantile(slopes,0.33);
qs(3) = quantile(slopes,0.66);
qs(4) = quantile(slopes,0.99);
boxes = {};
for i = 1:2
    load([SavePath files{i}])
    slopes = Res.slopes;
    diffs = Res.corramplout95-Res.corramplout05;
    boxes{i} = [];
    ls = [];
    for j = 1:3
        inds = find(and(slopes>qs(j),slopes<qs(j+1)));
        ls(j) = length(inds);
    end
    for j = 1:3
        inds = find(and(slopes>qs(j),slopes<qs(j+1)));
        boxes{i}(j,:) = diffs(inds(1:min(ls)));
    end
end


Colors = {[24;116;205]/255;[238;59;59]/255;[105;139;34]/255;[179;238;58]/255};
ColorsEdges = {};

for i = 1:2
    ColorsEdges{i} = Colors{i} - 20/255;
end
h = plot(100,100,'s','color',ColorsEdges{1},'MarkerFaceColor',Colors{1});
hold on
for i = 2:2
    plot(100,100,'s','color',ColorsEdges{i},'MarkerFaceColor',Colors{i});
end
h = boxplot(boxes{1}','colors',ColorsEdges{1},'positions',Positions','width',0.2,'symbol','ow');
h = h(5,:);
% h = findobj(gca,'Tag','Box');
set(gca,'XTickLabel',{''}); 
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),Colors{1}','FaceAlpha',.5);
end
for i = 2:2
   h = boxplot(boxes{i}','colors',ColorsEdges{i}','positions',(Positions+0.4*(i-1))','width',0.2,'symbol','ow')
   h = h(5,:);
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),Colors{i}','FaceAlpha',.2);
    end
    set(gca,'XTickLabel',{''}); 
end
hold off
ylim([0 1])
xlim([0 12])
text('Position',[1.6,-0.07],'HorizontalAlignment','center','String',{'Strong decreasing'; 'observed prevalence'},'FontSize',11) 
text('Position',[5.6,-0.07],'HorizontalAlignment','center','String',{'Moderately decreasing'; 'observed prevalence'},'FontSize',11)
text('Position',[9.6,-0.07],'HorizontalAlignment','center','String',{'Stable or slightly increasing'; 'observed prevalence'},'FontSize',11)
ylabel({'Width of the 90%','confidence interval     ','of \Delta estimates'},'FontSize',12)
% ylabel({'Width of the 95%','confidence interval','of \Delta estimates'},'FontSize',12)
title('To what level of accuracy can  \Delta be estimated?  ','FontSize',14)
legend('Assumed situation','Negative scenario')%,'Positive scenario 1','Positive scenario 1')
% legend boxoff 



subplot(1,2,2)

Positions = 1:4:9;

load([SavePath files{1}])
slopes = Res.slopes;
qs = [];
qs(1) = -10;
qs(2) = quantile(slopes,0.33);
qs(3) = quantile(slopes,0.66);
qs(4) = quantile(slopes,0.99);
boxes = {};
for i = 1:2
    load([SavePath files{i}])
    slopes = Res.slopes;
    diffs = Res.corrinflptout95-Res.corrinflptout05;
    boxes{i} = [];
    ls = [];
    for j = 1:3
        inds = find(and(slopes>qs(j),slopes<qs(j+1)));
        ls(j) = length(inds);
    end
    for j = 1:3
        inds = find(and(slopes>qs(j),slopes<qs(j+1)));
        boxes{i}(j,:) = diffs(inds(1:min(ls)));
    end
end



Colors = {[24;116;205]/255;[238;59;59]/255;[105;139;34]/255;[179;238;58]/255};
ColorsEdges = {};

for i = 1:2
    ColorsEdges{i} = Colors{i} - 20/255;
end
h = plot(100,100,'s','color',ColorsEdges{1},'MarkerFaceColor',Colors{1});
hold on
for i = 2:2
    plot(100,100,'s','color',ColorsEdges{i},'MarkerFaceColor',Colors{i});
end
h = boxplot(boxes{1}','colors',ColorsEdges{1},'positions',Positions','width',0.2,'symbol','ow');
h = h(5,:);
% h = findobj(gca,'Tag','Box');
set(gca,'XTickLabel',{''}); 
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),Colors{1}','FaceAlpha',.5);
end
for i = 2:2
   h = boxplot(boxes{i}','colors',ColorsEdges{i}','positions',(Positions+0.4*(i-1))','width',0.2,'symbol','ow')
   h = h(5,:);
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),Colors{i}','FaceAlpha',.2);
    end
    set(gca,'XTickLabel',{''}); 
end
hold off
ylim([0 600])
xlim([0 12])
set(gca,'YTick',0:96:600)
set(gca,'yTickLabel',{'0 yrs.','4 yrs.','8 yrs.','12 yrs.','16 yrs.','20 yrs.','24 yrs.'})
text('Position',[1.6,-50],'HorizontalAlignment','center','String',{'Strong decreasing'; 'observed prevalence'},'FontSize',11) 
text('Position',[5.6,-50],'HorizontalAlignment','center','String',{'Moderately decreasing'; 'observed prevalence'},'FontSize',11)
text('Position',[9.6,-50],'HorizontalAlignment','center','String',{'Stable or slightly increasing'; 'observed prevalence'},'FontSize',11)
ylabel({'Width of the 90%','confidence interval      ','of T estimates'},'FontSize',12)
% ylabel({'Width of the 95%','confidence interval','of T estimates'},'FontSize',12)
title('To what level of accuracy can  T be estimated?  ','FontSize',14)
legend('Assumed situation','Negative scenario')%,'Positive scenario 1','Positive scenario 1')
% legend boxoff 










Positions = 1:4:9;

load([SavePath files{1}])
slopes = Res.slopes;
qs = [];
qs(1) = -10;
qs(2) = quantile(slopes,0.33);
qs(3) = quantile(slopes,0.66);
qs(4) = quantile(slopes,0.99);
boxes = {};
for i = 1:4
    load([SavePath files{i}])
    slopes = Res.slopes;
    diffs = Res.corrpropaftintout975-Res.corrpropaftintout025;
    boxes{i} = [];
    ls = [];
    for j = 1:3
        inds = find(and(slopes>qs(j),slopes<qs(j+1)));
        ls(j) = length(inds);
    end
    for j = 1:3
        inds = find(and(slopes>qs(j),slopes<qs(j+1)));
        boxes{i}(j,:) = diffs(inds(1:min(ls)));
    end
end


Colors = {[24;116;205]/255;[238;59;59]/255;[105;139;34]/255;[179;238;58]/255};
ColorsEdges = {};

for i = 1:4
    ColorsEdges{i} = Colors{i} - 20/255;
end
h = plot(100,100,'s','color',ColorsEdges{1},'MarkerFaceColor',Colors{1});
hold on
for i = 2:4
    plot(100,100,'s','color',ColorsEdges{i},'MarkerFaceColor',Colors{i});
end
h = boxplot(boxes{1}','colors',ColorsEdges{1},'positions',Positions','width',0.2,'symbol','ow');
h = h(5,:);
% h = findobj(gca,'Tag','Box');
set(gca,'XTickLabel',{''}); 
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),Colors{1}','FaceAlpha',.5);
end
for i = 2:4
   h = boxplot(boxes{i}','colors',ColorsEdges{i}','positions',(Positions+0.4*(i-1))','width',0.2,'symbol','ow')
   h = h(5,:);
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),Colors{i}','FaceAlpha',.2);
    end
    set(gca,'XTickLabel',{''}); 
end
hold off
ylim([0 1])
xlim([0 12])
text('Position',[1.6,-0.07],'HorizontalAlignment','center','String',{'Strong decreasing'; 'observed prevalence'},'FontSize',11) 
text('Position',[5.6,-0.07],'HorizontalAlignment','center','String',{'Moderately decreasing'; 'observed prevalence'},'FontSize',11)
text('Position',[9.6,-0.07],'HorizontalAlignment','center','String',{'Stable or slightly increasing'; 'observed prevalence'},'FontSize',11)
ylabel({'Width of the 95%','confidence interval       ',''},'FontSize',12)
title({'To what level of accuracy can the proportion of the CU   ';'increase happening after 2004 be estimated?  '},'FontSize',14)
legend('Assumed situation','Negative scenario','Positive scenario 1','Positive scenario 1')
legend boxoff 



