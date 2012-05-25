function [] = PlotPostForAll(Res,Mode)


Stats = {'ampl','ampl. aft. 1995','ampl aft 2003','asympt'};
Tabs = {};



% tmpinds = find(not(isnan(Res.ampls(1,:))));
% Res.ampls = Res.ampls(:,tmpinds);
% Res.partampls1995 = Res.partampls1995(:,tmpinds);
% Res.partampls2003 = Res.partampls2003(:,tmpinds);
% Res.asympts = Res.asympts(:,tmpinds);
% tmp = find(not(Res.NbSamples(:,1)==0));
% Res.NbSamples = Res.NbSamples(tmp,:);
% inds = find(and(Res.NbSamples(:,2)>50, Res.NbSamples(:,3)>20));
% length(inds)
inds = 1:size(Res.ampls,2);


if strcmp(Mode,'ROC')
    Tabs{1} = Res.amplsROC(:,inds);
    Tabs{2} = Res.partampls1995ROC(:,inds);
    Tabs{3} = Res.partampls2003ROC(:,inds);
    Tabs{4} = Res.asymptsROC(:,inds);
else
    Tabs{1} = Res.ampls(:,inds);
    Tabs{2} = Res.partampls1995(:,inds);
    Tabs{3} = Res.partampls2003(:,inds);
    Tabs{4} = Res.asympts(:,inds);
end





Corrs = [];
Biases = [];
MSEs = [];
Vars = [];
ROCs = [];
inds = [1 3 2];
Pts = {};
for i = 1:4
    for j = 1:3
        Corrs(inds(j),i) = corr(Tabs{i}(4,:)',Tabs{i}(j,:)');
        Biases(inds(j),i) = mean(Tabs{i}(j,:)-Tabs{i}(4,:)); 
        MSEs(inds(j),i) = mean((Tabs{i}(4,:)-Tabs{i}(j,:)).^2);
        Vars(inds(j),i) = MSEs(inds(j),i)-Biases(inds(j),i).^2;
    end
end


AllROCs = [];
inds = [1 3 2];
quants = 0:0.2:1;
qs = [];
for i = 1:4
    for j = 1:3
        for k = 1:length(quants)
            q = quantile(Tabs{i}(4,:),quants(k));
            qs(i,k) = q;
    %         med = median(Tabs{i}(4,:));
            indsinf = find(Tabs{i}(4,:)<=q);
            indssup = find(Tabs{i}(4,:)>q);
            if isempty(indssup)
                indsinf = find(Tabs{i}(4,:)<q);
                indssup = find(Tabs{i}(4,:)>=q);
            end
            xs = [];
            ys = [];
            all = [Tabs{i}(1,:) Tabs{i}(2,:) Tabs{i}(3,:)];
            m = min(all);
            M = max(all);
    %         deltap = sort((all(1:end-1)+all(2:end))/2);
            deltap = sort(all);
            deltap = (deltap(1:end-1)+deltap(2:end))/2;
            deltap = [(min(all)-mean(diff(deltap))) deltap (max(all)+mean(diff(deltap)))];
            for l = 1:length(deltap)
                xs(l) = sum(Tabs{i}(j,indsinf)>deltap(l))/length(indsinf);
                ys(l) = sum(Tabs{i}(j,indssup)>deltap(l))/length(indssup);
            end
            ROC= 0;
            for l = 1:length(deltap)-1
                ROC = ROC + (xs(l)-xs(l+1))*(ys(l+1)+ys(l))/2;
            end
            AllROCs(inds(j),i,k) = ROC-0.5;
    %         ROCs(inds(j),i) = -diff(xs)*(ys(1:end-1)'+ys(2:end)')/2-0.5;
        end
    end
end

subplot(5,4,1:4)
bar(Corrs')
set(gca,'XTickLabel',Stats)
ylim([0 1])
ylabel('Correlation')

subplot(5,4,5:8)
bar(Biases')
set(gca,'XTickLabel',Stats)
ylim([-0.4 0.2])
ylabel('Biases')

subplot(5,4,9:12)
bar(MSEs')
set(gca,'XTickLabel',Stats)
ylim([0 0.2])
ylabel('MSEs')
legend('Det. BR','Sto. BR',' BM')

subplot(5,4,13:16)
bar(Vars')
set(gca,'XTickLabel',Stats)
ylabel('Vars')


for i = 1:4
   subplot(5,4,16+i) 
   plot(qs(i,:),squeeze(AllROCs(1,i,:)),'b')
   hold on
   plot(qs(i,:),squeeze(AllROCs(2,i,:)),'g')
   plot(qs(i,:),squeeze(AllROCs(3,i,:)),'r')
   hold off
   ylim([0 0.5])
end

ResStats = struct();
ResStats.Corrs = Corrs;
ResStats.Biases = Biases;
ResStats.MSEs = MSEs;
ResStats.Vars = Vars;
ResStats.ROCs = ROCs;
