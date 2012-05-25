
SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';
ListESSs = [];
inds = [];
for i = 1:100
    try load([SavePath '/VIH_PlayingWithSigm_ResultsBertallanfy__1_' num2str(i) '.mat'])
        ResBer = Ress{1};
        TmpESSs = [];
        for j = 1:size(ResBer.Thetas,1)
              temp = AutoCorrelation(ResBer.Thetas(j,:),100);
              TmpESSs(j) = 1/(1+2*sum(temp));
        end
        ListESSs(end+1) = min(TmpESSs);
        inds(end+1) = i;
    end
end

[b,ind] = max(ListESSs);

load([SavePath '/VIH_PlayingWithSigm_ResultsBertallanfy__1_' num2str(inds(ind)) '.mat']);

Cov = cov(Ress{1}.TransfThetas);
save([SavePath '/BestBertCov.mat'],'Cov')




                