% Cluster outputs

SavePath = '/Users/dureaujoseph/Documents/PhD_Data/fBM/';

SimSeries = 'TestingHidentifiability';

DataSet = 1;


hs = [];
clf
Ress = {};
for ind = 1:7
    ind
    load([SavePath '/D' num2str(DataSet) '_Exp' num2str(ind)])
    Ress{ind} = Res;
%     hs(ind,:) = Res.Thetas(1,:);
%     [fi,xi] = ksdensity(Res.Thetas(1,:));
%     plot(xi,fi)
%     hold on
%     plot(Res.Data.Y)
%     hold on
    PlotfBMoutput(Res)
%     disp(length(unique(Res.out_Ls))/150);
    pause()
end
% hold off




ESSs = {};

for k = 1:7
    k
    Res = Ress{k};
%     Ress{k}.h
    %figure(11); for i=10:10:min(500,loop); plot(out_Q(i,:),'b');hold on; end; plot(out_Q(1,:),'r');hold off
    try
        paths = Res.out_Zs;
    catch
        paths = squeeze(Res.Paths);
    end
    ESS=zeros(size(paths,2),1);
    for i=1:size(paths,2)
        r=sum(autocorr(paths(:,i),100));
        ESS(i)=100/(1+2*r);
    end
    ESSs{k} = ESS;
    
    figure(12);plot(ESS);title('ESS (%) for each v over time')   
    disp([min(ESS),median(ESS),max(ESS)])   
    
    Names = Res.Par.Names.Estimated;
    for i = 1:length(Names)
        ind = Res.Par.(Names{i}).Index;
        r=sum(autocorr(Ress{k}.Thetas(ind,:),1200));
        disp([Names{i} '   ' num2str(100/(1+2*r))]);
    end
%     r=sum(autocorr(Ress{k}.out_Hs,1200));
%     disp(100/(1+2*r));
%     r=sum(autocorr(Ress{k}.out_sigs,1200));
%     disp(100/(1+2*r));
end
    
plot(ESSs{1},'k')
hold on
plot(ESSs{2},'g')
plot(ESSs{3},'b')
hold off
ylabel('ESS(Z_t)','FontSize',20)
xlabel('t','FontSize',20)
h = legend('Gibbs-RW: ESS_{min}=0.1%','Joint-MALA: ESS_{min}=0.1%','Joint-HMC: ESS_{min}=1.1%');
set(h,'FontSize',20)




Par = Res.Par;
Names = Par.Names.Estimated;
for i = 1:length(Names)
    ind = Par.(Names{i}).Index;
    Par.(Names{i}).Value = Res.Thetas(ind,20000); 
end
Par = NoTransfToTransf(Par);

Data = Res.Data;
Res2 = RunJointMCMC_Full(Data,Par);


files = {'15sep08_15sep09.csv','15mar07_15mar08.csv'};
DataSeries = 'SP500';
load([SavePath '/Data_' DataSeries '_' num2str(2) '_' num2str(1) '.mat'])
DataTrue = Data;
Par = DataTrue.ParTrue;
Z = DataTrue.Z;

Bh = Z_to_Bh(Z,N,step,Par);
X = Bh_to_X_Full(Bh,step,Par);
Y = SampleObs_Full(X,Bh,step,Vol,Par);
BhTrue = Z_to_Bh(DataTrue.Z,N,step,Par);

clf
subplot(3,1,1)
plot(DataTrue.Y-DataTrue.Y(1),'g')
hold on
plot(Y)
hold off

subplot(3,1,2)
plot(cumsum(BhTrue),'g')
hold on
plot(cumsum(Bh))
hold off

subplot(3,1,3)
plot(cumsum(BhTrue),'g')
hold on
plot(cumsum(Bh))
hold off



PlotfBMoutput(Res)

Par = 
Par.nsteps = 10;
Par.loop = 200;
Par.h = 0.001;
Res = RunJointMCMC_Full(Data,Par);




