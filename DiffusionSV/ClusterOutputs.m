% Cluster outputs

SavePath = '/Users/dureaujoseph/Documents/PhD_Data/fBM/';

SimSeries = 'TestingHidentifiability';

for i = 1:5
    load([SavePath '/Res_Hsims_' num2str(i)],'Res')
%     plot(Res.Data.Y)
%     hold on
    PlotfBMoutput(Res)
%     disp(length(unique(Res.out_Ls))/150);
    pause()
end
% hold off



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




