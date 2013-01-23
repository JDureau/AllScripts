% Cluster outputs

SavePath = '/Users/dureaujoseph/Documents/PhD_Data/fBM/';

SimSeries = 'TestingHidentifiability';

for i = 1:10
    load([SavePath '/' SimSeries '_' num2str(i)],'Res')
%     plot(Res.Data.Y)
%     hold on
    PlotfBMoutput(Res)
    pause()
end
% hold off



files = {'15sep08_15sep09.csv','15mar07_15mar08.csv'};
DataSeries = 'SP500';
load([SavePath '/' DataSeries '_' num2str(2)])
plot(Res.Data.Y)
PlotfBMoutput(Res)
