% Cluster outputs

SavePath = '/Users/dureaujoseph/Documents/PhD_Data/fBM/';

SimSeries = 'TestingHidentifiability';

for i = 1:10
    load([SavePath '/' SimSeries '_' num2str(i)],'Res')
    PlotfBMoutput(Res)
    pause()
end

