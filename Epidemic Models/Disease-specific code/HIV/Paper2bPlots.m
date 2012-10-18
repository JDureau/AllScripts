% Paper 2b plots

SavePath = '/Users/dureaujoseph/Documents/PhD_Data/Avahan/';


Districts = {'Mysore_3rounds','Belgaum_3rounds','Bellary','EastGodavry','Guntur','Hyderabad','Shimoga','Yevatmal'};


% BM
for i = 1:length(Districts)
    tmp = [SavePath 'HIV_' Districts{i}  '.mat'];
    load(tmp)
    PlotResHIV(Res,Res.Parameters);
    Districts{i} 
    pause()
end

% Sigm
for i = 1:length(Districts)
    tmp = [SavePath 'HIV_' Districts{i}  '_Sigm.mat'];
    load(tmp)
    PlotResHIV(Res,Res.Parameters);
    Districts{i} 
    pause()
end
