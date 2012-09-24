% Results analysis
cd('/Users/dureaujoseph/Dropbox/AllScriptsGit/')
addpath([pwd '/MCMCTests/'])
addpath([pwd '/General Tools/'])
addpath([pwd '/Toolboxes'])


SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/MCMCstudy/';

%% STEP 1: check if they all get to the same posteriors

Methods =  {'MALA','GMCovMALA','GMRand','GMLang','HMCGMCov','GMind'};
Densities = {'GMM','GMM2','Banana'};

IndDensity = 3;
dim = 2;

cols = {'k','c','y','b','r','g','c','k'};

clf
for IndMethod = 1:8
        clf
    title(Methods{IndMethod})
    for IndLogOrNot = 0:0
        for indeps = 1:10
            Epss = (1:14)*0.2;
            Parameters.Epsil = Epss(indeps);
            if and(IndMethod>=5,IndMethod<=8)
                Parameters.Epsil = Parameters.Epsil/10;
            end
            try
                disp([Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(indeps) '.mat'])
                load([SavePath Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(indeps) '.mat']);
                Parameters.Epsil = Res.Parameters.Epsil;
                for  i = 1:dim
                    [fi,xi] = ksdensity(Res.Vals(i,:));
                    subplot(dim+1,1,i)
                    hold on
                    plot(xi,fi,cols{IndMethod})
                    hold off
                end
                subplot(dim+1,1,dim+1)
                scattercloud(Res.Vals(1,:),Res.Vals(2,:))
                pause(0.001)
            end
        end
    end
    pause()
end

%% STEP 2: check ESS plots

cols = {'k',':k','y','b','r','--r',':r','g'};
clf
hold on
for IndMethod = 1:8
    for IndLogOrNot = 0:1
        tmpESS = [];
        for indeps = 1:14
            Epss = (1:14)*0.2;
            Parameters.Epsil = Epss(indeps);
            if and(IndMethod>=5,IndMethod<=8)
                Parameters.Epsil = Parameters.Epsil/10;
            end
            try
                disp([Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat'])
                load([SavePath Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat']);
                tmpESS(indeps) = min(Res.RelESS);
            end
        end
        try
            plot(Epss,tmpESS,cols{IndMethod})
        end
        
%         plot(Epss,tmpESS,cols{IndMethod})
%         max(tmpESS)
    end
end
hold off
% ylim([0 80])

clf
for i = 1:2
    subplot(2,1,i)
    plot(Res.Vals(i,:))
end

