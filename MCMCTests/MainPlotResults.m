% Results analysis
cd('/Users/dureaujoseph/Dropbox/AllScriptsGit/')
addpath([pwd '/MCMCTests/'])
addpath([pwd '/General Tools/'])
addpath([pwd '/Toolboxes'])


SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/MCMCstudy/';

%% STEP 1: check if they all get to the same posteriors

Methods =  {'MALA','GMCovMALA','GMRand','GMLang','HMCGMCov','GMind'};
Densities = {'GMM','GMM2','Banana15','Banana5','Banana1'};

IndDensity = 3;
dim = 2;

cols = {'k','c','y','b','r','g','c','k'};

clf
for IndMethod = 1:6
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
                load([SavePath Densities{IndDensity} num2str(IndDensity) '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(indeps-1) '.mat']);
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

ARs = [];
ESSs = [];
qs = 0.1:0.1:0.9;
EQs = [];
cols = {'k',':k','y','b','r','--r',':r','g'};
clf
hold on
for IndMethod = 1:6
    clf
    for IndLogOrNot = 0:0
        tmpESS = [];
        Epss = [];
        for indeps = 1:10
%             Epss = (1:14)*0.2;
%             Parameters.Epsil = Epss(indeps);
%             if and(IndMethod>=5,IndMethod<=8)
%                 Parameters.Epsil = Parameters.Epsil/10;
%             end
            try
                disp([Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat'])
                load([SavePath Densities{IndDensity} num2str(IndDensity)  '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(indeps-1) '.mat']);
                tmpESS(indeps) = min(Res.RelESS);
                Epss(indeps) = Res.Parameters.Epsil;
                ARs(IndMethod,indeps) = Res.AccRate;
                ESSs(IndMethod,indeps) =  min(Res.RelESS);
            end
        end
        try
            plot(Epss,tmpESS,cols{IndMethod})
        end
        
%         plot(Epss,tmpESS,cols{IndMethod})
        [b,imax] = max(tmpESS);
        for i = 1:length(qs)
            load([SavePath Densities{IndDensity} num2str(IndDensity)  '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(imax-1) '.mat']);
            EQs(IndMethod,i) = quantile(Res.Vals(2,:),qs(i));
        end
        title(max(tmpESS))
    end
%     pause()
end
hold off

clf
for IndMethod = 1:6
    plot(ARs(IndMethod,:),ESSs(IndMethod,:),cols{IndMethod})
    hold on
end
hold off
ylim([0 2])
xlim([0 1])


% ylim([0 80])

clf
for i = 1:2
    subplot(2,1,i)
    plot(Res.Vals(i,:))
end

