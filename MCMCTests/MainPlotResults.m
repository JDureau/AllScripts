% Results analysis

SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/MCMCstudy/';

%% STEP 1: check if they all get to the same posteriors

Methods =  {'MALA','LocalMALA','RMALA','GMMRand','GMMLang'};
Densities = {'GMM2','Banana'};

IndDensity = 1;
dim = 2;

cols = {'k','','','b','r'};

clf
for IndMethod = [1 43 4 5]
    for IndLogOrNot = 0:1
        for indeps = 1:14
            Epss = (1:14)*0.2;
            Parameters.Epsil = Epss(indeps);
            try
                disp([Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat'])
                load([SavePath Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat']);
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
end

%% STEP 2: check ESS plots

cols = {'k','','','b','r'};
clf
hold on
for IndMethod = [1 4 5]
    for IndLogOrNot = 0:1
        tmpESS = [];
        for indeps = 1:14
            Epss = (1:14)*0.2;
            Parameters.Epsil = Epss(indeps);
            try
                disp([Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat'])
                load([SavePath Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat']);
                tmpESS(indeps) = min(Res.RelESS);
            end
        end
        try
            plot(Epss,tmpESS,cols{IndMethod})
        end
    end
end
hold off
ylim([0 2])



