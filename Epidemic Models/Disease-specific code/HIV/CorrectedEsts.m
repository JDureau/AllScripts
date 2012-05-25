function [] = CorrectedEsts(indRes)


cd('/users/ecologie/dureau/src/AllScripts/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR'])


SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan/';

load([SavePath 'TreatedRessRightPriors.mat'])
ResNormal = Res;
load([SavePath 'TreatedRessRightPriors_Biased.mat'])
ResBiased = Res;
load([SavePath 'TreatedRessRightPriors_AddingObs.mat'])
ResAddingObs = Res;
load([SavePath 'TreatedRessRightPriors_Restricted.mat'])
ResRestricted = Res;


Ress = {};
Ress{1} = ResNormal;
Ress{2} = ResBiased;
Ress{3} = ResAddingObs;
Ress{4} = ResRestricted;




slopess = {};
for indType = 1:4
    Res = Ress{indType};
    slopes = [];
    for i = 1:length(Res.PrevObs)
        if not(indType==3)
            p = polyfit([472,528,584]/24,Res.PrevObs{i}(7,[2 3 5]),1);
        else
            p = polyfit([472,528,584]/24,Res.PrevObs{i}(7,[3 4 6]),1);
        end
        slopes(i) = p(1);
    end
    slopess{indType} = slopes;
    hist(slopes)
    xlim([-5 5])
%     pause()
end

slopes = slopess{indRes};
Res = Ress{indRes};

X = [slopes; Res.amplout; Res.inflptout; Res.propaftintout; Res.amplin; Res.inflptin; Res.propaftintin];
TRX = [X(1,:) ; X(2,:) ; log(X(3,:)./(600-X(3,:))) ; log(max(0.001,X(4,:))./(1-min(0.999,X(4,:)))); log(X(5,:)./(1-X(5,:))) ; log(X(6,:)./(600-X(6,:))); log(max(0.001,X(7,:))./(1-min(0.999,X(7,:))))];
mean(TRX')
Ns = 1:10;
BICS = [];
for i = 1:length(Ns)
    i
    TempDenss = {};
    TempBICS = [];
    for j = 1:10
        try
            TempDenss{end+1} =  gmdistribution.fit(TRX',Ns(i));
            TempBICS(end+1) = TempDenss{end}.BIC;
        end
    end
    try
        [b,ind] = min(TempBICS);
        Denss{i} = TempDenss{ind};
        BICS(i) = TempBICS(ind);
    end
end
BICS
plot(BICS)
[b,ind] = min(BICS)
Dens = Denss{i};


corramplout = [];
corrinflptout = [];
corrpropaftintout = [];
corramplout025 = [];
corrinflptout025 = [];
corrpropaftintout025 = [];
corramplout975 = [];
corrinflptout975 = [];
corrpropaftintout975 = [];
for i = 1:length(Res.amplin)
    i
    step = 0.01;
    xis = step:step:1-step;
    fis = [];
    for j = 1:length(xis)
        trf = 0;
        for indcomp = 1:Dens.NComponents
            tempx = [slopes(i); Res.amplout(i); log(Res.inflptout(i)/(600-Res.inflptout(i))); log(max(0.001,Res.propaftintout(i))/(1-min(0.999,Res.propaftintout(i)))); log(xis(j)/(600-xis(j)))];
            tempinds = [1 2 3 4 5];
            trf = trf + Dens.PComponents(indcomp)*mvnpdf(tempx',Dens.mu(indcomp,tempinds),Dens.Sigma(tempinds,tempinds,indcomp));
        end
        %fis(j) = trf*(600/(Pars.Par.Inflpt*(600-Pars.Par.Inflpt)))*(1/(xis(j)*(1-xis(j))));
        fis(j) = trf*(1/(xis(j)*(1-xis(j))));
    end
    cdf = cumsum((fis*mean(diff(xis)))/sum(fis*mean(diff(xis))));
    corramplout(i) = xis(find(cdf>0.5,1));
    corramplout025(i) = xis(find(cdf>0.025,1));
    corramplout975(i) = xis(find(cdf>0.975,1));
    
    step = 6;
    xis = step:step:600-step;
    fis = [];
    for j = 1:length(xis)
        trf = 0;
        for indcomp = 1:Dens.NComponents
            tempx = [slopes(i); Res.amplout(i); log(Res.inflptout(i)/(600-Res.inflptout(i))); log(max(0.001,Res.propaftintout(i))/(1-min(0.999,Res.propaftintout(i)))); log(xis(j)/(1-xis(j)))];
            tempinds = [1 2 3 4 6];
            trf = trf + Dens.PComponents(indcomp)*mvnpdf(tempx',Dens.mu(indcomp,tempinds),Dens.Sigma(tempinds,tempinds,indcomp));
        end
        %fis(j) = trf*(600/(Pars.Par.Inflpt*(600-Pars.Par.Inflpt)))*(600/(xis(j)*(600-xis(j))));
        fis(j) = trf*(600/(xis(j)*(600-xis(j))));
    end
    cdf = cumsum((fis*mean(diff(xis)))/sum(fis*mean(diff(xis))));
    corrinflptout(i) = xis(find(cdf>0.5,1));
    corrinflout025(i) = xis(find(cdf>0.025,1));
    corrinflout975(i) = xis(find(cdf>0.975,1));
    
    step = 0.01;
    xis = step:step:1-step;
    fis = [];
    for j = 1:length(xis)
        trf = 0;
        for indcomp = 1:Dens.NComponents
            tempx = [slopes(i); Res.amplout(i); log(Res.inflptout(i)/(600-Res.inflptout(i))); log(max(0.001,Res.propaftintout(i))/(1-min(0.999,Res.propaftintout(i)))); log(xis(j)/(1-xis(j)))];
            tempinds = [1 2 3 4 7];
            trf = trf + Dens.PComponents(indcomp)*mvnpdf(tempx',Dens.mu(indcomp,tempinds),Dens.Sigma(tempinds,tempinds,indcomp));
        end
        %fis(j) = trf*(600/(Pars.Par.Inflpt*(600-Pars.Par.Inflpt)))*(1/(xis(j)*(1-xis(j))));
        fis(j) = trf*(1/(xis(j)*(1-xis(j))));
    end
    cdf = cumsum((fis*mean(diff(xis)))/sum(fis*mean(diff(xis))));
    corrpropaftintout(i) = xis(find(cdf>0.5,1));
    corrpropaftintout025(i) = xis(find(cdf>0.025,1));
    corrpropaftintout975(i) = xis(find(cdf>0.975,1));
end

Res.corrinflptout = corrinflptout;
Res.corramplout = corramplout;
Res.corrinflptout = corrinflptout;
Res.corrpropaftintout025 =corrpropaftintout025;
Res.corramplout025 = corramplout025;
Res.corrinflptout025 = corrinflptout025;
Res.corrpropaftintout975 =corrpropaftintout975;
Res.corramplout975 = corramplout975;
Res.corrinflptout975 = corrinflptout975;

if indRes == 1
    save('TreatedRessRightPriorsWithCorrs.mat','Res')
elseif indRes == 2
    save('TreatedRessRightPriorsWithCorrs_Biased.mat','Res')
elseif indRes == 3
    save('TreatedRessRightPriorsWithCorrs_AddingObs.mat','Res')
elseif indRes == 4
    save('TreatedRessRightPriorsWithCorrs_Restricted.mat','Res')
end

% load(['TreatedRessRightPriors.mat'])
% ResNormal = Res;
% load(['TreatedRessRightPriors_Biased.mat'])
% ResBiased = Res;
% load(['TreatedRessRightPriors_AddingObs.mat'])
% ResAddingObs = Res;
% load(['TreatedRessRightPriors_Restricted.mat'])
% ResRestricted = Res;

