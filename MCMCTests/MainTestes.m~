%% All tests on MCMC

pwd
cd('/Users/dureaujoseph/Dropbox/AllScriptsGit/')
addpath([pwd '/MCMCTests/'])
addpath([pwd '/General Tools/'])
addpath([pwd '/Toolboxes/'])
% SavePath = '/Users/dureaujoseph/Documents/Taf/These/Results/GMM_MCMC';

% cd('/Users/dureaujoseph/Dropbox/Taf/AllScripts/')
% addpath([pwd '\MCMCTests\'])
% addpath([pwd 'General Tools\'])
% addpath([pwd '\Toolboxes\'])
% SavePath = 'S:\ResultsMCMC';




Methods =  {'MALA','LocalMALA','RMALA','GMMRand','GMMLang'};
Densities = {'GMM2','Banana'};


Epss = (1:14)*0.2;
Parameters.Epsil = Epss(ind);


switch IndDensity
    case 1
        Parameters.f = @fGMM;
        mu1 = zeros(dim,1);
        mu2 = zeros(dim,1);
        mu2(1) = 4;
        sigma  = eye(dim);
        X = [mvnrnd(mu1,sigma,20000)' mvnrnd(mu2,sigma,10000)']';
        Parameters.RealDens = gmdistribution.fit(X,2);
        Parameters.Dim = dim;
        Parameters = FindFisherInfMat(zeros(1,dim),Parameters) ;
        Parameters.OptDens = Parameters.RealDens;
        Parameters.Dens = gmdistribution.fit(X,2);
        Parameters.G = @GGMM;
        Parameters.GDerivs = @GDerivsGMM;
    case 2
        Parameters.f = @fBanana;
        B = 0.1;
        X = mvnrnd(zeros(dim,1),eye(dim),10000);
        X(:,1) = 10*X(:,1);
        X(:,2) = X(:,2)-B*X(:,1).^2+100*B;
        Parameters.RealDens = gmdistribution.fit(X,15);
        Parameters.OptDens = Parameters.RealDens;
        Parameters.B = B;
        Parameters.Dim = dim;
        Parameters = FindFisherInfMat(zeros(1,dim),Parameters) ;
        test = 0;
        while not(test)
            try
                Parameters.Dens = gmdistribution.fit(X,14);
                test = 1;
            end
        end
end
       
switch IndMethod
    case 1
        Parameters.LogRatioFun = @LogRatioMALA;
        Parameters.SampleFun = @SampleMALA;
        Parameters.ScalingCov = -(Parameters.Hess^-1);
    case 2
        Parameters.LogRatioFun = @LogRatioLocalMALA;
        Parameters.SampleFun = @SampleLocalMALA;
    case 3
        Parameters.LogRatioFun = @LogRatioRMALA;
        Parameters.SampleFun = @SampleRMALA;
    case 4
        Parameters.LogRatioFun = @LogRatioGMMRand;
        Parameters.SampleFun = @SampleGMMRand;
    case 5
        Parameters.LogRatioFun = @LogRatioGMMLang;
        Parameters.SampleFun = @SampleGMMLang;
end

switch IndLogOrNot 
    case 0
        Parameters.Mode = 'L';
    case 1
        Parameters.Mode = 'Log';
end


Samples = {};
AccRates = [];
RelESSs = [];
Res = RunMCMC(Parameters.ArgMax',Parameters,5000);
Res.Samples = Res.Vals;
Res.Eps = Parameters.Epsil;

SavePath = '/users/ecologie/dureau/src/AllData/GMM/';
save([SavePath Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat'],'Res');

      





















% load and plot
Methods =  {'MALA','GMMRand','GMMLang','GMMLang2'};
Densities = {'GMM2','Banana'};

IndDensity = 1;
IndsMethods = 1:4;
IndsLogOrNot = 0:1;
Dims = 2;

Epss = (1:13)*0.2;
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/GMM_MCMC/';
AccRates = [];
RelEsss = [];
Samples = [];
cpt = 0;
legends = {};
for dim = Dims
    for IndMethod = IndsMethods
        for IndLogOrNot = IndsLogOrNot;
            cpt = cpt+1;
            for i = 1:length(Epss)
                load([SavePath Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Epss(i)) '.mat'
]);
                AccRates(cpt,i) = Res.AccRate;
                RelEsss(cpt,i) = min(Res.RelESS);
            end
            [b,ind] = max(RelEsss(cpt,:));
            load([SavePath Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Epss(ind)) '.mat']);
            Samples(cpt,:,:) = Res.Samples;
            legends{cpt} = [Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot)];
        end
    end
end
clf
% for i = 1:dim
%     subplot(2,dim,i)
%     for j = 1:cpt
%         hold on
%         [fi,xi] = ksdensity(squeeze(Samples(j,i,:)));
%         plot(xi,fi);
%     end
% end
% subplot(2,dim,dim+1:2*dim)
plot(AccRates',RelEsss')
legend(legends)


ESSs = {};
AccRates = {};
RelEsss = {};
for i = 1:length(Dims)
    ESSs{i} = {};
    for j = 1:length(Names)
        load([SavePath Names{j} num2str(Dims(i)) '.mat'])
        ESSs{i}{j} = Res.RelESSs;
        AccRates{i}{j} = Res.AccRates;
    end
end

clf
ms = [];
inds = [];
xms = [];
for i = 1:length(Dims)
    subplot(3,4,4+i)
    xs = AccRates{i}{1};
    ys = min(ESSs{i}{1}');
    [ms(1),inds(1)] = max(ys);
    xms(1) = xs(inds(1));
    plot(xs,ys,'g','LineWidth',2)
    hold on
    xs = AccRates{i}{2};
    ys = min(ESSs{i}{2}');
    [ms(2),inds(2)] = max(ys);
    xms(2) = xs(inds(2));
    plot(xs,ys,'b','LineWidth',2)
    xs = AccRates{i}{3};
    ys = min(ESSs{i}{3}');
    [ms(3),inds(3)] = max(ys);
    xms(3) = xs(inds(3));
    plot(xs,ys,'k','LineWidth',2)
    for j = 1:3
        plot(xms(j)*ones(100,1),ms(j)/100:ms(j)/100:ms(j),'--k')
    end
    hold off
    %     legend('MALA','GMMRand','GMMLang')
    title(['Dim = ' num2str(Dims(i))])
    ylim([0 30])
    set(gca,'XTick',0:0.2:1)
    xlabel('Acceptance rate')
    ylabel('min(RelESS)')
    
    subplot(3,4,8+i)
    xs = Epss;
    ys = min(ESSs{i}{1}');
    [ms(1),inds(1)] = max(ys);
    plot(xs,ys,'g','LineWidth',2)
    hold on
    ys = min(ESSs{i}{2}');
    [ms(2),inds(2)] = max(ys);
    plot(xs,ys,'b','LineWidth',2)
    ys = min(ESSs{i}{3}');
    [ms(3),inds(3)] = max(ys);
    plot(xs,ys,'k','LineWidth',2)
    for j = 1:3
        plot(xs(inds(j))*ones(100,1),ms(j)/100:ms(j)/100:ms(j),'--k')
    end
    hold off
    %     legend('MALA','GMMRand','GMMLang')
    title(['Dim = ' num2str(Dims(i))])
    xlabel('\lambda')
    ylim([0 30])
    ylabel('min(RelESS)')
end   
legend('Simple MALA','RW (or PGMM(1))','AGMM (1)')
subplot(3,4,2:3)
scattercloud(Res.Samples{inds(3)}(1,:),Res.Samples{inds(3)}(2,:),25)
title('2D projection of the targeted density')


%% GMM 1

% Names = {'/Gaussian_MALATests_','/Gaussian_GMMRandTests_','/Gaussian_GMMLangTests_'};
Names = {'\Gaussian_MALATests_','\Gaussian_GMMRandTests_','\Gaussian_GMMLangTests_'};
Epss = (1:9)*0.2;
Dims = [2 5 10 20];

for IndName = 1:length(Names)
    for IndDims = 1:length(Dims)
        dim = Dims(IndDims);
        Parameters.f = @fGMM;
        mu = zeros(dim,1);
        sigma  = eye(dim);
        X = mvnrnd(mu,sigma,20000);
        Parameters.RealDens = gmdistribution.fit(X,1);
        Parameters.Dim = dim;
        Parameters = FindFisherInfMat(zeros(1,dim),Parameters) ;
        Parameters.Mode = 'Log';
        
        switch IndName
            case 1
                Parameters.LogRatioFun = @LogRatioMALA;
                Parameters.SampleFun = @SampleMALA;
                Parameters.ScalingCov = -(Parameters.Hess^(-1));
            case 2
                Parameters.LogRatioFun = @LogRatioGMMrand;
                Parameters.SampleFun = @SampleGMMrand;
            case 3
                Parameters.LogRatioFun = @LogRatioGMMLang;
                Parameters.SampleFun = @SampleGMMLang;
        end
        
        Samples = {};
        AccRates = [];
        RelESSs = [];
        Parameters.OptDens = Parameters.RealDens;
        Parameters.Dens = Parameters.RealDens;
        for j = 1:length(Epss)
            disp([Names{IndName} ' ' num2str(Dims(IndDims)) ' ' num2str(Epss(j))])
            Parameters.Epsil = Epss(j);
            if or(dim==2, dim==5)
                Res = RunMCMC(Parameters.ArgMax',Parameters,400);
            else
                Res = RunMCMC(Parameters.ArgMax',Parameters,200);
            end
            AccRates(j) = Res.AccRate;
            RelESSs(j,:) = Res.RelESS;
            Samples{j} =  Res.Vals;
        end
        Res.AccRates = AccRates;
        Res.RelESSs = RelESSs;
        Res.Samples = Samples;
        Res.Epss = Epss;
        save([SavePath Names{IndName} num2str(dim) '.mat'],'Res');
    end
end
        

% plot

ESSs = {};
AccRates = {};
for i = 1:length(Dims)
    ESSs{i} = {};
    for j = 1:length(Names)
        load([SavePath Names{j} num2str(Dims(i)) '.mat'])
        ESSs{i}{j} = Res.RelESSs;
        AccRates{i}{j} = Res.AccRates;
    end
end

clf
ms = [];
inds = [];
xms = [];
for i = 1:length(Dims)
    subplot(3,4,4+i)
    xs = AccRates{i}{1};
    ys = min(ESSs{i}{1}');
    [ms(1),inds(1)] = max(ys);
    xms(1) = xs(inds(1));
    plot(xs,ys,'g','LineWidth',2)
    hold on
    xs = AccRates{i}{2};
    ys = min(ESSs{i}{2}');
    [ms(2),inds(2)] = max(ys);
    xms(2) = xs(inds(2));
    plot(xs,ys,'b','LineWidth',2)
    xs = AccRates{i}{3};
    ys = min(ESSs{i}{3}');
    [ms(3),inds(3)] = max(ys);
    xms(3) = xs(inds(3));
    plot(xs,ys,'k','LineWidth',2)
    for j = 1:3
        plot(xms(j)*ones(100,1),ms(j)/100:ms(j)/100:ms(j),'--k')
    end
    hold off
    %     legend('MALA','GMMRand','GMMLang')
    title(['Dim = ' num2str(Dims(i))])
    ylim([0 30])
    set(gca,'XTick',0:0.2:1)
    xlabel('Acceptance rate')
    ylabel('min(RelESS)')
    
    subplot(3,4,8+i)
    xs = Epss;
    ys = min(ESSs{i}{1}');
    [ms(1),inds(1)] = max(ys);
    plot(xs,ys,'g','LineWidth',2)
    hold on
    ys = min(ESSs{i}{2}');
    [ms(2),inds(2)] = max(ys);
    plot(xs,ys,'b','LineWidth',2)
    ys = min(ESSs{i}{3}');
    [ms(3),inds(3)] = max(ys);
    plot(xs,ys,'k','LineWidth',2)
    for j = 1:3
        plot(xs(inds(j))*ones(100,1),ms(j)/100:ms(j)/100:ms(j),'--k')
    end
    hold off
    %     legend('MALA','GMMRand','GMMLang')
    title(['Dim = ' num2str(Dims(i))])
    xlabel('\lambda')
    ylim([0 30])
    ylabel('min(RelESS)')
end   
legend('Simple MALA','RW (or PGMM(1))','AGMM (1)')
subplot(3,4,2:3)
scattercloud(Res.Samples{inds(3)}(1,:),Res.Samples{inds(3)}(2,:),25)
title('2D projection of the targeted density')

        
%% GMM 2

         

Names = {'\2Gaussian_MALATests_','\2Gaussian_LogGMMRandTests_','\2Gaussian_LogGMMLangTests_','\2Gaussian_LogGMMAndrieuRandTests_'};%,'\2Gaussian_LGMMLangTests_'};
Epss = (1:9)*0.2;
Dims = [2 5];

for IndName = 1:length(Names)
    for IndDims = 1:length(Dims)
        dim = Dims(IndDims);
        Parameters.f = @fGMM;
           
        mu1 = zeros(dim,1);
        mu2 = zeros(dim,1);
        mu2(1) = 4;
        sigma  = eye(dim);
        X = [mvnrnd(mu1,sigma,20000)' mvnrnd(mu2,sigma,10000)']';
        Parameters.RealDens = gmdistribution.fit(X,2);
        Parameters.Dim = dim;
        Parameters = FindFisherInfMat(zeros(1,dim),Parameters) ;
        
        
        switch IndName
            case 1
                Parameters.LogRatioFun = @LogRatioMALA;
                Parameters.SampleFun = @SampleMALA;
                Parameters.ScalingCov = -(Parameters.Hess^-1);
            case 2
                Parameters.Mode = 'Log';
                Parameters.LogRatioFun = @LogRatioGMMRand;
                Parameters.SampleFun = @SampleGMMRand;
            case 3
                Parameters.Mode = 'Log';
                Parameters.LogRatioFun = @LogRatioGMMLang;
                Parameters.SampleFun = @SampleGMMLang;
            case 4
                Parameters.Mode = 'Log';
                Parameters.LogRatioFun = @LogRatioGMMAndrieuRand;
                Parameters.SampleFun = @SampleGMMAndrieuRand;
            case 5
                Parameters.Mode = 'L';
                Parameters.LogRatioFun = @LogRatioGMMLang;
                Parameters.SampleFun = @SampleGMMLang;
        end
        
        Samples = {};
        AccRates = [];
        RelESSs = [];
        Parameters.OptDens = Parameters.RealDens;
        Parameters.Dens = gmdistribution.fit(X,5);
        for j = 1:length(Epss)
            disp([Names{IndName} ' ' num2str(Dims(IndDims)) ' ' num2str(Epss(j))])
            Parameters.Epsil = Epss(j);
            if or(dim==2,dim==5)
                Res = RunMCMC(Parameters.ArgMax',Parameters,40000);
            else
                Res = RunMCMC(Parameters.ArgMax',Parameters,20000);
            end
            AccRates(j) = Res.AccRate;
            RelESSs(j,:) = Res.RelESS;
            Samples{j} =  Res.Vals;
        end
        Res.AccRates = AccRates;
        Res.RelESSs = RelESSs;
        Res.Samples = Samples;
        Res.Epss = Epss;
        save([SavePath Names{IndName} num2str(dim) '_3GMMtest.mat'],'Res');
    end
end
        

% plot

ESSs = {};
AccRates = {};
Samples = {};
for i = 1:length(Dims)
    ESSs{i} = {};
    for j = 1:4%1:length(Names)
        load([SavePath Names{j} num2str(Dims(i)) '.mat'])
        [b,ind] = max(min(Res.RelESSs'));
        clf
        plot(Res.Samples{ind}(1,:),Res.Samples{ind}(2,:),'.')
        disp([num2str(Dims(i)) Names{j}])
%         pause()
        ESSs{i}{j} = Res.RelESSs;
        AccRates{i}{j} = Res.AccRates;
        Samples{i}{j} = Res.Vals;
    end
end

clf
ms = [];
inds = [];
xms = [];
for i = 1:length(Dims)
    subplot(3,4,4+i)
    xs = AccRates{i}{1};
    ys = min(ESSs{i}{1}');
    [ms(1),inds(1)] = max(ys);
    xms(1) = xs(inds(1));
    plot(xs,ys,'g','LineWidth',2)
    hold on
    xs = AccRates{i}{2};
    ys = min(ESSs{i}{2}');
    [ms(2),inds(2)] = max(ys);
    xms(2) = xs(inds(2));
    plot(xs,ys,'b','LineWidth',2)
    xs = AccRates{i}{3};
    ys = min(ESSs{i}{3}');
    [ms(3),inds(3)] = max(ys);
    xms(3) = xs(inds(3));
    plot(xs,ys,'k','LineWidth',2)
    for j = 1:3
        plot(xms(j)*ones(100,1),ms(j)/100:ms(j)/100:ms(j),'--k')
    end
    hold off
    %     legend('MALA','GMMRand','GMMLang')
    title(['Dim = ' num2str(Dims(i))])
    ylim([0 5.5])
    set(gca,'XTick',0:0.2:1)
    xlabel('Acceptance rate')
    ylabel('min(RelESS)')
    
    subplot(3,4,8+i)
    xs = Epss;
    ys = min(ESSs{i}{1}');
    [ms(1),inds(1)] = max(ys);
    plot(xs,ys,'g','LineWidth',2)
    hold on
    ys = min(ESSs{i}{2}');
    [ms(2),inds(2)] = max(ys);
    plot(xs,ys,'b','LineWidth',2)
    ys = min(ESSs{i}{3}');
    [ms(3),inds(3)] = max(ys);
    plot(xs,ys,'k','LineWidth',2)
    for j = 1:3
        plot(xs(inds(j))*ones(100,1),ms(j)/100:ms(j)/100:ms(j),'--k')
    end
    hold off
    %     legend('MALA','GMMRand','GMMLang')
    title(['Dim = ' num2str(Dims(i))])
    ylim([0 5.5])
    xlabel('\lambda')
    ylabel('min(RelESS)')
end   
legend('Simple MALA','PGMM (2)','AGMM (2)')
subplot(3,4,2:3)
scattercloud(Res.Samples{inds(3)}(1,:),Res.Samples{inds(3)}(2,:),25)
title('2D projection of the targeted density')




clf
ms = [];
inds = [];
xms = [];
for i = 1:length(Dims)
    subplot(length(Dims),1,i)
    dim = i;
    mu1 = zeros(dim,1);
    mu2 = zeros(dim,1);
    mu2(1) = 4;
    sigma  = eye(dim);
    X = [mvnrnd(mu1,sigma,20000)' mvnrnd(mu2,sigma,10000)']';
    [f,xi] = ksdensity(X(:,1));
    plot(xi,f,'r')
    hold on
    [f,xi] = ksdensity(Samples{i}{1}(1,:));
    plot(xi,f,'g')
    [f,xi] = ksdensity(Samples{i}{2}(1,:));
    plot(xi,f,'b')
    [f,xi] = ksdensity(Samples{i}{4}(1,:));
    plot(xi,f,'c')
    title(dim)
    hold off
end   
legend('aim','Simple MALA','PGMM (2)','AGMM (2)')

        
%% GMM 5

         

Names = {'\5Gaussian_MALATests_','\5Gaussian_LogGMMRandTests_','\5Gaussian_LogGMMLangTests_','\5Gaussian_LogGMMAndrieuRandTests_'};%,'\2Gaussian_LGMMLangTests_'};
Epss = (1:9)*0.2;
Dims = [2 5];

for IndName = 1:length(Names)
    for IndDims = 1:length(Dims)
        dim = Dims(IndDims);
        Parameters.f = @fGMM;
           
        mu1 = zeros(dim,1);
        mu2 = mu1;
        mu3 = mu1;
        mu4 = mu1;
        mu5 = mu1;
        
        mu2(1) = 1;
        mu3(1) = 2;
        mu4(1) = 3;
        mu5(1) = 4;
        sigma  = eye(dim);
        X = [mvnrnd(mu1,sigma,10000)' mvnrnd(mu2,sigma,10000)' mvnrnd(mu3,sigma,10000)' mvnrnd(mu4,sigma,10000)' mvnrnd(mu5,sigma,10000)']';
        Parameters.RealDens = gmdistribution.fit(X,5);
        Parameters.Dim = dim;
        Parameters = FindFisherInfMat(zeros(1,dim),Parameters) ;
        
        
        switch IndName
            case 1
                Parameters.LogRatioFun = @LogRatioMALA;
                Parameters.SampleFun = @SampleMALA;
                Parameters.ScalingCov = -(Parameters.Hess^-1);
            case 2
                Parameters.Mode = 'Log';
                Parameters.LogRatioFun = @LogRatioGMMRand;
                Parameters.SampleFun = @SampleGMMRand;
            case 3
                Parameters.Mode = 'Log';
                Parameters.LogRatioFun = @LogRatioGMMLang;
                Parameters.SampleFun = @SampleGMMLang;
            case 4
                Parameters.Mode = 'Log';
                Parameters.LogRatioFun = @LogRatioGMMAndrieuRand;
                Parameters.SampleFun = @SampleGMMAndrieuRand;
            case 5
                Parameters.Mode = 'L';
                Parameters.LogRatioFun = @LogRatioGMMLang;
                Parameters.SampleFun = @SampleGMMLang;
        end
        
        Samples = {};
        AccRates = [];
        RelESSs = [];
        Parameters.OptDens = Parameters.RealDens;
        Parameters.Dens = gmdistribution.fit(X,5);
        for j = 1:length(Epss)
            disp([Names{IndName} ' ' num2str(Dims(IndDims)) ' ' num2str(Epss(j))])
            Parameters.Epsil = Epss(j);
            if or(dim==2,dim==5)
                Res = RunMCMC(Parameters.ArgMax',Parameters,40000);
            else
                Res = RunMCMC(Parameters.ArgMax',Parameters,20000);
            end
            AccRates(j) = Res.AccRate;
            RelESSs(j,:) = Res.RelESS;
            Samples{j} =  Res.Vals;
        end
        Res.AccRates = AccRates;
        Res.RelESSs = RelESSs;
        Res.Samples = Samples;
        Res.Epss = Epss;
        save([SavePath Names{IndName} num2str(dim) '_5GMMtest.mat'],'Res');
    end
end
        

% plot

ESSs = {};
AccRates = {};
Samples = {};
for i = 1:length(Dims)
    ESSs{i} = {};
    for j = 1:4%1:length(Names)
        load([SavePath Names{j} num2str(Dims(i)) '_5GMMtest.mat'])
        [b,ind] = max(min(Res.RelESSs'));
        clf
        plot(Res.Samples{ind}(1,:),Res.Samples{ind}(2,:),'.')
        disp([num2str(Dims(i)) Names{j}])
%         pause()
        ESSs{i}{j} = Res.RelESSs;
        AccRates{i}{j} = Res.AccRates;
        Samples{i}{j} = Res.Vals;
    end
end

clf
ms = [];
inds = [];
xms = [];
for i = 1:length(Dims)
    subplot(3,4,4+i)
    xs = AccRates{i}{1};
    ys = min(ESSs{i}{1}');
    [ms(1),inds(1)] = max(ys);
    xms(1) = xs(inds(1));
    plot(xs,ys,'g','LineWidth',2)
    hold on
    xs = AccRates{i}{2};
    ys = min(ESSs{i}{2}');
    [ms(2),inds(2)] = max(ys);
    xms(2) = xs(inds(2));
    plot(xs,ys,'b','LineWidth',2)
    xs = AccRates{i}{3};
    ys = min(ESSs{i}{3}');
    [ms(3),inds(3)] = max(ys);
    xms(3) = xs(inds(3));
    plot(xs,ys,'k','LineWidth',2)
    for j = 1:3
        plot(xms(j)*ones(100,1),ms(j)/100:ms(j)/100:ms(j),'--k')
    end
    hold off
    %     legend('MALA','GMMRand','GMMLang')
    title(['Dim = ' num2str(Dims(i))])
    ylim([0 5.5])
    set(gca,'XTick',0:0.2:1)
    xlabel('Acceptance rate')
    ylabel('min(RelESS)')
    
    subplot(3,4,8+i)
    xs = Epss;
    ys = min(ESSs{i}{1}');
    [ms(1),inds(1)] = max(ys);
    plot(xs,ys,'g','LineWidth',2)
    hold on
    ys = min(ESSs{i}{2}');
    [ms(2),inds(2)] = max(ys);
    plot(xs,ys,'b','LineWidth',2)
    ys = min(ESSs{i}{3}');
    [ms(3),inds(3)] = max(ys);
    plot(xs,ys,'k','LineWidth',2)
    for j = 1:3
        plot(xs(inds(j))*ones(100,1),ms(j)/100:ms(j)/100:ms(j),'--k')
    end
    hold off
    %     legend('MALA','GMMRand','GMMLang')
    title(['Dim = ' num2str(Dims(i))])
    ylim([0 5.5])
    xlabel('\lambda')
    ylabel('min(RelESS)')
end   
legend('Simple MALA','PGMM (2)','AGMM (2)')
subplot(3,4,2:3)
scattercloud(Res.Samples{inds(3)}(1,:),Res.Samples{inds(3)}(2,:),25)
title('2D projection of the targeted density')




clf
ms = [];
inds = [];
xms = [];
for i = 1:length(Dims)
    subplot(length(Dims),1,i)
    dim = i;
    mu1 = zeros(dim,1);
    mu2 = zeros(dim,1);
    mu2(1) = 4;
    sigma  = eye(dim);
    X = [mvnrnd(mu1,sigma,20000)' mvnrnd(mu2,sigma,10000)']';
    [f,xi] = ksdensity(X(:,1));
    plot(xi,f,'r')
    hold on
    [f,xi] = ksdensity(Samples{i}{1}(1,:));
    plot(xi,f,'g')
    [f,xi] = ksdensity(Samples{i}{2}(1,:));
    plot(xi,f,'b')
    [f,xi] = ksdensity(Samples{i}{4}(1,:));
    plot(xi,f,'c')
    title(dim)
    hold off
end   
legend('aim','Simple MALA','PGMM (2)','AGMM (2)')

        
      
%% Banana

         
% Names = {'/Banana_MALATests_','/Banana_3GMMRandTests_','/Banana_3GMMLangTests_','/Banana_7GMMRandTests_','/Banana_7GMMLangTests_'};
Names = {'\Banana_MALATests_','\Banana_14LogGMMRandTests_','\Banana_14LogGMMLangTests_','\Banana_20LogGMMRandTests_','\Banana_20LogGMMLangTests_'};%,'\Banana_3LGMMRandTests_','\Banana_3LGMMLangTests_','\Banana_7LGMMRandTests_','\Banana_7LGMMLangTests_'};
Epss = (1:9)*0.2;
Dims = [2 5 10 20];
    
for IndName = 1:length(Names)
    for IndDims = 1:length(Dims)
        dim = Dims(IndDims);
        Parameters.f = @fBanana;
           
        B = 0.1;
        X = mvnrnd(zeros(dim,1),eye(dim),10000);
        X(:,1) = 10*X(:,1);
        X(:,2) = X(:,2)-B*X(:,1).^2+100*B;
        Parameters.RealDens = gmdistribution.fit(X,15);
        Parameters.OptDens = Parameters.RealDens;
        Parameters.B = B;
        
        Parameters.Dim = dim;
        Parameters = FindFisherInfMat(zeros(1,dim),Parameters) ;
        
        test = 0;
        while not(test)
            try
                switch IndName
                    case 1
                        Parameters.LogRatioFun = @LogRatioMALA;
                        Parameters.SampleFun = @SampleMALA;
                        Parameters.ScalingCov = -(Parameters.Hess^-1);
                    case 2
                        Parameters.Mode = 'Log';
                        Parameters.LogRatioFun = @LogRatioGMMRand;
                        Parameters.SampleFun = @SampleGMMRand;
                        Parameters.Dens = gmdistribution.fit(X,14);
                    case 3
                        Parameters.Mode = 'Log';
                        Parameters.LogRatioFun = @LogRatioGMMLang;
                        Parameters.SampleFun = @SampleGMMLang;
                        Parameters.Dens = gmdistribution.fit(X,14);
                    case 4
                        Parameters.Mode = 'Log';
                        Parameters.LogRatioFun = @LogRatioGMMRand;
                        Parameters.SampleFun = @SampleGMMRand;
                        Parameters.Dens = gmdistribution.fit(X,20);
                    case 5
                        Parameters.Mode = 'Log';
                        Parameters.LogRatioFun = @LogRatioGMMLang;
                        Parameters.SampleFun = @SampleGMMLang;
                        Parameters.Dens = gmdistribution.fit(X,20);
                        hold on
                            ezcontour(@(x,y)pdf(Parameters.Dens,[x y]),[-20 20],[-20 20]);
                        hold off
%                     case 6
%                         Parameters.Mode = 'L';
%                         Parameters.LogRatioFun = @LogRatioGMMRand;
%                         Parameters.SampleFun = @SampleGMMRand;
%                         Parameters.Dens = gmdistribution.fit(X,3);
%                     case 7
%                         Parameters.Mode = 'L';
%                         Parameters.LogRatioFun = @LogRatioGMMLang;
%                         Parameters.SampleFun = @SampleGMMLang;
%                         Parameters.Dens = gmdistribution.fit(X,3);
%                     case 8
%                         Parameters.Mode = 'L';
%                         Parameters.LogRatioFun = @LogRatioGMMRand;
%                         Parameters.SampleFun = @SampleGMMRand;
%                         Parameters.Dens = gmdistribution.fit(X,7);
%                     case 9
%                         Parameters.Mode = 'L';
%                         Parameters.LogRatioFun = @LogRatioGMMLang;
%                         Parameters.SampleFun = @SampleGMMLang;
%                         Parameters.Dens = gmdistribution.fit(X,7);
                end
            test = 1;    
            end
        end
        Samples = {};
        AccRates = [];
        RelESSs = [];
        Parameters.OptDens = Parameters.RealDens;
        Parameters.Dens = Parameters.RealDens;
        for j = 1:length(Epss)
            disp([Names{IndName} ' ' num2str(Dims(IndDims)) ' ' num2str(Epss(j))])
            Parameters.Epsil = Epss(j);
            if or(dim==2,dim==5)
                Res = RunMCMC(Parameters.ArgMax',Parameters,80000);
            else
                Res = RunMCMC(Parameters.ArgMax',Parameters,40000);
            end
            AccRates(j) = Res.AccRate;
            RelESSs(j,:) = Res.RelESS;
            Samples{j} =  Res.Vals;
        end
        Res.AccRates = AccRates;
        Res.RelESSs = RelESSs;
        Res.Samples = Samples;
        Res.Epss = Epss;
        save([SavePath Names{IndName} num2str(dim) '.mat'],'Res');
    end
end
        

% plot

ESSs = {};
AccRates = {};
for i = 1:length(Dims)
    ESSs{i} = {};
    for j = 1:5%length(Names)
        load([SavePath Names{j} num2str(Dims(i)) '.mat'])
        ESSs{i}{j} = Res.RelESSs;
        AccRates{i}{j} = Res.AccRates;
    end
end

clf
ms = [];
inds = [];
xms = [];
for i = 1:length(Dims)
    subplot(3,4,4+i)
    xs = AccRates{i}{1};
    ys = min(ESSs{i}{1}');
    [ms(1),inds(1)] = max(ys);
    xms(1) = xs(inds(1));
    plot(xs,ys,'g','LineWidth',2)
    hold on
    xs = AccRates{i}{2};
    ys = min(ESSs{i}{2}');
    [ms(2),inds(2)] = max(ys);
    xms(2) = xs(inds(2));
    plot(xs,ys,'b','LineWidth',2)
    xs = AccRates{i}{3};
    ys = min(ESSs{i}{3}');
    [ms(3),inds(3)] = max(ys);
    xms(3) = xs(inds(3));
    plot(xs,ys,'k','LineWidth',2)
    xs = AccRates{i}{4};
    ys = min(ESSs{i}{4}');
    [ms(4),inds(4)] = max(ys);
    xms(4) = xs(inds(4));
    plot(xs,ys,'b--','LineWidth',2)
    xs = AccRates{i}{5};
    ys = min(ESSs{i}{5}');
    [ms(5),inds(5)] = max(ys);
    xms(5) = xs(inds(5));
    plot(xs,ys,'--k','LineWidth',2)
    for j = 1:3
        plot(xms(j)*ones(100,1),ms(j)/100:ms(j)/100:ms(j),'--k')
    end
    hold off
    %     legend('MALA','GMMRand','GMMLang')
    title(['Dim = ' num2str(Dims(i))])
%     ylim([0 5.5])
    set(gca,'XTick',0:0.2:1)
    xlabel('Acceptance rate')
    ylabel('min(ESS)')
    
    subplot(3,4,8+i)
    xs = Epss;
    ys = min(ESSs{i}{1}');
    [ms(1),inds(1)] = max(ys);
    plot(xs,ys,'g','LineWidth',2)
    hold on
    ys = min(ESSs{i}{2}');
    [ms(2),inds(2)] = max(ys);
    plot(xs,ys,'b','LineWidth',2)
    ys = min(ESSs{i}{3}');
    [ms(3),inds(3)] = max(ys);
    plot(xs,ys,'k','LineWidth',2)
    for j = 1:3
        plot(xs(inds(j))*ones(100,1),ms(j)/100:ms(j)/100:ms(j),'--k')
    end
    hold off
    %     legend('MALA','GMMRand','GMMLang')
    title(['Dim = ' num2str(Dims(i))])
%     ylim([0 5.5])
    xlabel('\lambda')
    ylabel('min(ESS)')
end   
legend('Simple MALA','PGMM (2)','AGMM (2)')
subplot(3,4,2:3)
scattercloud(Res.Samples{inds(3)}(1,:),Res.Samples{inds(3)}(2,:),25)
title('2D projection of the targeted density')
        
        
        
      