%% AGMM : learning the densities

pwd
cd('/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/')
addpath([pwd '/MCMCTests/'])
addpath([pwd '/General Tools/'])
addpath([pwd '/Toolboxes/mcmcdiag/'])
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Results/GMM_MCMC';

cd('H:\My Documents\PhD Work\Matlab Scripts From Macbook')
addpath([pwd '\MCMCTests\'])
addpath([pwd '\General Tools\'])
addpath([pwd '\Toolboxes\'])
SavePath = 'S:\ResultsMCMC';

        
%% GMM 2

         

Names = {'\Learning_2Gaussian_LogGMMRandTests_','\Learning_2Gaussian_LogGMMLangTests_','\Learning_2Gaussian_LogGMMRandAndLangTests_'};
Epss = (1:9)*0.2;
Dims = [5];
NbSteps = 10;
StepSize = 2000;
converged = 0; 
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
        Parameters.Epsil = 2.38^2/dim;
        
       
        Samples = {};
        AccRates = [];
        RelESSs = [];
%         Parameters.OptDens = Parameters.RealDens;
        X = [mvnrnd(Parameters.ArgMax,-Parameters.Hess^(-1),20000)];
        Parameters.Dens = gmdistribution.fit(X,1);
        
        TempPos = Parameters.ArgMax';
        CurrentModel = Parameters.Dens;
        NbComps = [];
        Ress = {};
        samples = [];
        for j = 1:NbSteps
            disp([Names{IndName} ' ' num2str(Dims(IndDims)) ' ' num2str(j)])
            
            switch IndName
                case 1
                    Parameters.Mode = 'Log';
                    Parameters.LogRatioFun = @LogRatioGMMRand;
                    Parameters.SampleFun = @SampleGMMRand;
                    Parameters.TargetAR = 0.23;
                case 2
                    Parameters.Mode = 'Log';
                    Parameters.LogRatioFun = @LogRatioGMMLang;
                    Parameters.SampleFun = @SampleGMMLang;
                    Parameters.TargetAR = 0.55;
                case 3
                    if not(converged)
                        Parameters.Mode = 'Log';
                        Parameters.LogRatioFun = @LogRatioGMMRand;
                        Parameters.SampleFun = @SampleGMMRand;
                        Parameters.TargetAR = 0.23;
                    else
                        Parameters.Mode = 'Log';
                        Parameters.LogRatioFun = @LogRatioGMMLang;
                        Parameters.SampleFun = @SampleGMMLang;
                        Parameters.TargetAR = 0.55;
                    end
            end
            
            Res = RunMCMC(TempPos,Parameters,StepSize);
            Parameters.Epsil = exp(log(Parameters.Epsil) + 0.95^j*(Res.AccRate-Parameters.TargetAR));
            samples = [Res.Vals];
            TempModels = {};
            BICs = [];
            for k  = 1:3
                try
                    TempModels{k} = gmdistribution.fit(samples',k);
                    BICs(k) = MyBIC(TempModels{k},samples');
                catch
                    BICs(k) = Inf;
                end
            end
            [NewBIC,ind] = min(BICs);
            NewModel = TempModels{ind};
            ncomps = CurrentModel.NComponents;
            CurrentBIC = MyBIC(CurrentModel,samples');
            if NewBIC<CurrentBIC-2
                CurrentModel = NewModel;
                Parameters.Dens = CurrentModel;
            else
                converged = 1;
            end
            NbComps(j) = CurrentModel.NComponents;
            Ress{j} = Res;
        end
        save([SavePath Names{IndName} num2str(dim) '.mat'],'Ress');
    end
end
        




% plot
IndName = 1;
RelESSs = {};
AccRates = {};
Lambdas = {};
Models = {};
for i = 1:length(Dims)
    load([SavePath Names{IndName} num2str(Dims(i)) '.mat']);
    RelESSs{i} = [];
    AccRates{i} = [];
    Lambdas{i} = [];
    Models{i} = {};
    for j = 1:NbSteps
        Parameters = Ress{j}.Parameters;
        TempRes = RunMCMC(Parameters.ArgMax',Parameters,50000);
        RelESSs{i}(j,:) =  TempRes.RelESS;
        AccRates{i}(j) = TempRes.AccRate;
        Lambdas{i}(j) = Parameters.Epsil;
        Models{i}{j} = Parameters.Dens;
    end
end
Descriptors = struct();
Descriptors.RelESSs = RelESSs;
Descriptors.AccRates = AccRates;
Descriptors.Lambdas = Lambdas;
Descriptors.Models = Models;
save([SavePath Names{IndName} '_descriptors.mat'],'Descriptors');



load([SavePath Names{3} '_descriptors.mat']);
dim = 1;

clf
subplot(3,1,1)
plot(0:length(Descriptors.Lambdas{dim})-1,min(Descriptors.RelESSs{dim}'),'LineWidth',2)
xlim([0 NbSteps-1])
set(gca,'XTick',0:2:NbSteps-1)
set(gca,'XTickLabel',0:2*StepSize:NbSteps*(StepSize-1))
xlabel('Iterations')
ylabel('min(RelESS)')
title(['Mixed Learning phase for the AGMM on the GMM(2) example (dim = ' num2str(Dims(dim)) ')'],'FontWeight','bold')
subplot(3,1,2)
plot(0:length(Descriptors.Lambdas{dim})-1,Descriptors.Lambdas{dim},'LineWidth',2)
xlim([0 NbSteps-1])
set(gca,'XTick',0:2:NbSteps-1)
set(gca,'XTickLabel',0:2*StepSize:NbSteps*(StepSize-1))
xlabel('Iterations')
ylabel('\lambda')
subplot(3,1,3)
xs = [];
ys = [];
for j = 1:NbSteps
    model = Descriptors.Models{dim}{j};
    ncomps = model.NComponents;
    for k = 1:ncomps
        ys = [ys randn(1,round(1000*model.PComponents(k)))*sqrt(model.Sigma(1,1,k))+model.mu(k,1)];
    end
    xs = [xs j+rand(1,length(ys)-length(xs))-1];
end
scattercloud(xs,ys)
hold on
for j = 1:NbSteps
    model = Descriptors.Models{dim}{j};
    ncomps = model.NComponents;
    for k = 1:ncomps        
        plot(j+0.001*randn(1,1)-1,model.mu(k,1),'og','MarkerEdgeColor','g','MarkerFaceColor','g')
    end
end
hold off
xlim([-0.05 NbSteps-1])
set(gca,'XTick',0:2:NbSteps-1)
set(gca,'XTickLabel',0:2*StepSize:NbSteps*(StepSize-1))
xlabel('Iterations')
ylabel('Projected mixture models')
legend('Model marginal posterior density', 'Mixture projected centers')


      
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
        Parameters.Dens = Parameters.RealDens;3
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
        
        
        
      