function [] = RunGMMTestsOnCluster(IndDensity,IndMethod,IndLogOrNot,dim,ind)

ind = ind-floor(ind/14)*14;

Methods =  {'MALA','GMCovMALA','GMRand','GMLang','HMC','HMCGMCov','HMCGMCovGrad','GMind'};
Densities = {'GMM','GMM2','Banana'};



SavePath = '/users/ecologie/dureau/src/AllData/GMM/';
s = RandStream('mcg16807','Seed',sum(fix(clock)))
RandStream.setDefaultStream(s)

%IndMethod = 4;
%IndDensity = 1;
%IndLogOrNot = 1;
%cdim = 2; % min 2


cd('/users/ecologie/dureau/src/AllScripts/')
addpath([pwd '/MCMCTests/'])
addpath([pwd '/General Tools/'])
addpath([pwd '/Toolboxes'])





switch IndDensity
    case 1
        Parameters.f = @fGMM;
        mu1 = zeros(dim,1);
        sigma  = eye(dim);
        X = [mvnrnd(mu1,sigma,20000)']';
        Parameters.RealDens = gmdistribution.fit(X,1);
        Parameters.Dim = dim;
        Parameters = FindFisherInfMat(zeros(1,dim),Parameters) ;
        Parameters.OptDens = Parameters.RealDens;
        Parameters.Dens = gmdistribution.fit(X,1);
        Parameters.G = @GGMM;
        Parameters.GDerivs = @GDerivsGMM;
    case 2
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
    case 3
        load([SavePath '/BananaModelParameters15.mat'])
        Parameters.f = @fBanana;
        Parameters.fGrad = @ComputeGMMGrad;
%         Parameters.f = @fGMM;

%         B = 0.1;
%         X = mvnrnd(zeros(dim,1),eye(dim),50000);
%         X(:,1) = 10*X(:,1);
%         X(:,2) = X(:,2)-B*X(:,1).^2+100*B;
%         Parameters.TrueSamples = X;
%         Parameters.RealDens = gmdistribution.fit(X,15);
%         scattercloudGM(X(:,1),X(:,2),Parameters.RealDens)
%         Parameters.OptDens = Parameters.RealDens;
%         Parameters.B = B;
%         Parameters.Dim = dim;
%         [b,ind] = max(Parameters.RealDens.PComponents);
        Parameters = FindFisherInfMat(Parameters.RealDens.mu(ind,:),Parameters) ;
%         test = 0;
%         Parameters.Epsil = Parameters.Epsil/10;
%         while not(test)
%             try
%                 Parameters.Dens = Parameters.RealDens;%gmdistribution.fit(X,3);
%                 test = 1;
%             end
%         end
end
       
switch IndMethod
    case 1
        Parameters.LogRatioFun = @LogRatioMALA;
        Parameters.SampleFun = @SampleMALA;
        Parameters.ScalingCov = -(Parameters.Hess^-1);
%         Parameters.OptAR = 0.7;
        EpsMin = 0.1;
        EpsMax = 0.5;
        NbIterations = 50000;
    case 2
        Parameters.LogRatioFun = @LogRatioGMCovMALA;
        Parameters.SampleFun = @SampleGMCovMALA;
%         Parameters.OptAR = 0.7;
        EpsMin = 1;
        EpsMax = 1.5;
        NbIterations = 50000;
    case 3
        Parameters.LogRatioFun = @LogRatioGMMRand;
        Parameters.SampleFun = @SampleGMMRand;
%         Parameters.OptAR = 0.23;
        EpsMin = 2;
        EpsMax = 5;
        NbIterations = 50000;
    case 4
        Parameters.LogRatioFun = @LogRatioGMMLang;
        Parameters.SampleFun = @SampleGMMLang;
%         Parameters.OptAR = 0.7;
        EpsMin = 1;
        EpsMax = 1.5;
        NbIterations = 50000;
%     case 5
%         Parameters.LogRatioFun = @LogRatioHMC;
%         Parameters.SampleFun = @SampleHMC;
%         Parameters.ScalingCov = -(Parameters.Hess^-1);
%         Parameters.ArgMax = [Parameters.ArgMax Parameters.ArgMax];
%         Parameters.OptAR = 0.8;
%         if IndDensity == 3
%             Parameters.fGrad = @ComputeBananaGrad;
% %             Parameters.fGrad = @ComputeGMMGrad;
%         else
%             die
%         end
%         if IndLogOrNot
%             die
%         end
%         NbIterations = 30000;
%         Parameters.Epsil = Parameters.Epsil/20;
    case 5
        Parameters.LogRatioFun = @LogRatioGMCovHMC;
        Parameters.SampleFun = @SampleGMCovHMC;
        Parameters.ScalingCov = -(Parameters.Hess^-1);
        Parameters.ArgMax = [Parameters.ArgMax Parameters.ArgMax];        
        NbIterations = 20000;
        EpsMin = 0.05;
        EpsMax = 0.2;
    case 7
        Parameters.LogRatioFun = @LogRatioGMMind;
        Parameters.SampleFun = @SampleGMMind;
        NbIterations = 100000;
        Parameters.Epsil = Parameters.Epsil;
        if or(IndLogOrNot, ind>1)
            die
        end
end

Epss = EpsMin:(EpsMax-EpsMax)/10:EpsMax;
Parameters.Epsil = Epss(ind+1);

switch IndLogOrNot 
    case 0
        Parameters.Mode = 'L';
    case 1
        Parameters.Mode = 'Log';
end


Samples = {};
AccRates = [];
RelESSs = [];
Res = RunMCMC(Parameters.ArgMax',Parameters,NbIterations);

SavePath = '/users/ecologie/dureau/src/AllData/GMM/';
save([SavePath Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat'],'Res');

  



