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



Epss = (1:14)*0.2;
Parameters.Epsil = Epss(ind+1);


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
        load([SavePath '/BananaModelParameters100.mat'])
%         Parameters.f = @fBanana;
%         B = 0.1;
%         X = mvnrnd(zeros(dim,1),eye(dim),100000);
%         X(:,1) = 10*X(:,1);
%         X(:,2) = X(:,2)-B*X(:,1).^2+100*B;
%         Parameters.TrueSamples = X;
%         Parameters.RealDens = gmdistribution.fit(X,100);
%         scattercloudGM(X(:,1),X(:,2),Parameters.RealDens)
%         Parameters.OptDens = Parameters.RealDens;
%         Parameters.B = B;
%         Parameters.Dim = dim;
%         [b,ind] = max(Parameters.RealDens.PComponents);
%         Parameters = FindFisherInfMat(Parameters.RealDens.mu(ind,:),Parameters) ;
%         test = 0;
% %         Parameters.Epsil = Parameters.Epsil/10;
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
        NbIterations = 100000;
        if IndLogOrNot
            die
        end
    case 2
        Parameters.LogRatioFun = @LogRatioGMCovMALA;
        Parameters.SampleFun = @SampleGMCovMALA;
        NbIterations = 100000;
        if IndLogOrNot
            die
        end
    case 3
        Parameters.LogRatioFun = @LogRatioGMMRand;
        Parameters.SampleFun = @SampleGMMRand;
        NbIterations = 100000;
    case 4
        Parameters.LogRatioFun = @LogRatioGMMLang;
        Parameters.SampleFun = @SampleGMMLang;
        NbIterations = 100000;
    case 5
        Parameters.LogRatioFun = @LogRatioHMC;
        Parameters.SampleFun = @SampleHMC;
        Parameters.ScalingCov = -(Parameters.Hess^-1);
        Parameters.ArgMax = [Parameters.ArgMax Parameters.ArgMax];
        if IndDensity == 3
            Parameters.fGrad = @ComputeBananaGrad;
        else
            die
        end
        if IndLogOrNot
            die
        end
        NbIterations = 30000;
        Parameters.Epsil = Parameters.Epsil/20;
    case 6
        Parameters.LogRatioFun = @LogRatioGMCovHMC;
        Parameters.SampleFun = @SampleGMCovHMC;
        Parameters.ScalingCov = -(Parameters.Hess^-1);
        Parameters.ArgMax = [Parameters.ArgMax Parameters.ArgMax];
        if IndDensity == 3
            Parameters.fGrad = @ComputeBananaGrad;
        else
            die
        end
        if IndLogOrNot
            die
        end
        NbIterations = 30000;
        Parameters.Epsil = Parameters.Epsil/20;
    case 7
        Parameters.LogRatioFun = @LogRatioGMCovHMC;
        Parameters.SampleFun = @SampleGMCovHMC;
        Parameters.ScalingCov = -(Parameters.Hess^-1);
        Parameters.ArgMax = [Parameters.ArgMax Parameters.ArgMax];
        if IndDensity == 3
            Parameters.fGrad = @ComputeGMMGrad;
        else
            die
        end
        if IndLogOrNot
            die
        end
        NbIterations = 30000;
        Parameters.Epsil = Parameters.Epsil/20;
    case 8
        Parameters.LogRatioFun = @LogRatioGMMind;
        Parameters.SampleFun = @SampleGMMind;
        NbIterations = 100000;
        Parameters.Epsil = Parameters.Epsil;
        if or(IndLogOrNot, ind>1)
            die
        end
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
Res = RunMCMC(Parameters.ArgMax',Parameters,NbIterations);

SavePath = '/users/ecologie/dureau/src/AllData/GMM/';
save([SavePath Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat'],'Res');

  



