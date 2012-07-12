function [] = RunGMMTestsOnCluster(IndDensity,IndMethod,IndLogOrNot,dim,ind)

ind = ind-floor(ind/14)*14;

Methods =  {'MALA','LocalMALA','GMCovMALA','GMMRand','GMMLang','GMHMC','GMCovHMC','GMind'};
Densities = {'GMM','GMM2','Banana'};

s = RandStream('mcg16807','Seed',sum(fix(clock)))
RandStream.setDefaultStream(s)

%IndMethod = 4;
%IndDensity = 1;
%IndLogOrNot = 1;
%cdim = 2; % min 2


cd('/users/ecologie/dureau/src/AllScripts/')
addpath([pwd '/MCMCTests/'])
addpath([pwd '/General Tools/'])
addpath([pwd '/Toolboxes/mcmcdiag/'])
addpath([pwd '/Toolboxes/'])


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
        Parameters.f = @fBanana;
        B = 0.1;
        X = mvnrnd(zeros(dim,1),eye(dim),100000);
        X(:,1) = 10*X(:,1);
        X(:,2) = X(:,2)-B*X(:,1).^2+100*B;
        Parameters.TrueSamples = X;
        Parameters.RealDens = gmdistribution.fit(X,60);
        scattercloudGM(X(:,1),X(:,2),Parameters.RealDens)
        Parameters.OptDens = Parameters.RealDens;
        Parameters.B = B;
        Parameters.Dim = dim;
        [b,ind] = max(Parameters.RealDens.PComponents);
        Parameters = FindFisherInfMat(Parameters.RealDens.mu(ind,:),Parameters) ;
        test = 0;
%         Parameters.Epsil = Parameters.Epsil/10;
        while not(test)
            try
                Parameters.Dens = Parameters.RealDens;%gmdistribution.fit(X,3);
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
        Parameters.LogRatioFun = @LogRatioGMCovMALA;
        Parameters.SampleFun = @SampleGMCovMALA;
    case 4
        Parameters.LogRatioFun = @LogRatioGMMRand;
        Parameters.SampleFun = @SampleGMMRand;
    case 5
        Parameters.LogRatioFun = @LogRatioGMMLang;
        Parameters.SampleFun = @SampleGMMLang;
    case 6
        Parameters.LogRatioFun = @LogRatioGMHMC;
        Parameters.SampleFun = @SampleGMHMC;
        Parameters.ScalingCov = -(Parameters.Hess^-1);
        Parameters.ArgMax = [Parameters.ArgMax Parameters.ArgMax];
    case 7
        Parameters.LogRatioFun = @LogRatioGMCovHMC;
        Parameters.SampleFun = @SampleGMCovHMC;
        Parameters.ScalingCov = -(Parameters.Hess^-1);
        Parameters.ArgMax = [Parameters.ArgMax Parameters.ArgMax];
    case 8
        Parameters.LogRatioFun = @LogRatioGMMind;
        Parameters.SampleFun = @SampleGMMind;
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
Res = RunMCMC(Parameters.ArgMax',Parameters,50);
Res.Samples = Res.Vals;
Res.Eps = Parameters.Epsil;

SavePath = '/users/ecologie/dureau/src/AllData/GMM/';
save([SavePath Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat'],'Res');

  



