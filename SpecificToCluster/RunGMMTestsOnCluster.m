function [] = RunGMMTestsOnCluster(IndDensity,IndMethod,IndLogOrNot,dim,ind)

ind = ind-floor(ind/14)*14;

Methods =  {'MALA','LocalMALA','RMALA','GMMRand','GMMLang'};
Densities = {'GMM2','Banana'};

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
Res = RunMCMC(Parameters.ArgMax',Parameters,1000);
Res.Samples = Res.Vals;
Res.Eps = Parameters.Epsil;

SavePath = '/users/ecologie/dureau/src/AllData/GMM/';
save([SavePath Densities{IndDensity} '_' Methods{IndMethod} '_dim' num2str(dim) '_Log' num2str(IndLogOrNot) '_eps' num2str(Parameters.Epsil) '.mat'],'Res');

  



