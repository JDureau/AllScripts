function Result = ParameterEstimation_GMM(Data,Parameters)




% Initialise TempPar
TempPar = ProposeInitialParameter(Data, Parameters);

% Calibrate the method
[Parameters, TempPar] = CalibrateMethod( Data, Parameters, TempPar);
[Parameters, TempPar] = CalibrateMethod( Data, Parameters, TempPar);
[Parameters, TempPar] = CalibrateMethod( Data, Parameters, TempPar);

% Burn-In
Res = RunEstimationMethod(Data,Parameters,TempPar,Parameters.NbItBurnIn);
TempPar = Res.TempPar;

% Calibrate the method after BurnIn
[Parameters, TempPar] = CalibrateMethod( Data, Parameters, TempPar);
[Parameters, TempPar] = CalibrateMethod( Data, Parameters, TempPar);
[Parameters, TempPar] = CalibrateMethod( Data, Parameters, TempPar);

% Run The algorithm
Res = RunEstimationMethod(Data,Parameters,TempPar,Parameters.NbIterations);

% The results from previous algorithms might not be perfect. We'll exploit
% them to do something better.


SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
load([SavePath '\2010_12_06_Result_RWSigmaRWBeta.mat'])
ThetasSamples = Result_RWSigmaRWBeta.Thetas;
LogThetasSamples = log(Result_RWSigmaRWBeta.Thetas);
LogThetasSamples = Res.LogThetas;
% Result_RWSigmaRWBeta


% LogThetasSamples = Res.LogThetas;
Parameters = DefineScalingPars(LogThetasSamples,Parameters);
ScaledLogThetasSamples = Scale(LogThetasSamples,Parameters);


Parameters.GMeth = 'GMM';
Parameters.DensityModel = FindBestDensityModel(ScaledLogThetasSamples',Parameters,1);
Parameters.DensityModel = FindBestDensityModel(ScaledLogThetasSamples');

Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'GMM';


Parameters.aim = 0.23;
[Parameters, TempPar] = CalibrateMethod( Data, Parameters, TempPar);
Res = RunEstimationMethod(Data,Parameters,TempPar,10000);

PlotRun(Res)



ndims = 2;
Parameters.MCMCType = 'Inde';
Parameters.PosteriorDensityModel = gmdistribution.fit(ScaledLogThetasSamples',1);
VPs = [];
LogLiksOld = [];
LogLiksNew = [];
for IndDim = 1:ndims
    Condition = 0;
    it = 0;
    while not(Condition)
        it = it+1;
        disp(it)
        Cov = Parameters.PosteriorDensityModel.Sigma;
        DensityModel = Parameters.PosteriorDensityModel;
        [V,D] = eig(Cov);
        EigValues = [];
        for i = 1:ndims
            EigValues(i) = D(i,i);
        end
        [EigValues,inds] = sort(EigValues,'descend');
        VPs(IndDim,it,:) = EigValues;
        V = V(:,inds);
        NewEigValues = EigValues;
        NewEigValues(IndDim) = NewEigValues(IndDim)*1.21;
        NewD = diag(NewEigValues);
        NewCov = V*NewD*V^-1;
        mu = DensityModel.mu;
        NewDensityModel = gmdistribution(mu,NewCov,1);
        NewParameters = Parameters;
        NewParameters.PosteriorDensityModel = NewDensityModel ;
        NewRes = RunEstimationMethod(Data,NewParameters,TempPar,1000);
        ScaledNewThetas = Scale(NewRes.Thetas,Parameters);
        scatter(ScaledNewThetas(1,:),ScaledNewThetas(2,:))
        hold on
        ezcontour(@(x,y)pdf(NewDensityModel,[x y]),[-8 6],[-8 6])
        hold off
        pause(0.1)
        LogLiksOld(it) = sum(log(pdf(DensityModel,ScaledNewThetas')));
        LogLiksNew(it) = sum(log(pdf(NewDensityModel,ScaledNewThetas')));
        
        if sum(log(pdf(DensityModel,ScaledNewThetas'))) > sum(log(pdf(NewDensityModel,ScaledNewThetas')))
            Condition = 1;            
        end
        Parameters = NewParameters;  
        Parameters.PosteriorDensityModel = gmdistribution.fit(ScaledNewThetas',1);
    end
end

Res = RunEstimationMethod(Data,Parameters,TempPar,6000);
ScaledData = Scale(Res.Thetas,Parameters)';
Parameters.DensityModel = FindBestDensityModel(ScaledData);

Res2 = RunEstimationMethod(Data,Parameters,TempPar,6000);

ScaledData = Scale(Res2.Thetas,Parameters)';
Parameters2 = Parameters;
% Parameters2.DensityModel = FindBestDensityModel(ScaledData);

Res3 = RunEstimationMethod(Data,Parameters2,TempPar,20000);

Parameters2.DensityModel = gmdistribution.fit(ScaledData,2);
Covs = Parameters2.DensityModel.Sigma;
Mus = Parameters2.DensityModel.mu;

NewMus = [];
NewCovs = [];
for i = 1:Parameters2.DensityModel.NComponents
    NewMus(i,:) = Mus(i,:)*diag(Parameters2.ScalingStds) + Parameters2.ScalingMeans;
    temp = squeeze(Covs(:,:,i));
    NewCovs(:,:,i) = diag(Parameters2.ScalingStds)*temp*diag(Parameters2.ScalingStds);
end
NewDensityModel = gmdistribution(NewMus,NewCovs,Parameters2.DensityModel.PComponents);

Parameters2.RealScaleDensityModel = NewDensityModel ;
ezsurf(@(x,y)pdf(NewDensityModel,[x y]),[0 20],[0.5 1.5]*10^-5)

x = Res2.Thetas(1,:);
y = Res2.Thetas(2,:);
u = ksdensity(x,x,'function','cdf');
v = ksdensity(y,y,'function','cdf');



NbIts = 1000;

Parameters2.MCMCType = 'Lang';
Parameters2.G = 'PseudoCst';
Parameters2.Epsil = 1;

Parameters2.NoDerivatives = 1;
TempPar = ProposeInitialParameter(Data, Parameters2);
[Parameters2, TempPar] = CalibrateMethod( Data, Parameters2, TempPar);
ResPseudo = RunEstimationMethod(Data,Parameters2,TempPar,NbIts);

Parameters2.MCMCType = 'Rand';
Parameters2.Epsil = 1;

Parameters2.RWsamplCov = cov(Res2.Thetas');
TempPar = ProposeInitialParameter(Data, Parameters2);
Parameters2.ComputeRWsamplCov = 0;
[Parameters2, TempPar] = CalibrateMethod( Data, Parameters2, TempPar);
ResOne = RunEstimationMethod(Data,Parameters2,TempPar,NbIts);


mu = [1 2];
sigma = diag([2,3]);
p = 1;
dis = gmdistribution(mu,sigma,p);
Pars = Parameters2;
Pars.RealScaleDensityModel = dis;

clf
x1 = -2:0.1:6;
x2 = 0*ones(size(x1));
y = log(pdf(dis,[x1' x2']));
plot(x1,y)
hold on
Pars.Pars = [5 x2(1)]';
grad = ComputeFirstDerivative(Pars);
xis = [Pars.Pars(1)-1 Pars.Pars(1) Pars.Pars(1)+1];
yis = [log(pdf(dis,Pars.Pars'))-grad(1) log(pdf(dis,Pars.Pars')) log(pdf(dis,Pars.Pars'))+grad(1)];  
plot(xis,yis,'r')
hold off

clf
i = 1;
j = 2;
x1 = -6:0.1:6;
x2 = 0*ones(size(x1));
y = [] ;
for k = 1:length(x1)
    Pars.Pars = [x1(k) x2(1)]';
    y(k,:) = ComputeFirstDerivative(Pars);
end
plot(x1,y(:,j))
hold on
ind =83;
Pars.Pars = [x1(ind) x2(1)]';
grad = ComputeSecondDerivative(Pars);
xis = [Pars.Pars(1)-1 Pars.Pars(1) Pars.Pars(1)+1];
yis = [y(ind,j)-grad(i,j) y(ind,j) y(ind,j)+grad(i,j)];  
plot(xis,yis,'r')
hold off


