 function Result = EstimationSRUKFGen(Data, Model, Parameters)

 % Sarkka 2007 On Unscented Kalman Filtering for State Estimation of Continuous-Time Nonlinear Systems

 
Parameters = Model.InitializeParameters(Parameters);
 
L = length(Parameters.InitialState);
Parameters.UKFalpha.Value = 1;%1;
Parameters.UKFkappa.Value = 0;
Parameters.UKFbeta.Value = 0;




CheckParametersGen(Parameters)

if not(isfield(Parameters,'RWinEKF'))
    Parameters.RWinEKF = 0;
end


tic
NbVarsTot = Parameters.NbVariables;
Observations = Data.Observations;
ObservationInstants = Data.Instants;
NbComputingSteps = Data.NbComputingSteps;
Parameters.NbObs = length(ObservationInstants) ;

TStep     = Parameters.ComputationTStep;

InitialState = Parameters.InitialState;

L = length(InitialState);
msigm = repmat(InitialState,1,2*L+1);
initcov = 0.001;
kappa = Parameters.UKFkappa.Value;
msigm = msigm + [zeros(L,1) sqrt((L+kappa)*initcov)*eye(L) -sqrt((L+kappa)*initcov)*eye(L)];

% msigm = msigm + [zeros(L,1) 0*sqrt((L+kappa)*initcov)*eye(L) -0*sqrt((L+kappa)*initcov)*eye(L)];


Liks = 1;
LogLik = 0;
LiksMeanTraj = 1;
LogLikMeanTraj = 0;
LiksTheoret = 1;
LogLikTheoret = 0;
ms = InitialState;
W = zeros(2*L+1,1);
W(1) = kappa / (L+kappa);
for k = 1:2*L
    W(k+1) = 1/(2*(L+kappa));
end
[xm Cov] = UT(msigm,W);
Covs = [];
Covs(1,:,:)=Cov;

% disp(m(3))

for IndTime = 2:length(ObservationInstants)    
    if not(mean(Cov == Cov')==1)
        disp('pb: not sym')
    end
%     Res = Model.SRUKF_projection(Data,Model,msigm,Cov,NbComputingSteps(IndTime),IndTime,Parameters);
    Res = Model.SRUKF_projection(Data,Model,msigm,Cov,NbComputingSteps(IndTime),IndTime,Parameters);
    Model = Res.Model;
    msigm = Res.msigm;

    % p 116 Van der Merwe:
    IndObservedVar = Data.ObservedVariables(IndTime);

    [zhat,Pz] = UT(msigm(IndObservedVar,:),W);
     
    sgns = {'-','+','+'};
    [Qy,Sy] = qr([(sqrt(W(2)).*(msigm(IndObservedVar,2:end)-zhat)')' chol(Model.ObservationMeasurementNoise{IndTime})]',0);
    Sy = cholupdate(Sy,(msigm(IndObservedVar,1)-zhat)*sqrt((W(1))),sgns{sign(W(1))+2});

    [xm,Pk] = UT(msigm,W);

    n = length(xm);
    Pxz = zeros(n,1);
    for k = 1:2*n+1
       Pxz = Pxz +   W(k)*(msigm(:,k)-xm)*(msigm(IndObservedVar,k) - zhat);
    end
    
   
    
    
    
    Kk = Pxz*((Sy^-1)')*((Sy^-1));
    
   
    % same Kk
    % same Pk

    m = xm + Kk*(Observations(IndObservedVar,IndTime)-zhat);
    
   
    
    U = Kk*Sy;
    Sx = Pk;
    Sx = cholupdate(Sx,U,'-');
 
    msigm = repmat(m,1,2*L+1) + ([zeros(L,1) sqrt(L+kappa)*Sx -sqrt(L+kappa)*Sx]);
    [m,Cov] = UT(msigm,W);
    Cov
    Sx*(Sx')
    die
    
    
    Cov = Sx*(Sx');
    
    
    
    [msigm W] = SigmaPoints(m,Cov,kappa);
   
   
    if IndTime == 5
        die
    end
     
    ms(:,IndTime) = m;
    Covs(IndTime,:,:) = Cov;
    

end

LogPrior = 0;
LogCorr = 0;
Names = Parameters.Names.Estimated;

for i = 1:length(Names)
    temp = Parameters.(Names{i}).Prior(Names{i},Parameters);
    LogPrior = LogPrior +log(temp);
    temp = Parameters.(Names{i}).CorrFunct(Names{i},Parameters);
    LogCorr = LogCorr + log(temp);
%     if strcmp(Names{i},'BRmm1')
%         disp(['m' ' ' num2str(Parameters.(Names{i}).Value) ' ' num2str(log(temp))])
%     end
%     if strcmp(Names{i},'SigmaRW')
%         disp(['sigma' ' ' num2str(Parameters.(Names{i}).Value) ' ' num2str(log(temp))])
%     end
end
Result.LogPrior = LogPrior;
Result.LogCorr = LogCorr;

Result.Data = Data;
Result.Parameters = Parameters;
Result.PosteriorMeans = ms;
Result.PosteriorCovs = Covs;
Sigmas = [];
for i = 1:NbVarsTot
    Sigmas(i,:) = sqrt(Covs(:,i,i));
end
Result.Posterior975 = ms + 2*Sigmas;
Result.Posterior025 = ms - 2*Sigmas;
Result.Likelihood = prod(Liks);
Result.Liks = Liks;
if not(isreal(LogLik))
    LogLik = 0;
end
% LogLik = LogLik + log(normpdf(exp(ms(6,4))-exp(ms(6,17)),0,0.1));
Result.LogLik = LogLik;
% disp(Parameters.SigmaRW.Value)

