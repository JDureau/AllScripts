 function Result = EstimationUKFGen(Data, Model, Parameters)

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


m = InitialState;
Cov = zeros(NbVarsTot,NbVarsTot);
try  
    Ps = 0.001*diag(max(0.00001,abs(m)));
%     Ps = Parameters.InitialCov;
catch
    Ps = Cov;
end
Cov = Ps;

Liks = 1;
LogLik = 0;
LiksMeanTraj = 1;
LogLikMeanTraj = 0;
LiksTheoret = 1;
LogLikTheoret = 0;
ms = m;
Covs = [];
Covs(1,:,:)=Cov;

% disp(m(3))

for IndTime = 2:length(ObservationInstants)    
    if not(mean(Cov == Cov')==1)
        disp('pb: not sym')
    end
    Res = Model.UKF_projection(Data,Model,m,Cov,NbComputingSteps(IndTime),IndTime,Parameters);
    Model = Res.Model;
    xmtmp = Res.m;
    Pmtmp = Res.Cov;
%     if not(mean(Cov1 == Cov1')==1)
%         disp('pb: not sym')
%     end
    kappa = Parameters.UKFkappa.Value;
    [fkappi W] = SigmaPoints(xmtmp,Pmtmp,kappa)   ;
    [xkm,Pkm] = UT(fkappi,W);
    
    
%     Res.m
%     Res.Cov
%     die
    
%     
%     W = zeros(2*n+1,1);
%     W(1) = kappa / (n+kappa);
%     for k = 1:2*n
%         W(k+1) = 1/(2*(n+kappa));
%     end


    IndObservedVar = Data.ObservedVariables(IndTime);
    [zhat,Pz] = UT(fkappi(IndObservedVar,:),W,Model.ObservationMeasurementNoise{IndTime}');

    n = length(xmtmp);
    Pxz = zeros(n,1);
    for k = 1:2*n+1
       Pxz = Pxz +   W(k)*(fkappi(:,k)-xkm)*(fkappi(IndObservedVar,k) - zhat);
    end
    Kk = Pxz*(Pz^(-1));
    
  
    % Same Kk
    % Same Pmtmp
    
    m = xkm + Kk*(Observations(IndObservedVar,IndTime)-zhat);
    Cov = Pmtmp - Kk*Pz*(Kk');
    
    
    
    if IndTime == 4 
        die
    end
    
    
    
%     [V,D] = eig(Cov);
%     temp = real(eig(Cov));
%     temp = max(temp,0.001);
%     D = diag(temp);
%     vp = temp;
%     temp = ( V*D*V' + (V*D*V')')/2;
%     Cov = temp;
    
%     Cov = (Cov+Cov')/2;
%     Res.m(1)
%     Res.m(2)
%     Cov(1,1)
%     Cov(1,2)
%     Cov(2,1)
%     Cov(2,2)
%     eig(Cov)
%     Cov - Cov'
%     die
%     'xz'
%     Pxz
%     'z'
%     Pz
%     sqrt(Res.Cov(3,3))
%     Cov
%     die

    
%     if not(mean(eig(Cov)>0)==1)
%         IndTime
%         disp('pb: not def pos')
%     end
   
    [V,D] = eig(Cov);
    temp = real(eig(Cov));
    temp = max(temp,0.001);
    D = diag(temp);
    temp = ( V*D*V' + (V*D*V')')/2;
    Cov = temp;
   
     
    ms(:,IndTime) = m;
    Covs(IndTime,:,:) = Cov;
    
%  
 
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

