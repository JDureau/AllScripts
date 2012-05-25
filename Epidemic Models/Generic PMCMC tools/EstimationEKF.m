 function Result = EstimationEKF(Data, Parameters)

CheckParameters(Parameters)
    

tic
NbVarsTot = size(Data.RealData,1);
Observations = Data.Observations;
ObservationInstants = Data.Instants;
ObservedVariables = Parameters.ObservedVariables;

TStep     = Parameters.ComputationTStep;
if std(diff(ObservationInstants))>mean(diff(ObservationInstants))/100000
    disp('Attention: IrregularTSteps')
    die % For the moment, we don't deal with this. 
    pause(2)
end
ObservationTStep = mean(diff(ObservationInstants));
ComputationTStep = Parameters.ComputationTStep;
if ComputationTStep>ObservationTStep
    disp('Problem: CTStep>OTStep')
    ComputationTStep = min(ComputationTStep,ObservationTStep);
end


NbTSteps = ceil(ObservationTStep/ComputationTStep);

InitialState = Data.RealData(:,1);

m = InitialState;
Cov = zeros(4,4);
Ps = eye(4);
% initilization of Cov:
% mpred = m;
% for IndDiscr = 1:NbTSteps
%     mpred(1) = mpred(1) + (-mpred(4)*mpred(1)*mpred(2)+Parameters.Beta*mpred(3))*TStep;
%     mpred(2) = mpred(2) + ( mpred(4)*mpred(1)*mpred(2)-Parameters.Gamma*mpred(2))*TStep;
%     mpred(3) = mpred(3) + (Parameters.Gamma*mpred(2) - Parameters.Beta*mpred(3))*TStep;
%     mpred(4) = mpred(4);
%     Jacobian = zeros(4,4);
%     Jacobian(1,1) = -mpred(4)*mpred(2);
%     Jacobian(1,2) = -mpred(4)*mpred(1);
%     Jacobian(1,3) = Parameters.Beta;
%     Jacobian(1,4) = -mpred(1)*mpred(2);
%     Jacobian(2,1) =  mpred(4)*mpred(2);
%     Jacobian(2,2) =  mpred(4)*mpred(1)- Parameters.Gamma;
%     Jacobian(2,3) =  0;
%     Jacobian(2,4) =  mpred(1)*mpred(2);
%     Jacobian(3,1) =  0;
%     Jacobian(3,2) =  Parameters.Gamma;
%     Jacobian(3,3) =  - Parameters.Beta;
%     Q = zeros(4,4);
%     Q(4,4) = (Parameters.SigmaDiffusion(4))^2;
%     exp = expm(IndDiscr*TStep*Jacobian);
%     Cov = Cov + exp*Q*exp'*TStep;
% end
Cov = Ps*Cov;

Liks = 1;
LogLik = 0;
LiksMeanTraj = 1;
LogLikMeanTraj = 0;
LiksTheoret = 1;
LogLikTheoret = 0;
ms = m;
Covs = [];
Covs(1,:,:)=Cov;
for IndTime = 2:length(ObservationInstants)
    
    mpred = m;
    for IndDiscr = 1:NbTSteps
        mtemp = mpred;
        mpred(1) = mpred(1) + (-mtemp(4)*mtemp(1)*mtemp(2)+ Parameters.Beta*mtemp(3))*TStep;
        mpred(2) = mpred(2) + ( mtemp(4)*mtemp(1)*mtemp(2)- Parameters.Gamma*mtemp(2))*TStep;
        mpred(3) = mpred(3) + ( Parameters.Gamma*mtemp(2) - Parameters.Beta*mtemp(3))*TStep;
        if strcmp(Parameters.DiffusionType,'Add')
            mpred(4) = mtemp(4);
        elseif strcmp(Parameters.DiffusionType,'OUD')
            mpred(4) = mtemp(4) + ( Parameters.KappaOU*(Parameters.LogExpMuOU-mtemp(4)))*TStep;
        end
        
        Jacobian = zeros(4,4);
        Jacobian(1,1) = -mpred(4)*mpred(2);
        Jacobian(1,2) = -mpred(4)*mpred(1);
        Jacobian(1,3) = Parameters.Beta;
        Jacobian(1,4) = -mpred(1)*mpred(2);
        Jacobian(2,1) =  mpred(4)*mpred(2);
        Jacobian(2,2) =  mpred(4)*mpred(1)- Parameters.Gamma;
        Jacobian(2,3) =  0;
        Jacobian(2,4) =  mpred(1)*mpred(2);
        Jacobian(3,1) =  0;
        Jacobian(3,2) =  Parameters.Gamma;
        Jacobian(3,3) =  - Parameters.Beta;
        if strcmp(Parameters.DiffusionType,'OUD')
            Jacobian(4,4) = -Parameters.KappaOU;
        end
        Q = zeros(4,4);
        if strcmp(Parameters.DiffusionType,'Add')
            Q(4,4) = (Parameters.SigmaRW)^2;
        elseif strcmp(Parameters.DiffusionType,'OUD')
            Q(4,4) = (Parameters.SigmaOU)^2;
        end
        Cov = Cov + (Jacobian*Cov+Cov*Jacobian'+Q)*TStep;
    end
    Ck = zeros(1,4); % Jacobian of obs
    Ck(1,2) = 1;
    
    % new attempt
    predCov = Cov;
    ypred = mpred(2);
    vk = Observations(2,IndTime) - ypred; 
    Rk = (Parameters.SigmaObs(2)*Observations(2,IndTime))^2;
    Sk = Ck*Cov*Ck' + Rk;
    Kk = Cov*Ck'*Sk^-1;
    m = mpred + Kk*vk;
    Cov = Cov - Kk*Sk*Kk';
   
%     % old
%     for IndObservedVar = 1:length(Parameters.ObservedVariables) 
%         Rk_km1 = Ck*Cov*Ck'+(Parameters.SigmaObs(ObservedVariables(IndObservedVar))*Observations(ObservedVariables(IndObservedVar),IndTime))^2;
%     end
%     Kk = Cov*Ck'*Rk_km1^-1;
%     ypred = mpred(2);
%     res = Observations(2,IndTime) - ypred; 
%     m = mpred + Kk*res;
%     Cov = Cov - Kk*Rk_km1*Kk';
%     
    ms(:,IndTime) = m;
    Covs(IndTime,:,:) = Cov;
    
    
    % Likelihood of the measurement, in order to finally compute p(y|pars)
    % (formula taken from Sarkka's documentation for his EKF library.
%     % old
%     tempLik = normpdf(res,0,Rk_km1);
    
    % new
%     tempsamples = Observations(2,IndTime) + Observations(2,IndTime)*Parameters.SigmaObs(2)*randn(1,10000);
%     tempLik = mean(normpdf(tempsamples,mpred(2),sqrt(predCov(2,2))));

    c1 = 1/(2*(Observations(2,IndTime)*Parameters.SigmaObs(2))^2)+1/(2*predCov(2,2));
    c2 = Observations(2,IndTime)/(2*(Observations(2,IndTime)*Parameters.SigmaObs(2))^2)+mpred(2)/(2*predCov(2,2));
    c3 = Observations(2,IndTime)^2/(2*(Observations(2,IndTime)*Parameters.SigmaObs(2))^2)+mpred(2)^2/(2*predCov(2,2));

    tempLik = 1/(sqrt(2*pi)*Observations(2,IndTime)*Parameters.SigmaObs(2)*sqrt(predCov(2,2))*sqrt(2*c1))*exp(c2^2/c1-c3);

%     
% %     old
%     tempLik = normpdf(vk,0,Sk);
    
    Liks(IndTime) = tempLik;
    LogLik = LogLik + log(tempLik);
    tempLik = normpdf(m(2),Observations(2,IndTime),Observations(2,IndTime)*Parameters.SigmaObs(2));
    LiksMeanTraj(IndTime) = tempLik;
    LogLikMeanTraj = LogLikMeanTraj + log(tempLik);
    
    tempLik = normpdf(Observations(2,IndTime),Observations(2,IndTime),Observations(2,IndTime)*Parameters.SigmaObs(2));
    LiksTheoret(IndTime) = tempLik;
    LogLikTheoret = LogLikTheoret + log(tempLik);
    
%     InvCov = Cov^-1;
%     InvCov2 = InvCov;
%     o = Observations(2,IndTime);
%     InvCov2(5,5) = 1/(Parameters.SigmaObs(2)*o)^2;
%     InvCov2(2,2) = InvCov2(2,2) + 1/(Parameters.SigmaObs(2)*o)^2;
%     InvCov2(2,5) = -1/(Parameters.SigmaObs(2)*o)^2;
%     InvCov2(5,2) = -1/(Parameters.SigmaObs(2)*o)^2;
%     
%     Cov2 = InvCov2^-1;
%     Liks(IndTime) = normpdf(o,mpred(2),sqrt(Cov2(5,5)));
%     if Liks == 0
%         IndTime
%     end
%     LogLik = LogLik + log(normpdf(o,mpred(2),sqrt(Cov2(5,5))));
%     if isnan(LogLik)
%         IndTime
%     end
end

Result.Data = Data;
Result.Parameters = Parameters;
Result.PosteriorMeans = ms;
Result.PosteriorCovs = Covs;
Sigmas = [];
for i = 1:4
    Sigmas(i,:) = sqrt(Covs(:,i,i));
end
Result.Posterior975 = ms + 2*Sigmas;
Result.Posterior025 = ms - 2*Sigmas;
Result.Likelihood = prod(Liks);
Result.Liks = Liks;
Result.LiksMeanTraj = LiksMeanTraj;
Result.LiksTheoret = LiksTheoret;
if not(isreal(LogLik))
    LogLik = 0;
end
Result.LogLik = LogLik;
Result.RealTime = toc;