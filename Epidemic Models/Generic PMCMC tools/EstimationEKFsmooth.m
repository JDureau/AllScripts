function Result = EstimationEKFsmooth(Data, Parameters)

% Extended Kalman smoother
% Formulas come from p95 of Sarkka's phd.
% To work, this requires cov to be invertible, so the dynamic of S I R to
% be stochastic also.


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

m_x = InitialState ;
Cov = zeros(4,4);
Ps = eye(4);
Q = diag(Parameters.SigmaDiffusion.^2);
% initilization of Cov:
% mpred = m_x;
% CovsPred(1,:,:) = Cov;
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
%     
%     exp = expm(IndDiscr*TStep*Jacobian);
%     Cov = Cov + exp*Q*exp'*TStep;
%     CovsPreds(:,:,1,IndDiscr) = Cov;
% end
% Cov = Ps*Cov;


ms_x = m_x;
mpreds_x = m_x;
Covs = [];
Covs(1,:,:) = Cov;

f_rec = [];
dets = [];
LogLikTheoret = 0;
LogLikMeanTrajFilter = 0;
LiksTheoret = normpdf(Observations(2,1),Observations(2,1),Observations(2,1)*Parameters.SigmaObs(2));
LiksMeanTrajFilter = normpdf(Observations(2,1),Observations(2,1),Observations(2,1)*Parameters.SigmaObs(2));
for IndTime = 1:length(ObservationInstants)-1
    
    mpred = m_x;
    for IndDiscr = 1:NbTSteps
        f=[];
        f(1) = (-mpred(4)*mpred(1)*mpred(2)+Parameters.Beta*mpred(3));
        f(2) = ( mpred(4)*mpred(1)*mpred(2)-Parameters.Gamma*mpred(2));
        f(3) = (Parameters.Gamma*mpred(2) - Parameters.Beta*mpred(3));
        f(4) = 0;
        mpred(1) = mpred(1) + f(1)*TStep;
        mpred(2) = mpred(2) + f(2)*TStep;
        mpred(3) = mpred(3) + f(3)*TStep;
        mpred(4) = mpred(4);
        f_rec(:,IndTime,IndDiscr) = f;
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
        Jacs_rec(:,:,IndTime,IndDiscr) = Jacobian;
        
        Cov = Cov + (Jacobian*Cov+Cov*Jacobian'+Q)*TStep;
        mpreds_x(:,IndTime,IndDiscr) = mpred;
        CovsPreds(:,:,IndTime,IndDiscr) = Cov;
        
    end
    Ck = zeros(1,4); % Jacobian of obs
    Ck(1,2) = 1;
    for IndObservedVar = 1:length(Parameters.ObservedVariables)    
        Rk_km1 = Ck*Cov*Ck'+(Parameters.SigmaObs(ObservedVariables(IndObservedVar))*Observations(ObservedVariables(IndObservedVar),IndTime+1))^2;
    end
    Kk = Cov*Ck'*Rk_km1^-1;
    ypred = mpred(2);
    res = Observations(2,IndTime+1) - ypred; 
    m_x = mpred + Kk*res;
    Cov = Cov - Kk*Rk_km1*Kk';
    dets(IndTime) = det(Cov);
    ms_x(:,IndTime+1) = m_x;
    Covs(IndTime+1,:,:) = Cov;
    
    tempLik = normpdf(Observations(2,IndTime+1),Observations(2,IndTime+1),Observations(2,IndTime+1)*Parameters.SigmaObs(2));
    LiksTheoret(IndTime+1) = tempLik;
    LogLikTheoret = LogLikTheoret + log(tempLik);
    
    tempLik = normpdf(m_x(2),Observations(2,IndTime+1),Observations(2,IndTime+1)*Parameters.SigmaObs(2));
    LiksMeanTrajFilter(IndTime+1) = tempLik;
    LogLikMeanTrajFilter = LogLikMeanTrajFilter + log(tempLik);
    
end

ms_s(:,size(ms_x,2)) = ms_x(:,end);
m_s = ms_x(:,end);
Covs_s = [];
Covs_s(size(ms_x,2),:,:) = Covs(size(ms_x,2),:,:);
Q = diag(Parameters.SigmaDiffusion.^2);
ReconstructedLambdas = [];
LogLikMeanTrajSmoother = 0;
LiksMeanTrajSmoother(length(ObservationInstants)) = LiksMeanTrajFilter(end); 

for IndTime = length(ObservationInstants)-1:-1:1
    for IndDiscr = NbTSteps:-1:1
        P = squeeze(CovsPreds(:,:,IndTime,IndDiscr));
        F = squeeze(Jacs_rec(:,:,IndTime,IndDiscr));
        mpred = mpred - TStep*( f_rec(:,IndTime,IndDiscr) + (F*P + Q)*P^-1*(mpred-mpreds_x(:,IndTime,IndDiscr)));
        Cov = Cov - TStep*( (F*P+Q)*P^-1*Cov + Cov*P^-1*(P*F'+Q) - Q );
        ReconstructedLambdas = [mpred(4) ReconstructedLambdas]; 
    end
    m_s = mpred;
    ms_s(:,IndTime) = m_s;
    Covs_s(IndTime,:,:) = Cov;
    
    tempLik = normpdf(m_s(2),Observations(2,IndTime),Observations(2,IndTime)*Parameters.SigmaObs(2));
    LiksMeanTrajSmoother(IndTime) = tempLik;
    LogLikMeanTrajSmoother = LogLikMeanTrajSmoother + log(tempLik);
end
LiksMeanTrajSmoother(1) = normpdf(Observations(2,1),Observations(2,1),Observations(2,1)*Parameters.SigmaObs(2));

ms_s(:,1) = ms_x(:,1);
Covs_s(1,:,:) = Covs(1,:,:);


Result.LiksTheoret = LiksTheoret;
Result.LiksMeanTrajFilter = LiksMeanTrajFilter;
Result.LiksMeanTrajSmoother = LiksMeanTrajSmoother;
Result.Data = Data;
Result.Parameters = Parameters;
Result.PosteriorMeansFilter = ms_x;
Result.PosteriorCovsFilter = Covs;
Result.PosteriorMeansSmoother = ms_s;
Result.PosteriorCovsSmoother = Covs_s;
Sigmas = [];
for i = 1:4
    Sigmas(i,:) = sqrt(Covs_s(:,i,i));
end
Result.Posterior975 = ms_s + 2*Sigmas;
Result.Posterior025 = ms_s - 2*Sigmas;
Result.ReconstructedLambdas = ReconstructedLambdas;
% plot(Data.RealData(4,:),'g')
% hold on
% plot(ms_s(4,1:end))
% plot(ms_x(4,1:end),'k')
% % plot(ResultSMC.PosteriorMeans(4,:))
% plot(ResultSMC.Posterior975(4,:),'r')
% plot(ResultSMC.Posterior025(4,:),'r')
% hold off

temp = 1;
for IndObservedVar = 1:length(Parameters.ObservedVariables)    
    Result.Likelihood = prod(normpdf(ms_s(2,:),Observations(ObservedVariables(IndObservedVar),:),Observations(ObservedVariables(IndObservedVar),:)*Parameters.SigmaObs(ObservedVariables(IndObservedVar))));
end
Result.RealTime = toc;