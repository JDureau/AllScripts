% Data
Parameters.ComputationTStep = 0.1;
Parameters.ObsLength = 40;
Parameters.Sigma.RealValue = 0.1;
Parameters.SigmaObs = 0.05;

n = Parameters.ObsLength/Parameters.ComputationTStep;
Data.Beta = cumsum(randn(1,n)*Parameters.Sigma.RealValue*sqrt(Parameters.ComputationTStep));
Data.NbComputingSteps = [0 ones(1,Parameters.ObsLength-1)/Parameters.ComputationTStep];
Data.Observations = Data.Beta(cumsum(ones(1,Parameters.ObsLength-1)/Parameters.ComputationTStep))+Parameters.SigmaObs*randn(1,Parameters.ObsLength-1);

plot(Data.Beta)
hold on
plot(cumsum(ones(1,Parameters.ObsLength-1)/Parameters.ComputationTStep),Data.Observations,'g')
hold off

NbIts = 5000;
NbParts = 200;

BetasRecord = zeros(NbIts,sum(Data.NbComputingSteps)+1);
SigmasRecord = Parameters.Sigma.RealValue;
Parameters.Sigma.Value = Parameters.Sigma.RealValue;
Accepted = [];
lambda = 0.0251;
C = 0.99;
% beta given sigma
BetasSamples = zeros(NbParts,sum(Data.NbComputingSteps));
ind = 1;
for IndStep = 2:Parameters.ObsLength
    for IndCStep = 1:1/Parameters.ComputationTStep
        BetasSamples(:,ind+1) = BetasSamples(:,ind)+randn(NbParts,1)*Parameters.Sigma.RealValue*sqrt(Parameters.ComputationTStep);
        ind = ind+1;
    end
    Weigths = normpdf(BetasSamples(:,ind),Data.Observations(IndStep-1),Parameters.SigmaObs);
    Weigths = Weigths/sum(Weigths);
    u = rand(1,1)/NbParts;
    s = 0;
    KeptInds = [];
    resind = 1;
    for ipart = 1:NbParts
        k = 0;
        s = s+Weigths(ipart);
        while s>u
            k=k+1;
            u = u+1/NbParts;
            KeptInds(resind) = ipart;
            resind = resind+1;
        end
    end
    BetasSamples = BetasSamples(KeptInds,:);
end
for IndIt = 2:NbIts
    % Sigma given Beta
    RandInd = ceil(rand(1,1)*NbParts);
    Beta = BetasSamples(RandInd,:);
    BetasRecord(IndIt,:) = Beta;
    Parameters.SigmaStar.Value = Parameters.Sigma.Value + lambda*randn(1,1);
    BetaStar = Beta/Parameters.Sigma.Value*Parameters.SigmaStar.Value;
    LogLik = sum(log(normpdf(Beta(cumsum(ones(1,Parameters.ObsLength-1)/Parameters.ComputationTStep)),Data.Observations,Parameters.SigmaObs)));
    LogLikStar = sum(log(normpdf(BetaStar(cumsum(ones(1,Parameters.ObsLength-1)/Parameters.ComputationTStep)),Data.Observations,Parameters.SigmaObs)));
    AccRate = LogLikStar-LogLik;
    if log(rand(1,1))<AccRate
        Parameters.Sigma.Value = Parameters.SigmaStar.Value;
        Beta = BetaStar;
        Accepted(IndIt-1) = 1;
    else
        Accepted(IndIt-1) = 0;
    end
    if IndIt>10
        lambda = exp(log(lambda)-C^IndIt*(0.23-mean(Accepted)));
    end
    SigmasRecord(IndIt) = Parameters.Sigma.Value;
    
    % beta given sigma
    BetasSamples = zeros(NbParts,sum(Data.NbComputingSteps));
    ind = 1;
    for IndStep = 2:Parameters.ObsLength
        for IndCStep = 1:1/Parameters.ComputationTStep
            BetasSamples(:,ind+1) = BetasSamples(:,ind)+randn(NbParts,1)*Parameters.Sigma.RealValue*sqrt(Parameters.ComputationTStep);
            BetasSamples(1,ind+1) = Beta(ind+1);
            ind = ind+1;
        end
        Weigths = normpdf(BetasSamples(:,ind),Data.Observations(IndStep-1),Parameters.SigmaObs);
        Weigths = Weigths/sum(Weigths);
        u = rand(1,1)/NbParts;
        s = 0;
        KeptInds = [];
        resind = 1;
        for ipart = 1:NbParts
            k = 0;
            s = s+Weigths(ipart);
            while s>u
                k=k+1;
                u = u+1/NbParts;
                KeptInds(resind) = ipart;
                resind = resind+1;
            end
        end
        KeptInds(1) = 1;  
        BetasSamples = BetasSamples(KeptInds,:);
    end
    disp([num2str(IndIt) '  ' num2str(mean(Accepted)*100,3) '%   ' num2str(lambda,3) ])
end


subplot(2,1,1)
plot(SigmasRecord)
subplot(2,1,2)

plot(Data.Beta)
hold on
plot(cumsum(ones(1,Parameters.ObsLength-1)/Parameters.ComputationTStep),Data.Observations,'g')
plot(mean(BetasRecord),'k')
plot(quantile(BetasRecord,0.025),'r')
plot(quantile(BetasRecord,0.975),'r')
hold off



