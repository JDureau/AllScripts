function ResGibbs = SimpleGibbs(Betas,Observations,Parameters,NbIterations)


n =  length(Observations);
Data.Beta = Betas;
Data.Observations = Observations;


NbIts = 2000;
NbParts = 500;

BetasRecord = zeros(NbIts,n);
SigmasRecord = Parameters.SigmaRW.Value;
Parameters.Sigma.Value = Parameters.SigmaRW.Value;

Accepted = [];
lambda = Parameters.Lambda;
C = 0.995;
% beta given sigma
BetasSamples = zeros(NbParts,n);
for IndStep = 2:n
    BetasSamples(:,IndStep) = BetasSamples(:,IndStep-1)+randn(NbParts,1)*Parameters.Sigma.Value;
    Weigths = normpdf(BetasSamples(:,IndStep),Data.Observations(IndStep),Parameters.SigmaObs);
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
plot(Data.Beta)
hold on
plot(Data.Observations,'g')
RandInd = ceil(rand(1,1)*NbParts);
Beta = BetasSamples(RandInd,:);
plot(Beta,'r')
hold off
alpha = 0.0001;
beta = 0.0001;
Accepted = [];
C = 0.99;


% TempPar.SigmaRW.TransfValue = log(0.1);
% Parameters.SigmaRW.Value = 0.1;
% Parameters = UpdateParsNoTransfToTransf(Parameters);
% TempPar.Paths = [Beta;Beta];
% ResGibbs = RunEstimationMethod(DataGen, DObs,Parameters,TempPar,1000);


for IndIt = 2:NbIterations
    disp(IndIt)
    % Sigma given Beta
    RandInd = ceil(rand(1,1)*NbParts);
    Beta = BetasSamples(RandInd,:);
    BetasRecord(IndIt,:) = Beta;
    
%     Parameters.Sigma.Value = sqrt(gamrnd((n-1)/2+alpha,(beta+sum(diff(Beta).^2)/2)^-1)^-1);
    Parameters.SigmaStar.Value = Parameters.Sigma.Value + randn(1,1)*lambda;
    LogLikStar = sum(log(normpdf((Beta/Parameters.Sigma.Value)*Parameters.SigmaStar.Value,Data.Observations,Parameters.SigmaObs)));
    LogLik = sum(log(normpdf(Beta,Data.Observations,Parameters.SigmaObs)));
    AccRate = LogLikStar-LogLik;
    if log(rand(1,1))<AccRate
        Beta = (Beta/Parameters.Sigma.Value)*Parameters.SigmaStar.Value;
        Parameters.Sigma.Value = Parameters.SigmaStar.Value;
        Accepted(end+1) = 1;
    else
        Accepted(end+1) = 0;
    end
%     if IndIt>10
%         lambda = exp(log(lambda)-C^IndIt*(0.23-mean(Accepted)));
%     end
    [mean(Accepted)]
    
    SigmasRecord(IndIt) = Parameters.Sigma.Value;
    
    % beta given sigma
    BetasSamples = zeros(NbParts,n);
    for IndStep = 2:n
        BetasSamples(:,IndStep) = BetasSamples(:,IndStep-1)+randn(NbParts,1)*Parameters.Sigma.Value;
        BetasSamples(1,IndStep) = Beta(IndStep);
        Weigths = normpdf(BetasSamples(:,IndStep),Data.Observations(IndStep),Parameters.SigmaObs);
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
end
mean(SigmasRecord.^2)
var(SigmasRecord.^2)

subplot(3,1,1)
plot(SigmasRecord)
subplot(3,1,2)
plot(Data.Beta)
hold on
plot(Data.Observations,'g')
plot(mean(BetasRecord),'k')
plot(quantile(BetasRecord,0.025),'r')
plot(quantile(BetasRecord,0.975),'r')
hold off

subplot(3,1,3)
hist(SigmasRecord)

ResGibbs.SigmasRecord = SigmasRecord;
ResGibbs.Parameters = Parameters;



