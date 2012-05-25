cd('\\STUDENT1\D_USERS\DUREAU\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Filtering')
addpath('\\STUDENT1\D_USERS\DUREAU\PhD Work\Matlab Scripts\General Tools')
addpath('\\STUDENT1\D_USERS\DUREAU\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Parameter Estimation')
addpath('\\STUDENT1\D_USERS\DUREAU\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Filtering')
addpath('\\STUDENT1\D_USERS\DUREAU\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Joint Sampling')

InitialState.S = 10000;
InitialState.I = 100;
InitialState.R = 0;

Parameters = struct();
Parameters.Beta = 0.01;
Parameters.Gamma = 0.1;
Parameters.SigmaObs = 0.01;
Parameters.ObsNoise = 0;
Parameters.SigmaDiffusion(1) = 0.0001;
Parameters.SigmaDiffusion(2) = 0.0001;
Parameters.SigmaDiffusion(3) = 0.0001;
Parameters.SigmaDiffusion(4) = 0.2/InitialState.S;
Parameters.BetaInf = 0.9*Parameters.Beta;
Parameters.BetaSup = 1.1*Parameters.Beta;
Parameters.GammaInf = 0.999999999999*Parameters.Gamma;
Parameters.GammaSup = 1.000000000001*Parameters.Gamma;
Parameters.ObservationTStep = 0.1;
Parameters.ComputationTStep = 0.001;
Parameters.ObservationLength = 15;

time = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
RealLambdas =1*(2-Sigmoid(time-7))/InitialState.S;
RealTrajectory = ComputeTrajectory(InitialState,RealLambdas,Parameters);

Data.RealData = [ RealTrajectory.Ss ; RealTrajectory.Is; RealTrajectory.Rs; RealTrajectory.LambdasTime];
Data.Observations = RealTrajectory.Os;
Data.Instants = RealTrajectory. ObsTime;


Parameters.NbParticules = 200;
Results = {};

NbIters = 30
for i = 1:NbIters
    i
    x0(1) = rand(1,1)*(Parameters.BetaSup - Parameters.BetaInf) + Parameters.BetaInf ; 
    x0(2) = rand(1,1)*(Parameters.GammaSup - Parameters.GammaInf) + Parameters.GammaInf ; 
    Results{i} = OptParEstimation(@EstimationSMC, Data, Parameters, x0);
end

plot(Results{1}.HistoryX(:,1),Results{1}.HistoryFval,'.')
xlim([Parameters.BetaInf Parameters.BetaSup])
hold on
plot(ResultsComp.SamplesTheta(:,1),ResultsComp.LogLikelihoods,'.k')
for i = 1:NbIters
    plot(Results{i}.HistoryX(1,1),Results{i}.HistoryFval(1),'or')
    plot(Results{i}.HistoryX(end,1),Results{i}.HistoryFval(end),'og')
    plot(Results{i}.HistoryX(:,1),Results{i}.HistoryFval,'--')
end
hold off


%%%%%%%%%%%%%%%%%%%%%%%
% Likelihood cuts

% Along beta
noise = 0.3;
v = -1:0.2:1;
values = Parameters.Beta*(1 + noise*v);
logliks = [];
for i = 1:length(values)
    i
    x = [values(i) Parameters.Gamma];
    logliks(i) = ParEst_EstimationSMC(x, Data, Parameters, 100);
end
plot(values,logliks)


% Along gamma
noise = 0.3;
v = -1:0.2:1;
values = Parameters.Gamma*(1 + noise*v);
logliks = [];
logliks2 = [];
for i = 1:length(values)
    i
    x = [Parameters.Beta values(i)];
    logliks(i) = ParEst_EstimationSMC(x, Data, Parameters, 100);
    logliks2(i) = ParEst_EstimationSMCSarkka(x, Data, Parameters, 100);
end
plot(values,logliks)
hold on
plot(values,logliks,'k')
hold off