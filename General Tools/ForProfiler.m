function [] = ForProfiler()

% Step 1: generate data
InitialParameters = struct();
InitialParameters.Beta = 1/100;
InitialParameters.Gamma = 1/10;
InitialState.S = 10000;
InitialState.I = 100;
InitialState.R = 0;
SigmaObs = 0;
SigmaLambda = 0.1/InitialState.S;

ObservationTimeStep = 1;
ComputationTimeStep = 0.01;
ObservationLength = 100;

time = 0:ComputationTimeStep:ObservationLength;
RealLambdas = (2-Sigmoid(time-50))/InitialState.S;
% RealLambdas = (2*ones(size(time)))/InitialState.S;

PreviousParameters = InitialParameters;
PreviousState = InitialState;
Stab = InitialState.S;
Otab = InitialState.I+ randn(1,1)*SigmaObs*InitialState.I;
Itab = InitialState.I;
Rtab = InitialState.R;
Time = 0;
RealTrajectory = struct();
for i = ObservationTimeStep:ObservationTimeStep:ObservationLength
    inds = find(and(time>=(i-ObservationTimeStep),time<i));
    ExtractedLambdas = RealLambdas(inds);
    [PreviousState,PreviousParameters] = ComputeNextStateGivenLambda(PreviousState,PreviousParameters,ExtractedLambdas,ObservationTimeStep,ComputationTimeStep);
    Stab(end + 1) = PreviousState.S;
    Itab(end + 1) = PreviousState.I;
    Otab(end + 1) = PreviousState.I + randn(1,1)*SigmaObs*PreviousState.I;
    Rtab(end + 1) = PreviousState.R;
    Time(end+1) = i;
end
RealTrajectory.Ss = Stab;
RealTrajectory.Is = Itab;
RealTrajectory.Os = Otab;
RealTrajectory.Rs = Rtab;
RealTrajectory.Lambdas = RealLambdas;
RealTrajectory.Parameters = InitialParameters;
subplot(2,1,1)
plot(Time,Itab)
xlabel('I')
subplot(2,1,2)
plot(time,RealLambdas)
xlabel('\lambda')
Data.Measures = Itab;
Data.Instants = Time;

% Step 2: estimate p(Theta_i|Y)
NbIters = 100;
NbBlocks = 10;
BlockLength=Time(end)/NbBlocks;
CurrentTrajectory = RealTrajectory;
LambdaRecord = [];
IRecord = [];
AcceptedParticles = 0;
TotalParticles = 0;
for IndIter = 1:NbIters
    disp(['Iter ' num2str(IndIter) ])
    disp(['Acceptance rate : ' num2str(AcceptedParticles/TotalParticles)])
    for IndBlock = 1:NbBlocks-1
        tk = (IndBlock-1)*BlockLength*ObservationTimeStep;
        tkp2 = min((IndBlock+1)*BlockLength*ObservationTimeStep,Time(end));
        Temp_time = tk:ComputationTimeStep:tkp2;
        rands = randn(1,length(Temp_time))*sqrt(ComputationTimeStep)*SigmaLambda;
        LambdasInd = (IndBlock-1)*BlockLength*ObservationTimeStep/ComputationTimeStep + 1 ;
        StartingLambda = CurrentTrajectory.Lambdas(LambdasInd);
        AimedLambda = CurrentTrajectory.Lambdas(min(LambdasInd+2*BlockLength*ObservationTimeStep/ComputationTimeStep,length(time)));
        Alphas = (Temp_time-tk)/(tkp2-tk);
        PotLambdas_ModifiedBlock = StartingLambda + cumsum(rands) + Alphas.*(AimedLambda-StartingLambda-sum(rands));
        PotLambdas = CurrentTrajectory.Lambdas;
        PotLambdas(LambdasInd:LambdasInd+length(PotLambdas_ModifiedBlock)-1) = PotLambdas_ModifiedBlock;
        PotTrajectory = ModifyTrajectory(CurrentTrajectory,PotLambdas,time,ObservationTimeStep,ComputationTimeStep);
        LogLikPot = 0;
        LogLikCur = 0;
        for IndTStep = (IndBlock-1)*BlockLength + 1 : length(Time)    
            LogLikPot = LogLikPot + log(normpdf(PotTrajectory.Is(IndTStep),Itab(IndTStep),max(SigmaObs,0.1)*Itab(IndTStep)));
            LogLikCur = LogLikCur + log(normpdf(CurrentTrajectory.Is(IndTStep),Itab(IndTStep),max(SigmaObs,0.1)*Itab(IndTStep)));
        end
        if LogLikPot>LogLikCur
            CurrentTrajectory = PotTrajectory;
            AcceptedParticles = AcceptedParticles + 1;
        else if log(rand(1,1))<LogLikPot-LogLikCur
                CurrentTrajectory = PotTrajectory;
                AcceptedParticles = AcceptedParticles + 1;
            end
        end
        TotalParticles = TotalParticles + 1;
        Inds = 1:ObservationTimeStep/ComputationTimeStep:length(CurrentTrajectory.Lambdas);
        LambdaRecord(end+1,:) = CurrentTrajectory.Lambdas(Inds);
        IRecord(end+1,:) = CurrentTrajectory.Is; 
    end
end

