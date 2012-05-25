function Result = MIF(Data,Parameters)

Result = struct();
Result.SamplesX = [];
Result.SamplesTheta = [];

ObservedVariables = Parameters.ObservedVariables;
Observations = Data.Observations;
ObservationInstants = Data.Instants;
ComputationTStep = Parameters.ComputationTStep;
ObservationTStep = mean(diff(ObservationInstants));
NbTSteps = ceil(ObservationTStep/ComputationTStep);
NbIterations = Parameters.MIFNbIterations;
NbParticules = Parameters.MIFNbParticules;
a = Parameters.MIFCoolingParameters;
b = Parameters.MIFb;

ThetaRecord = [];

% Initialize Theta
Thetas = ones(NbParticules,2);
Thetas(:,1) = Parameters.Beta*Thetas(:,1);
Thetas(:,2) = Parameters.Gamma*Thetas(:,2);
ThetasRecord(:,1) = mean(Thetas);
ThetasRecordTemp = zeros(length(ObservationInstants),2);
ThetasVarRecordTemp = zeros(length(ObservationInstants),2);
LogLiks = [];

for IndIt = 2:NbIterations

    disp(['Iteration ' num2str(IndIt)])
    Thetas(:,1) = ThetasRecord(1,IndIt-1);
    Thetas(:,2) = ThetasRecord(2,IndIt-1);
    for i = 1:NbParticules
        Xf(i,:) = Data.RealData(:,1);
    end
    
    Thetas(:,1) = Thetas(:,1) + a^(IndIt-1)*b*randn(NbParticules,1)*Parameters.MIFSigmaBeta;
    Thetas(:,2) = Thetas(:,2) + a^(IndIt-1)*b*randn(NbParticules,1)*Parameters.MIFSigmaGamma;
    ThetasRecordTemp(1,:) = mean(Thetas) ; 
    ThetasVarRecordTemp(1,:) = var(Thetas) ; 
    RecordStates = zeros(length(ObservationInstants),NbParticules,4);
    RecordThetas = zeros(length(ObservationInstants),NbParticules,2);
    LogLiksTemp = zeros(1,length(ObservationInstants));
    
    for IndTime = 2:length(ObservationInstants)
        
        TempXf = Xf;
        rands = randn(NbParticules,NbTSteps);
        for IndDiscr = 1:NbTSteps
            TempXf(:,1) =  TempXf(:,1) + (- Xf(:,1).*Xf(:,2).*Xf(:,4) + Thetas(:,1).*Xf(:,3))*ComputationTStep;
            TempXf(:,2) =  TempXf(:,2) + (  Xf(:,1).*Xf(:,2).*Xf(:,4) - Thetas(:,2).*Xf(:,2))*ComputationTStep;
            TempXf(:,3) =  TempXf(:,3) + (- Thetas(:,1).*Xf(:,3) + Thetas(:,2).*Xf(:,2))*ComputationTStep;
            TempXf(:,4) =  TempXf(:,4) + ( sqrt(ComputationTStep)*Parameters.SigmaDiffusion(4)*rands(:,IndDiscr));
            Xf = TempXf;
        end
        
        Weigths = ones(NbParticules,1);
        for IndObservedVar = 1:length(Parameters.ObservedVariables)    
            Weigths = Weigths.*normpdf(Xf(:,ObservedVariables(IndObservedVar)),Observations(ObservedVariables(IndObservedVar), IndTime),Parameters.SigmaObs(ObservedVariables(IndObservedVar))*Observations(ObservedVariables(IndObservedVar),IndTime));
        end
        ResampledIndexes = Resample(Weigths);
        if length(unique(ResampledIndexes))<NbParticules/10
            disp(['Nb Kept Particules : ' num2str(length(unique(ResampledIndexes)))]) 
        end
        Xf(:,1) = Xf(ResampledIndexes,1);
        Xf(:,2) = Xf(ResampledIndexes,2);
        Xf(:,3) = Xf(ResampledIndexes,3);
        Xf(:,4) = Xf(ResampledIndexes,4);
        ThetasRecordTemp(IndTime,:) = mean(Thetas(ResampledIndexes,:))';
        ThetasVarRecordTemp(IndTime,:) = var(Thetas(ResampledIndexes,:));
        Thetas(:,1) = Thetas(ResampledIndexes,1) + randn(NbParticules,1)*a^(IndIt-1)*ObservationTStep*Parameters.MIFSigmaBeta;
        Thetas(:,2) = Thetas(ResampledIndexes,2) + randn(NbParticules,1)*a^(IndIt-1)*ObservationTStep*Parameters.MIFSigmaGamma;
        RecordStates(IndTime,:,:) = Xf;
        RecordThetas(IndTime,:,:) = Thetas;
        LogLiksTemp(IndTime) = log(sum(Weigths)/NbParticules);
    end
    temp = 0;
    for IndTime = 2:length(ObservationInstants)
        temp = temp + ThetasVarRecordTemp(2,1).*(ThetasVarRecordTemp(IndTime,1).^-1).*(ThetasRecordTemp(IndTime,1)-ThetasRecordTemp(IndTime-1,1));
    end
    ThetasRecord(1,IndIt) = ThetasRecord(1,IndIt-1) + temp';
    ThetasRecord(2,IndIt) = ThetasRecord(2,IndIt-1);
    LogLiks(IndIt-1) = sum(LogLiksTemp);
end

Result.ThetasRecord = ThetasRecord;
Result.RecordStates = RecordStates;
Result.LogLiks = LogLiks;





