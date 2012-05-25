function ResultSMC = EstimationSMC(Pars, Data, Parameters, NbResampling, Task)

% If Task = Filter, we return ResultSMC with all filtering informaiton.
% If Task = ParEst, we just return the loglikelihood
if nargin == 4
    Task = 'Filter';
end

ObservedVariables = Parameters.ObservedVariables;
beta = Pars(1);
gamma = Pars(2);

tic 
NbVarsTot = size(Data.RealData,1);
Observations = Data.Observations;
ObservationInstants = Data.Instants;

% Bootstrap Filter adapted to SIR problem to be faster


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



PotentialStatesCellVect = ones(NbResampling,NbVarsTot);
NewPotentialStatesCellVect= [];

InitialState = Data.RealData(:,1);

PreviousState = InitialState;
PosteriorMeansRecord = diag(InitialState)*ones(NbVarsTot,length(ObservationInstants));
NbKeptSamples = NbResampling*ones(1,length(ObservationInstants));
Posterior975Record = InitialState; 
Posterior225Record = InitialState;

SamplesMat = PreviousState(end)*ones(NbResampling,1);
ResampledParticules = [];
ResampledPredictionParticules = [];
WeigthsMat = 1/NbResampling*ones(NbResampling,1);

for IndResampling = 1:NbResampling
    PotentialStatesCellVect(IndResampling,:)=PreviousState;
end
CumWeigths=(1:NbResampling)/NbResampling;

Liks = 1;
LogLik = 0;

for IndTime = 2:length(ObservationInstants)
%     disp(IndTime)
 
    % Variables initialisation
    Weigths = zeros(1,NbResampling);
    Variables = zeros(NbVarsTot,NbResampling);
    
    RandU1 = rand(1,1)/NbResampling;
    KeptSamples = zeros(1,NbResampling);
    
%      plot(1,Data.RealData(2,1),'or')
%     hold on
%     plot(10,Data.RealData(2,2),'or')
%     tmp = [];
    
    
    ResampledPredictionParticules(IndTime-1,:,:) = PotentialStatesCellVect;
    Current = 1;
    for IndResampling = 1:NbResampling
        % Resample.
%         Ind = find(CumWeigths>RandU1+(IndResampling-1)/NbResampling,1);
%         KeptSamples(IndResampling) = Ind;
        while CumWeigths(Current)<=RandU1+(IndResampling-1)/NbResampling
            Current = Current+1;
        end
        KeptSamples(IndResampling) = Current;
    end
    NewPotentialStatesCellVect = PotentialStatesCellVect(KeptSamples,:);
    ResampledParticules(IndTime-1,:,:) = PotentialStatesCellVect(KeptSamples,:);
    
    NbTSteps = ceil(ObservationTStep/ComputationTStep);
       
    rands = randn(NbResampling,NbTSteps);
    tempNewPotentialStatesCellVect = NewPotentialStatesCellVect;
    for IndDiscr = 1:NbTSteps
        tempNewPotentialStatesCellVect(:,1) = tempNewPotentialStatesCellVect(:,1) + (-NewPotentialStatesCellVect(:,4).*NewPotentialStatesCellVect(:,1).*NewPotentialStatesCellVect(:,2)+beta*NewPotentialStatesCellVect(:,3))*ComputationTStep;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands1(:,IndDiscr);
        tempNewPotentialStatesCellVect(:,2) = tempNewPotentialStatesCellVect(:,2) + (NewPotentialStatesCellVect(:,4).*NewPotentialStatesCellVect(:,1).*NewPotentialStatesCellVect(:,2)-gamma*NewPotentialStatesCellVect(:,2))*ComputationTStep;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands2(:,IndDiscr);
        tempNewPotentialStatesCellVect(:,3) = tempNewPotentialStatesCellVect(:,3) + (-beta*NewPotentialStatesCellVect(:,3) + gamma*NewPotentialStatesCellVect(:,2))*ComputationTStep;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands3(:,IndDiscr);
        if Parameters.DiffusionType =='GBM'
            tempNewPotentialStatesCellVect(:,4) = tempNewPotentialStatesCellVect(:,4) + sqrt(ComputationTStep)*Parameters.SigmaGBM*tempNewPotentialStatesCellVect(:,4)*rands(:,IndDiscr);
        elseif Parameters.DiffusionType =='OUD'
            tempNewPotentialStatesCellVect(:,4) = tempNewPotentialStatesCellVect(:,4) + (Parameters.MuOU - tempNewPotentialStatesCellVect(:,4))*Parameters.KappaOU*ComputationTStep + sqrt(ComputationTStep)*Parameters.SigmaOU*rands(:,IndDiscr);            
        elseif Parameters.DiffusionType =='OUP'
            tempNewPotentialStatesCellVect(:,4) = tempNewPotentialStatesCellVect(:,4) + (((Parameters.MuOU-tempNewPotentialStatesCellVect(:,4))/Parameters.SigmaPriorOUP).^Parameters.PowerOUP)*ComputationTStep + sqrt(ComputationTStep)*Parameters.SigmaBMOUP*rands(:,IndDiscr);            
        else             
            tempNewPotentialStatesCellVect(:,4) = tempNewPotentialStatesCellVect(:,4) + sqrt(ComputationTStep)*Parameters.SigmaRW*rands(:,IndDiscr);
        end
        NewPotentialStatesCellVect = tempNewPotentialStatesCellVect;
    end        
     
    Weigths = ones(NbResampling,1);
    for IndObservedVar = 1:length(Parameters.ObservedVariables)    
        Weigths = Weigths.*normpdf(NewPotentialStatesCellVect(:,ObservedVariables(IndObservedVar)),Observations(ObservedVariables(IndObservedVar), IndTime),Parameters.SigmaObs(ObservedVariables(IndObservedVar))*Observations(ObservedVariables(IndObservedVar),IndTime));
    end
    
    Liks(IndTime) = (mean(Weigths));
    LogLik = LogLik + log(mean(Weigths));
    
    Variables = NewPotentialStatesCellVect';
     
    Weigths = Weigths/sum(Weigths);
    WeigthsMat(:,IndTime) = Weigths;

    PotentialStatesCellVect = NewPotentialStatesCellVect;

    
    
    
    
%     
    
    CumWeigths=cumsum(Weigths)/sum(Weigths);
    NbKeptSamples(IndTime) = length(unique(KeptSamples));

end




PosteriorMeansRecord = zeros(4,length(ObservationInstants));
if Task == 'Filter'
    %Resample last step:
    Current = 1;
    for IndResampling = 1:NbResampling
        % Resample.
    %         Ind = find(CumWeigths>RandU1+(IndResampling-1)/NbResampling,1);
    %         KeptSamples(IndResampling) = Ind;
        while CumWeigths(Current)<=RandU1+(IndResampling-1)/NbResampling
            Current = Current+1;
        end
        KeptSamples(IndResampling) = Current;
    end
    ResampledParticules(length(ObservationInstants),:,:) = PotentialStatesCellVect(KeptSamples,:);
    ResampledPredictionParticules(length(ObservationInstants),:,:) = PotentialStatesCellVect;
    for IndTime = 1:length(ObservationInstants)
        for var =1:size(Variables,1)
            q = quantile(ResampledParticules(IndTime,:,var),[0.025,0.975]);
            Posterior975Record(var,IndTime) = q(1);
            Posterior225Record(var,IndTime) = q(2);
            PosteriorMeansRecord(var,IndTime) = mean(ResampledParticules(IndTime,:,var));
            PredictionMeansRecord(var,IndTime) = mean(ResampledPredictionParticules(IndTime,:,var));
            q = quantile(ResampledPredictionParticules(IndTime,:,var),[0.025,0.975]);
            Prediction975Record(var,IndTime) = q(1);
            Prediction225Record(var,IndTime) = q(2);
        end
    end
end


if Task == 'Filter'
    ResultSMC.PosteriorMeansRecord = PosteriorMeansRecord;
    ResultSMC.PredictionMeansRecord = PredictionMeansRecord;
    ResultSMC.NbKeptSamples = NbKeptSamples;
    ResultSMC.NbParticules = NbResampling;
    ResultSMC.Data = Data;
    ResultSMC.Parameters = Parameters;
    ResultSMC.Posterior975Record = Posterior975Record;
    ResultSMC.Posterior225Record = Posterior225Record;
    ResultSMC.Prediction975Record = Prediction975Record;
    ResultSMC.Prediction225Record = Prediction225Record;
    ResultSMC.ResampledParticules = ResampledParticules;
    ResultSMC.WeigthsRecord = WeigthsMat;
    ResultSMC.Likelihood = prod(Liks);
    ResultSMC.Liks = Liks;
    ResultSMC.LogLik = LogLik;
    ResultSMC.RealTime = toc;
else
    ResultSMC = LogLik;
end
