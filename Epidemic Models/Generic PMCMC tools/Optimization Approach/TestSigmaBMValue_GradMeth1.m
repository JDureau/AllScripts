function [f,g,h] = TestSigmaBMValue_GradMeth1(sigma, Parameters, Data, NbResampling)

% This function applies a filter and returns Loglikelihood and gradient

% sigma = exp(logsigma);

ObservedVariables = Parameters.ObservedVariables;
beta = Parameters.Beta;
gamma = Parameters.Gamma;

NbVarsTot = size(Data.RealData,1);
Observations = Data.Observations;
ObservationInstants = Data.Instants;

% Bootstrap Filter adapted to SIR problem to be faster
ObservationTStep = mean(diff(ObservationInstants));
ComputationTStep = Parameters.ComputationTStep;


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

Grads  = zeros(NbResampling,1);
Hesss  = zeros(NbResampling,1);


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
    
    Grads  = Grads(NbResampling,1);
    Hesss  = Hesss(NbResampling,1);

    
    
    NbTSteps = ceil(ObservationTStep/ComputationTStep);
       
    rands = randn(NbResampling,NbTSteps);
    tempNewPotentialStatesCellVect = NewPotentialStatesCellVect;
    for IndDiscr = 1:NbTSteps
        tempNewPotentialStatesCellVect(:,1) = tempNewPotentialStatesCellVect(:,1) + (-NewPotentialStatesCellVect(:,4).*NewPotentialStatesCellVect(:,1).*NewPotentialStatesCellVect(:,2)+beta*NewPotentialStatesCellVect(:,3))*ComputationTStep;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands1(:,IndDiscr);
        tempNewPotentialStatesCellVect(:,2) = tempNewPotentialStatesCellVect(:,2) + (NewPotentialStatesCellVect(:,4).*NewPotentialStatesCellVect(:,1).*NewPotentialStatesCellVect(:,2)-gamma*NewPotentialStatesCellVect(:,2))*ComputationTStep;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands2(:,IndDiscr);
        tempNewPotentialStatesCellVect(:,3) = tempNewPotentialStatesCellVect(:,3) + (-beta*NewPotentialStatesCellVect(:,3) + gamma*NewPotentialStatesCellVect(:,2))*ComputationTStep;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands3(:,IndDiscr);
              
        tempNewPotentialStatesCellVect(:,4) = tempNewPotentialStatesCellVect(:,4) + sqrt(ComputationTStep)*sigma*rands(:,IndDiscr);
        
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

    
    
    CumWeigths=cumsum(Weigths)/sum(Weigths);
    NbKeptSamples(IndTime) = length(unique(KeptSamples));
    Grads = Grads +  (PotentialStatesCellVect(:,4) - ResampledParticules(IndTime-1,:,4)').^2/(sigma^3*ObservationTStep) - 1/sigma;
    Hesss  = Hesss  +  3*(PotentialStatesCellVect(:,4) - ResampledParticules(IndTime-1,:,4)').^2/(sigma^4*ObservationTStep) + 1/sigma^2;

end

Grad = sum(Weigths.*Grads)/sum(Weigths);
Hess = sum(Weigths.*Hesss)/sum(Weigths);


f = LogLik;
if nargout == 2
    g = Grad;
end
if nargout == 3
    g = Grad;
    h = Hess;
end
disp(['Sigma = ' num2str(sigma)])
disp(['LogLik = ' num2str(-f)])
