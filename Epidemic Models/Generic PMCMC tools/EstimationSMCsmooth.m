function ResultSMC = EstimationSMCsmooth(Pars, Data, Parameters, NbResampling, Task, ForcedTrajectory)

% we took of noise on SIR transitions to male code faster. (end of august)
% If Task = Filter, we return ResultSMC with all filtering informaiton.
% If Task = ParEst, we just return the loglikelihood
if nargin == 4
    Task = 'Filter';
end

if nargin == 6
    ForcedTraj = 1;
else
    ForcedTraj = 0;
end

beta = Pars(1);
gamma = Pars(2);

tic
NbVarsTot = size(Data.RealData,1);
Observations = Data.Observations;
ObservationInstants = Data.Instants;
ObservedVariables = Parameters.ObservedVariables;

% Verifications
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
PotentialStatesCellVect = ones(NbResampling,NbVarsTot);
InitialState = Data.RealData(:,1);
PreviousState = InitialState;
PosteriorMeansRecord = diag(InitialState)*ones(NbVarsTot,length(ObservationInstants));
NbKeptSamples = NbResampling*ones(1,length(ObservationInstants));
LikelihoodsNew = ones(1,NbResampling);
ProportionKept = [];

ResampledParticules = zeros(length(ObservationInstants),4,NbResampling);
for IndTime = 1:length(ObservationInstants)
    ResampledParticules(IndTime,:,:) = diag(PreviousState)*ones(length(PreviousState),NbResampling);
end

CompleteLambdas = zeros(NbResampling,NbTSteps*(length(ObservationInstants)-1));
% CompleteLambdas2 = zeros(NbResampling,NbTSteps*(length(ObservationInstants)-1));

for IndResampling = 1:NbResampling
    PotentialStatesCellVect(IndResampling,:)=PreviousState;
    CompleteLambdas(IndResampling,1) = PreviousState(4);
    CompleteLambdas2(IndResampling,1) = PreviousState(4);
end
CumWeigths=(1:NbResampling)/NbResampling;


Liks = ones(1,length(ObservationInstants));
LogLik = 0;
FathersTab = [];

for IndTime = 2:length(ObservationInstants)
    
    RandU1 = rand(1,1)/NbResampling;
    KeptSamples = zeros(1,NbResampling);
    
    Current = 1;
    for IndResampling = 1:NbResampling
        while CumWeigths(Current)<=RandU1+(IndResampling-1)/NbResampling
            Current = Current+1;
        end
        KeptSamples(IndResampling) = Current;
    end
    if ForcedTraj
        KeptSamples(1) = 1;
    end
%     length(find(KeptSamples == 1))
    
    NewPotentialStatesCellVect = PotentialStatesCellVect(KeptSamples,:);
    ResampledParticules = ResampledParticules(:,:,KeptSamples);
%     CompleteLambdas2 = CompleteLambdas2(KeptSamples,:);
    LikelihoodsNew = LikelihoodsNew(:,KeptSamples);
    FathersTab(:,IndTime-1) = KeptSamples;
       
    rands = randn(NbResampling,NbTSteps);
    tempNewPotentialStatesCellVect = NewPotentialStatesCellVect ;
    for IndDiscr = 1:NbTSteps
        tempNewPotentialStatesCellVect(:,1) = tempNewPotentialStatesCellVect(:,1) + (-NewPotentialStatesCellVect(:,4).*NewPotentialStatesCellVect(:,1).*NewPotentialStatesCellVect(:,2)+beta*NewPotentialStatesCellVect(:,3))*ComputationTStep ;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
        tempNewPotentialStatesCellVect(:,2) = tempNewPotentialStatesCellVect(:,2) + ( NewPotentialStatesCellVect(:,4).*NewPotentialStatesCellVect(:,1).*NewPotentialStatesCellVect(:,2)-gamma*NewPotentialStatesCellVect(:,2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands(:,IndDiscr,2);
        tempNewPotentialStatesCellVect(:,3) = tempNewPotentialStatesCellVect(:,3) + (-beta*NewPotentialStatesCellVect(:,3) + gamma*NewPotentialStatesCellVect(:,2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
        if Parameters.DiffusionType =='GBM'
            tempNewPotentialStatesCellVect(:,4) = tempNewPotentialStatesCellVect(:,4) + sqrt(ComputationTStep)*Parameters.SigmaGBM*tempNewPotentialStatesCellVect(:,4).*rands(:,IndDiscr);
        elseif Parameters.DiffusionType =='OUD'
            tempNewPotentialStatesCellVect(:,4) = tempNewPotentialStatesCellVect(:,4) + (Parameters.MuOU-tempNewPotentialStatesCellVect(:,4))*Parameters.KappaOU*ComputationTStep + sqrt(ComputationTStep)*Parameters.SigmaOU*rands(:,IndDiscr);            
        elseif Parameters.DiffusionType =='OUP'
            tempNewPotentialStatesCellVect(:,4) = tempNewPotentialStatesCellVect(:,4) + (((Parameters.MuOU-tempNewPotentialStatesCellVect(:,4))/Parameters.SigmaPriorOUP).^Parameters.PowerOUP)*ComputationTStep + sqrt(ComputationTStep)*Parameters.SigmaBMOUP*rands(:,IndDiscr);            
        else             
            tempNewPotentialStatesCellVect(:,4) = tempNewPotentialStatesCellVect(:,4) + sqrt(ComputationTStep)*Parameters.SigmaRW*rands(:,IndDiscr);
        end
        if ForcedTraj
            tempNewPotentialStatesCellVect(1,4) = ForcedTrajectory((IndTime-2)*NbTSteps + IndDiscr);
        end            
        NewPotentialStatesCellVect = tempNewPotentialStatesCellVect;
        CompleteLambdas(:,(IndTime-2)*NbTSteps + IndDiscr) = tempNewPotentialStatesCellVect(:,4);
%         CompleteLambdas2(:,(IndTime-2)*NbTSteps + IndDiscr) = tempNewPotentialStatesCellVect(:,4);
    end 
    
        
                    
%     Weigths = max(eps,normpdf(NewPotentialStatesCellVect(:,2),Observations(IndTime),Parameters.SigmaObs*Observations(IndTime)));
    Weigths = ones(NbResampling,1);
    for IndObservedVar = 1:length(Parameters.ObservedVariables)    
        Weigths = Weigths.*normpdf(NewPotentialStatesCellVect(:,ObservedVariables(IndObservedVar)),Observations(ObservedVariables(IndObservedVar), IndTime),Parameters.SigmaObs(ObservedVariables(IndObservedVar))*Observations(ObservedVariables(IndObservedVar),IndTime));
    end
    LikelihoodsNew(IndTime,:) = LikelihoodsNew(IndTime-1,:).*normpdf(NewPotentialStatesCellVect(:,ObservedVariables(IndObservedVar)),Observations(ObservedVariables(IndObservedVar), IndTime),Parameters.SigmaObs(ObservedVariables(IndObservedVar))*Observations(ObservedVariables(IndObservedVar),IndTime))';
%     if sum(Weigths)==0
%         Weigths = max(eps,normpdf(NewPotentialStatesCellVect(:,2),Observations(IndTime),Parameters.SigmaObs*Observations(IndTime)));
%         'you'
%         ResultSMC.NulWeigths = 1;
%     end
    
    Liks(IndTime) = (mean(Weigths));
    LogLik = LogLik + log(mean(Weigths));
         
    Weigths = Weigths/sum(Weigths);
    ResampledParticules(IndTime,:,:) = NewPotentialStatesCellVect';
    
    PotentialStatesCellVect = NewPotentialStatesCellVect;
       
    CumWeigths=cumsum(Weigths)/sum(Weigths);
    ProportionKept(IndTime) = length(unique(KeptSamples))/length(KeptSamples);
%     hist(Weigths)
%     title(IndTime)
%     pause()
end

for var = 1:length(PreviousState)
    PosteriorMeansRecord(var,:) = mean(ResampledParticules(:,var,:),3);
end

NbSamplesForEstimation = ones(1,length(ObservationInstants));
for i = 1:length(ObservationInstants)
    NbSamplesForEstimation(i) = length(unique(ResampledParticules(i,1,:)));
end


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
    ResampledParticules = ResampledParticules(:,:,KeptSamples);
    LikelihoodsNew = LikelihoodsNew(:,KeptSamples);
%     CompleteLambdas2 = CompleteLambdas2(KeptSamples,:);
    
    for IndTime = 1:length(ObservationInstants)
        for var =1:4
            q = quantile(ResampledParticules(IndTime,var,:),[0.025,0.975]);
            Posterior975Record(var,IndTime) = q(1);
            Posterior225Record(var,IndTime) = q(2);
            PosteriorMeansRecord(var,IndTime) = mean(ResampledParticules(IndTime,var,:));
        end
        NbKeptSamples(IndTime) = length(unique(ResampledParticules(IndTime,1,:)));
    end
end

FathersTab(:,end+1) = KeptSamples;

FathersTab2 = FathersTab;
RebuiltLambdas = zeros(size(CompleteLambdas2));
test = [];
for i = size(FathersTab,2):-1:2
    RebuiltLambdas(:,(i-2)*NbTSteps + 1: (i-1)*NbTSteps) = CompleteLambdas(FathersTab2(:,i),(i-2)*NbTSteps + 1: (i-1)*NbTSteps);
    FathersTab2(:,i-1) = FathersTab2(FathersTab2(:,i),i-1);
%     test(i) = sum(sum(RebuiltLambdas(:,(i-2)*NbTSteps + 1: (i-1)*NbTSteps) == CompleteLambdas(:,(i-2)*NbTSteps + 1: (i-1)*NbTSteps)));
end



if Task == 'Filter'
    ResultSMC.PosteriorMeansRecord = PosteriorMeansRecord;
    ResultSMC.NbKeptSamples = NbKeptSamples;
    ResultSMC.NbParticules = NbResampling;
    ResultSMC.Data = Data;
    ResultSMC.Parameters = Parameters;
    ResultSMC.Posterior975Record = Posterior975Record;
    ResultSMC.Posterior225Record = Posterior225Record;
    ResultSMC.ResampledParticules = ResampledParticules;
    ResultSMC.Likelihood = prod(Liks);
    ResultSMC.Liks = Liks;
    ResultSMC.LogLik = LogLik;
    ResultSMC.RealTime = toc;
    ResultSMC.CompleteLambdas = RebuiltLambdas;
    ResultSMC.LikelihoodsNew = LikelihoodsNew;
    ResultSMC.ProportionKept = ProportionKept;
else
    ResultSMC = - LogLik;
end
