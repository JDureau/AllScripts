function ResultSMC = EstimationSMCsmoothGen(Data, Model, Parameters)


Parameters = Model.InitializeParameters(Parameters);

% CheckParametersGen(Parameters)

NbResampling = Parameters.NbParticules;


ind = 1;
Pars = zeros(1,Parameters.NbParsEstimated);
Names = Parameters.Names.Estimated;
for i = 1:length(Names);
    Parameters.(Names{i}).Ind = ind;
    Pars(ind) = Parameters.(Names{i}).TransfValue;
    ind = ind+1;
end
NbParsEst = ind - 1;
    
NbVarsTot = length(Parameters.InitialState);
ObservationInstants = zeros(1,length(Data.Y));
NbComputingSteps = Data.NbComputingSteps;

InitialState = Parameters.InitialState;
InitialStates = repmat(InitialState,1,NbResampling);

PosteriorMeansRecord = diag(InitialState)*ones(NbVarsTot,length(ObservationInstants));
NbKeptSamples = NbResampling*ones(1,length(ObservationInstants));
ProportionKept = [];

ResampledParticules = zeros(length(ObservationInstants),NbVarsTot,NbResampling);
Particules =  zeros(length(ObservationInstants),NbVarsTot,NbResampling);
Particules(1,:,:) = InitialStates;

 
Variables = InitialStates';

CumWeigths=(1:NbResampling)/NbResampling;
Weigths = ones(1,NbResampling);

Liks = ones(1,length(ObservationInstants));
LogLik = 0;
FathersTab = zeros(NbResampling,length(ObservationInstants)-1);

if not(Parameters.NoPaths)
    Paths = zeros(NbResampling,length(Parameters.PathsToKeep),sum(Data.NbComputingSteps));
else
    Paths = [];
end


for IndTime = 2:length(ObservationInstants)
    
    
    RandU1 = rand(1,1)/NbResampling;

    
    if and(not(sum(Weigths)==0),not(isnan(Weigths)))
        Weigths = Weigths/sum(Weigths);
    else
        Weigths = ones(size(Weigths))/NbResampling;
    end
    ESS = 1/sum(Weigths.^2);
    
    KeptSamples = 1:NbResampling;   
    if ESS<NbResampling/2
        u = rand(1,1)/Parameters.NbParticules;
        s = 0;
        KeptSamples = [];
        resind = 1;
        for ipart = 1:Parameters.NbParticules
            k = 0;
            s = s+Weigths(ipart);
            while s>u
                k=k+1;
                u = u+1/Parameters.NbParticules;
                KeptSamples(resind) = ipart;
                resind = resind+1;
            end
        end
    end
    
    if not(length(KeptSamples) == Parameters.NbParticules)
        'stop'
    end
    
    % Bug in there:
%     if ESS<NbResampling/2
%         
%         Current = 1;
%         for IndResampling = 1:NbResampling
%             while CumWeigths(Current)<=RandU1+(IndResampling-1)/NbResampling
%                 Current = Current+1;
%             end
%             KeptSamples(IndResampling) = Current;
%         end
%     else
%         KeptSamples = 1:NbResampling;
%     end
%     try
%         % This is for GIBBS PMCMC
%         if Parameters.ForceTraj
%             KeptSamples(1) = 1;
%         end
%     end
    
    if length(unique(KeptSamples))<NbResampling/10
%         disp('Crash (<10%Parts)')
    end
    
    tmp = mean(Variables)';
    
    Variables = Variables(KeptSamples,:);
   
%     if IndTime == 3
%         [tmp mean(Variables)']
%         die
%     end
    
    
    
    try
        FathersTab(:,IndTime-1) = KeptSamples;
    catch
        'stop'
    end
    Res = Model.SMC_projection(Variables,IndTime,NbComputingSteps(IndTime),Data,Parameters,Paths);
        
%     mean(Res.Variables)'
%     diag(cov(Variables))
%     die
        
    Variables = Res.Variables;
    Paths = Res.Paths;
    
    
    
%     if IndTime == 2
%         disp(mean(Variables)')
%         die
%     end


    try
        coeff = Parameters.MultCoeff.Value/10;
    catch
        coeff = 1;
    end
    
    
    if 0%and(strcmp(Parameters.Problem,'ImperialHIV2'),length(Data.ObservedVariables{IndTime}) == 2)
        Weigths = eval(Model.LikFunction1);
        Weigths = Weigths.*eval(Model.LikFunction2);
    else
        Weigths = eval(Model.LikFunction);
    end
%     if IndTime == 2
%         disp(IndTime)
%         disp(Data.Y(IndTime)-Data.Y(IndTime-1))
%         disp(mean(Variables(:,2)))
%         disp(mean(Weigths))
%    end
    
    
    Liks(IndTime) = (mean(Weigths));
    LogLik = LogLik + max(-700,log(mean(Weigths)));
    if isinf(log(mean(Weigths)))
%         disp(IndTime)
%         disp('Strong Crash (-Inf)')
    end
    
%     if mean(Weigths) == 0
%         ResultSMC.LogLik = -Inf;
%         ResultSMC.Grad = [];
%         ResultSMC.G = [];
%         ResultSMC.CompletePaths = [];
%         ResultSMC.Cov = [];
%         ResultSMC.CompInd = 1;
%         return
%     end
        
         
    Weigths = Weigths/sum(Weigths);
    Particules(IndTime,:,:) = Variables';
     
    CumWeigths=cumsum(Weigths)/sum(Weigths);
%     ProportionKept(IndTime) = length(unique(KeptSamples))/length(KeptSamples);
end

RandU1 = rand(1,1)/NbResampling;
KeptSamples = zeros(1,NbResampling);

Current = 1;
for IndResampling = 1:NbResampling
    while CumWeigths(Current)<=RandU1+(IndResampling-1)/NbResampling
        Current = Current+1;
    end
    KeptSamples(IndResampling) = Current;
end
try
    % This is for GIBBS PMCMC
    if Parameters.ForceTraj
        KeptSamples(1) = 1;
    end
end    
FathersTab(:,end+1) = KeptSamples;

FathersTab2 = FathersTab;
if not(Parameters.NoPaths)
    RebuiltPaths = zeros(size(Paths));
end
% test = [];
for i = size(FathersTab,2):-1:2
    if not(Parameters.NoPaths)
        RebuiltPaths(:,:,sum(Data.NbComputingSteps(1:i-1)) + 1: sum(Data.NbComputingSteps(1:i))) = Paths(FathersTab2(:,i),:,sum(Data.NbComputingSteps(1:i-1)) + 1: sum(Data.NbComputingSteps(1:i)));
    end
    ResampledParticules(i,:,:) = Particules(i,:,FathersTab2(:,i));
    FathersTab2(:,i-1) = FathersTab2(FathersTab2(:,i),i-1);
%     test(i) = sum(sum(RebuiltLambdas(:,(i-2)*NbTSteps + 1: (i-1)*NbTSteps) == CompleteLambdas(:,(i-2)*NbTSteps + 1: (i-1)*NbTSteps)));
end
ResampledParticules(1,:,:) = Particules(1,:,FathersTab2(:,1));

for IndTime = 1:length(ObservationInstants)
    for ivar =1:Parameters.NbVariables
%         q = quantile(ResampledParticules(IndTime,ivar,:),[0.025,0.975]);
%         Posterior975Record(ivar,IndTime) = q(1);
%         Posterior225Record(ivar,IndTime) = q(2);
        PosteriorMeansRecord(ivar,IndTime) = mean(ResampledParticules(IndTime,ivar,:));
    end
    NbKeptSamples(IndTime) = length(unique(ResampledParticules(IndTime,1,:)));
end



if strcmp(Parameters.MCMCType, 'Lang')
    if strcmp(Parameters.GMeth, 'PseudoCst')
        NbParsEst = length(Parameters.Pars);
        Grad = ComputeFirstDerivative(Parameters);
        G = -ComputeSecondDerivative(Parameters);
    end
    if strcmp(Parameters.GradMeth,'Kalman')
        ResGrad = KalmanNumericDerivativesGen(Data,Parameters);
        Grad = ResGrad.Grad;
        if not(isreal(Grad))
            disp('Problem: complex value')
        end
    end
end


ResultSMC.PosteriorMeansRecord = PosteriorMeansRecord;
% plot(PosteriorMeansRecord(5,:))
% title(LogLik)
% pause(0.01)
ResultSMC.NbKeptSamples = NbKeptSamples;
ResultSMC.NbParticules = NbResampling;
ResultSMC.Data = Data;
ResultSMC.Parameters = Parameters;
% ResultSMC.Posterior975Record = Posterior975Record;
% ResultSMC.Posterior225Record = Posterior225Record;
ResultSMC.ResampledParticules = ResampledParticules;
ResultSMC.Likelihood = prod(Liks);
ResultSMC.Liks = Liks;
if isnan(LogLik)
    LogLik = -Inf;
end
if isinf(LogLik)
    disp('let s see');
end
ResultSMC.LogLik = LogLik;

if not(Parameters.NoPaths)
    ResultSMC.CompletePaths(:,:,1) = InitialStates(Parameters.PathsToKeep,:)';
    ResultSMC.CompletePaths(:,:,2:size(RebuiltPaths,3)+1) = RebuiltPaths;
end
ResultSMC.ProportionKept = ProportionKept;

try
    ResultSMC.Grad = Grad;
catch
    ResultSMC.Grad = [];
end
try        
    ResultSMC.Hess = Hess;
catch
    ResultSMC.Hess = [];
end

    
Thetas = zeros(Parameters.NbParsEstimated,1);
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Thetas(i,1) = Parameters.(Names{i}).Value;
end
ResultSMC.Thetas = Thetas;
    
if strcmp(Parameters.MCMCType, 'Lang')
    if strcmp(Parameters.GMeth, 'one')
        ResultSMC.G = eye(NbParsEst);
        ResultSMC.GParDer = zeros(NbParsEst,NbParsEst);
    elseif strcmp(Parameters.GMeth, 'cst given')
        ResultSMC.G = Parameters.G;
        ResultSMC.GParDer = zeros(NbParsEst,NbParsEst);
    elseif strcmp(Parameters.GMeth, 'Kal')
        
    elseif strcmp(Parameters.GMeth, 'ful')
        ResultSMC.G = G;
        ResultSMC.GParDer = GParDer;
    elseif strcmp(Parameters.GMeth, 'PseudoCst')
        ResultSMC.G = G;
        ResultSMC.GParDer = zeros(NbParsEst,NbParsEst);
    end
else
    ResultSMC.G = [];
    ResultSMC.GParDer = [];
end
if strcmp(Parameters.MCMCType, 'GMM')
    if strcmp(Parameters.GMeth, 'GMM')     
        ind = PickRandComp(Scale((Pars)',Parameters),Parameters.DensityModel, Parameters);
        ResultSMC.Cov = diag(Parameters.ScalingStds)*squeeze(Parameters.DensityModel.Sigma(:,:,ind))*diag(Parameters.ScalingStds);
        ResultSMC.CompInd = ind;
    end
end
if strcmp(Parameters.MCMCType, 'Rand')
    if strcmp(Parameters.GMeth, 'cst given')
        ResultSMC.Cov = Parameters.G^-1;
        ResultSMC.CompInd = 1;
    end
    if strcmp(Parameters.GMeth, 'GMM')     
        ind = PickRandComp(Scale((Pars)',Parameters),Parameters.DensityModel, Parameters);
        ResultSMC.Cov = diag(Parameters.ScalingStds)*squeeze(Parameters.DensityModel.Sigma(:,:,ind))*diag(Parameters.ScalingStds);
        ResultSMC.CompInd = ind;
    end
end
    
LogPrior = 0;
LogCorr = 0;
Names = Parameters.Names.Estimated;

for i = 1:length(Names)
    temp = Parameters.(Names{i}).Prior(Names{i},Parameters);
    LogPrior = LogPrior +log(temp);
    temp = Parameters.(Names{i}).CorrFunct(Names{i},Parameters);
    LogCorr = LogCorr - log(temp);
end
ResultSMC.LogPrior = LogPrior;
ResultSMC.LogCorr = LogCorr;


% ResultSMC.Coalescence = mean(sum((squeeze(ResultSMC.CompletePaths(:,1,:)==repmat(ResultSMC.CompletePaths(1,1,:),NbResampling,1)))'));
    
    