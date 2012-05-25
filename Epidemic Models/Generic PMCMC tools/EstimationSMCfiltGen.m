function ResultSMC = EstimationSMCfiltGen(Data, Model, Parameters)

Parameters = Model.InitializeParameters(Parameters);

CheckParametersGen(Parameters)

NbResampling = Parameters.NbParticules;


ind = 1;
Pars = [];
Names = Parameters.Names.Estimated;
for i = 1:length(Names);
    Parameters.(Names{i}).Ind = ind;
    Pars(ind) = Parameters.(Names{i}).TransfValue;
    ind = ind+1;
end
NbParsEst = ind - 1;
    
NbVarsTot = Parameters.NbVariables;
Observations = Data.Observations;
ObservationInstants = Data.Instants;
NbComputingSteps = Data.NbComputingSteps;

Variables = ones(NbResampling,NbVarsTot);
InitialState = Parameters.InitialState;
PreviousState = InitialState;
PosteriorMeansRecord = diag(InitialState)*ones(NbVarsTot,length(ObservationInstants));
NbKeptSamples = NbResampling*ones(1,length(ObservationInstants));
ProportionKept = [];

ResampledParticules = zeros(length(ObservationInstants),NbVarsTot,NbResampling);


if not(Parameters.NoPaths)
    CompletePaths= zeros(NbResampling,length(Parameters.PathsToKeep),sum(Data.NbComputingSteps));
end
    
for IndResampling = 1:NbResampling
    Variables(IndResampling,:)=PreviousState;
    if not(Parameters.NoPaths)
        CompletePaths(IndResampling,:,1) = PreviousState(Parameters.PathsToKeep);
    end
end
currentind = 1;
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
    KeptSamples = zeros(1,NbResampling);
    
    Weigths = Weigths/sum(Weigths);
    ESS = 1/sum(Weigths.^2);
    
    if ESS<NbResampling/2
        Current = 1;
        for IndResampling = 1:NbResampling
            while CumWeigths(Current)<=RandU1+(IndResampling-1)/NbResampling
                Current = Current+1;
            end
            KeptSamples(IndResampling) = Current;
        end
    else
        KeptSamples = 1:NbResampling;
    end
    
    Variables = Variables(KeptSamples,:);
    
    FathersTab(:,IndTime-1) = KeptSamples;

%     if IndTime>2
%         Paths(:,:,sum(Data.NbComputingSteps(1:IndTime-1))) = Paths(KeptSamples,:,sum(Data.NbComputingSteps(1:IndTime-1)));
%     end

    ResampledParticules(IndTime-1,:,:) = Variables';

    Res = Model.SMC_projection(Variables,IndTime,NbComputingSteps(IndTime),Data,Parameters,Paths);
    
    Variables = Res.Variables;
    Paths = Res.Paths;
    

%     for i = size(CompletePaths,3)+1:size(CompletePaths,3)+size(Paths,3)
%         CompletePaths(:,:,i) = Paths(:,:,i);
%     end
    try
        coeff = Parameters.MultCoeff.Value/10;
    catch
        coeff = 1;
    end

    try
        Weigths = eval(Model.LikFunction);
    catch
        disp('prob')
    end
        
    Liks(IndTime) = (mean(Weigths));
    LogLik = LogLik + log(mean(Weigths));
    if isinf(LogLik)
%         disp('InfLogLik')
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
    
     
    CumWeigths=cumsum(Weigths)/sum(Weigths);
    ProportionKept(IndTime) = length(unique(KeptSamples))/length(KeptSamples);
end


for IndTime = 1:length(ObservationInstants)
    for ivar =1:Parameters.NbVariables
        q = quantile(ResampledParticules(IndTime,ivar,:),[0.025,0.975]);
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

Thetas = zeros(Parameters.NbParsEstimated,Parameters.NbParticules);
Names = Parameters.Names.Estimated;
for i = 1:Parameters.NbParsEstimated
    Thetas(i,:) = repmat(Parameters.(Names{i}).Value,1,Parameters.NbParticules);
end
    
ResultSMC.Parameters = Parameters;
ResultSMC.Thetas = Thetas;
ResultSMC.PosteriorMeansRecord = PosteriorMeansRecord;
ResultSMC.Paths = Paths;
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
ResultSMC.LogLik = LogLik;

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
    
    
    