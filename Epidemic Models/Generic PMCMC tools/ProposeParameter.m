function StarPar = ProposeParameter(Data, Model, Parameters, TempPar)

% The new parameter needs to have the following fields:
%  - its value
%  - its LogLik
%  - its LogAccRate
%  - an associated trajectory for lambda
%  - potential information to propose a new one



% Step 1: sample star parameter
ParametersStar = Parameters;
LogPrior = 0;
LogCorr = 0;
Names = Parameters.Names.Estimated;
if strcmp(Parameters.MCMCType,'Lang')
    epsil = Parameters.Epsil;
    CholCov = chol(2.38^2/length(Names)*TempPar.Cov);
    StarPar.Par = mvnrnd(TempPar.Mu,TempPar.Cov);
    for i = 1:length(Names)
        StarPar.(Names{i}).TransfValue =  StarPar.Par(Parameters.(Names{i}).Index);
        ParametersStar.(Names{i}).TransfValue = StarPar.(Names{i}).TransfValue;
        ParametersForPriors = UpdateParsTransfToNoTransf(ParametersStar);
        LogPrior = LogPrior + log(eval(Parameters.(Names{i}).Prior));
        LogCorr = LogCorr + log(eval(Parameters.(Names{i}).Corr));
    end
elseif strcmp(Parameters.MCMCType,'GMM')
    epsil = Parameters.Epsil;
    ind = PickRandComp(Scale((TempPar.Par),Parameters),Parameters.DensityModel, Parameters);
    Sigma = diag(Parameters.ScalingStds)*squeeze(Parameters.DensityModel.Sigma(:,:,ind))*diag(Parameters.ScalingStds);
    StarPar.Par = mvnrnd(TempPar.Par + Parameters.Epsil^2/2*(ScaleBack(Parameters.DensityModel.mu(ind,:)',Parameters)-TempPar.Par),epsil^2*Sigma);
    for i = 1:length(Names)
        StarPar.(Names{i}).TransfValue =  StarPar.Par(Parameters.(Names{i}).Index);
        ParametersStar.(Names{i}).TransfValue = StarPar.(Names{i}).TransfValue;
        ParametersForPriors = UpdateParsTransfToNoTransf(ParametersStar);
        LogPrior = LogPrior + log(eval(Parameters.(Names{i}).Prior));
    end  
elseif strcmp(Parameters.MCMCType,'Rand')
    Cov = TempPar.Cov;
    CholCov = chol(Cov);
    try
        StarPar.Par = mvnrnd(TempPar.Par,Parameters.Epsil^2*2.38^2/length(Names)*Cov);
    catch
        'sampling problem in Propose parameter'
        StarPar.Par = [0 ];
    end
    for i = 1:length(Names)
        StarPar.(Names{i}).TransfValue =  StarPar.Par(Parameters.(Names{i}).Index);
        ParametersStar.(Names{i}).TransfValue = StarPar.(Names{i}).TransfValue;
        ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);
        %tmp = Parameters.(Names{i}).Prior(Names{i},ParametersStar);
        %LogPrior = LogPrior + log(tmp);
        %LogCorr = LogCorr + log(Parameters.(Names{i}).CorrFunct(Names{i},ParametersStar));
    end       
elseif strcmp(Parameters.MCMCType, 'Inde')
     StarPar.Par = random(Parameters.DensityModel , 1);
%      Cov = TempPar.Cov;
%     CholCov = chol(Cov);
%     StarPar.Par = mvnrnd(TempPar.Par,Parameters.Epsil^2*2.38^2/length(Names)*Cov);

     for i = 1:length(Names)
        StarPar.(Names{i}).TransfValue =  StarPar.Par(Parameters.(Names{i}).Index);
        ParametersStar.(Names{i}).TransfValue = StarPar.(Names{i}).TransfValue;
        ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);
        %tmp = Parameters.(Names{i}).Prior(Names{i},ParametersStar);
        %LogPrior = LogPrior + log(tmp);
        %LogCorr = LogCorr + log(Parameters.(Names{i}).CorrFunct(Names{i},ParametersStar));
    end 
end
ParametersStar.Pars = StarPar.Par;
ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);
    
    
% Step 2: test it
if strcmp(Parameters.ModelType,'SMC')   
    ResStar = EstimationSMCsmoothGen(Data,Model,ParametersStar);
elseif strcmp(Parameters.ModelType,'Kalman')   
    try
        ResStar = EstimationEKFGen(Data, Model, ParametersStar);
%     catch
%         'ouch'
%         ResStar = EstimationSMCsmoothGen(Data,Model,ParametersStar);
    end
    ResStar.Cov = TempPar.Cov;
end
ParametersStar = ResStar.Parameters;
for i = 1:length(Names)
    StarPar.(Names{i}).TransfValue = ParametersStar.(Names{i}).TransfValue;
    StarPar.Par(Parameters.(Names{i}).Index) = StarPar.(Names{i}).TransfValue;
end
% GradStar = ResStar.Grad;
% GStar = ResStar.G;

StarPar.Paths = [];
RandInd = ceil(rand(1,1)*Parameters.NbParticules);
if not(Parameters.NoPaths)
    try
        StarPar.Paths = ResStar.CompletePaths(RandInd,:,:);
        if Parameters.KeepAll
            inds = ceil(Parameters.NbParticules*rand(1,1));
            StarPar.AllPaths = ResStar.CompletePaths(inds,:,:);
        end
%     catch
%         StarPar.Paths = [];
    end
end
LogCorr = ResStar.LogCorr;
LogPrior = ResStar.LogPrior;
StarPar.LogLik = ResStar.LogLik;
StarPar.LogPost = ResStar.LogLik + LogPrior - LogCorr;
StarPar.LogCorr = LogCorr;
StarPar.LogPrior = LogPrior;


% if (StarPar.LogLik>-415)
%     disp('yo')
% end


% Acceptance Ratios
if strcmp(Parameters.MCMCType, 'Lang')
    StarPar.Grad = GradStar;
    StarPar.Mu = StarPar.Par + epsil^2/2*GStar^-1*GradStar' ;%  - epsil^2/2*GStar^-2*GParDerStar ;
    if not(isreal(StarPar.Mu))
        disp('Problem: complex value')
    end
    StarPar.Cov = epsil^2*GStar^-1;
    if prod(eig(StarPar.Cov)) <= 0
        disp('Problem with cov')
        [V,D] = eig(StarPar.Cov) ;
        EigValues = [];
        for i = 1:size(D,1)
            EigValues(i) = D(i,i);
        end
        EigValues = max(eps,EigValues);
        StarPar.Cov = V*diag(EigValues)*V^-1;
    end
    StarPar.CompInd = ResStar.CompInd;
    StarPar.MuSampledFrom = TempPar.Mu;
    StarPar.CovSampledFrom = TempPar.Cov;
    qStarTemp = mvnpdf(StarPar.Par, TempPar.Mu, TempPar.Cov);
    qTempStar = mvnpdf(TempPar.Par, StarPar.Mu, StarPar.Cov);
    StarPar.LogAccRate = StarPar.LogPost + log(qTempStar) - TempPar.LogPost - log(qStarTemp);
elseif strcmp(Parameters.MCMCType, 'GMM')
    % Temp to Star 
    temp = 0;
    post = posterior(Parameters.DensityModel,TempPar.Par');
    if not(sum(post) == 1)
        disp('pb')
    end
    for IndComp = 1:Parameters.DensityModel.NComponents
        temp = temp + (mvnpdf(StarPar.Par,TempPar.Par + Parameters.Epsil^2/2*(ScaleBack(Parameters.DensityModel.mu(IndComp,:)',Parameters)-TempPar.Par),Parameters.Epsil^2*diag(Parameters.ScalingStds)*squeeze(Parameters.DensityModel.Sigma(:,:,IndComp))*diag(Parameters.ScalingStds))*post(IndComp));
    end
    LogqTempStar = log(temp);
    % Star to Temp 
    temp = 0;
    post = posterior(Parameters.DensityModel,StarPar.Par');
    if not(sum(post) == 1)
        disp('pb')
    end
    for IndComp = 1:Parameters.DensityModel.NComponents
        temp = temp + (mvnpdf(TempPar.Par,StarPar.Par + Parameters.Epsil^2/2*(ScaleBack(Parameters.DensityModel.mu(IndComp,:)',Parameters)-StarPar.Par),Parameters.Epsil^2*diag(Parameters.ScalingStds)*squeeze(Parameters.DensityModel.Sigma(:,:,IndComp))*diag(Parameters.ScalingStds))*post(IndComp));
    end
    LogqStarTemp = log(temp);
    StarPar.LogAccRate = StarPar.LogPost + LogqStarTemp - TempPar.LogPost - LogqTempStar;
    StarPar.Cov = ResStar.Cov;
    StarPar.CovSampledFrom = Parameters.Epsil^2*ResStar.Cov;
    StarPar.MuSampledFrom = TempPar.Par;
    StarPar.CompInd = ResStar.CompInd;
elseif strcmp(Parameters.MCMCType, 'Rand')    
    StarPar.LogAccRate = StarPar.LogPost  - ( TempPar.LogPost  ); 
    StarPar.Cov = ResStar.Cov;
    StarPar.CovSampledFrom = Parameters.Epsil^2*Cov;
    StarPar.MuSampledFrom = TempPar.Par;
%     StarPar.CompInd = ResStar.CompInd;
%     StarPar.Coalescence = ResStar.Coalescence;
elseif strcmp(Parameters.MCMCType, 'Inde')
    StarPar.LogAccRate = StarPar.LogPost  - ( TempPar.LogPost  ); 
%     StarPar.Cov = ResStar.Cov;
%     StarPar.CovSampledFrom = Parameters.Epsil^2*Cov;
%     StarPar.MuSampledFrom = TempPar.Par;
    
    
    qStarTemp = pdf(Parameters.DensityModel, StarPar.Par);
    qTempStar = pdf(Parameters.DensityModel, TempPar.Par);
%     'loglik'
%     StarPar.LogLik  - TempPar.LogLik
%     'logpost'
%     StarPar.LogPost  - TempPar.LogPost
%     'logq'
%     log(qTempStar)  - log(qStarTemp)
    StarPar.LogAccRate = StarPar.LogPost + log(qTempStar) - TempPar.LogPost - log(qStarTemp);
    StarPar.Ratio = qTempStar/exp(StarPar.LogLik);
end