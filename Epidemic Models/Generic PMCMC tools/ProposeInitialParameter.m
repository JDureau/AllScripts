function StarPar = ProposeInitialParameter(Data, Model, Parameters)


% Step 1: sample star parameter
ParametersStar = Parameters;
epsil = Parameters.Epsil;
Prior = 1;
Names = Parameters.Names.Estimated;
if strcmp(Parameters.MCMCType, 'Lang')
    for i = 1:length(Names)
        StarPar.(Names{i}).TransfValue =  Parameters.(Names{i}).TransfValue;
        ParametersStar.(Names{i}).TransfValue = StarPar.(Names{i}).TransfValue;
        StarPar.Par(Parameters.(Names{i}).Index,1) = StarPar.(Names{i}).TransfValue;
    end
elseif strcmp(Parameters.MCMCType, 'GMM')
    for i = 1:length(Names)
        ind = Parameters.(Names{i}).Index;
        mu = ScaleBack(Parameters.DensityModel.mu(1,:)',Parameters);
        StarPar.(Names{i}).TransfValue =  mu(ind);
        ParametersStar.(Names{i}).TransfValue = StarPar.(Names{i}).TransfValue;
        StarPar.Par(Parameters.(Names{i}).Index,1) = StarPar.(Names{i}).TransfValue;
    end
elseif strcmp(Parameters.MCMCType, 'Rand')
    Cov = Parameters.G^-1;
    StarPar.Par = [0];
    for i = 1:length(Names)
        StarPar.(Names{i}).TransfValue =  Parameters.(Names{i}).TransfValue;
        ParametersStar.(Names{i}).TransfValue = StarPar.(Names{i}).TransfValue;
        StarPar.Par(Parameters.(Names{i}).Index,1) = StarPar.(Names{i}).TransfValue;
    end
elseif strcmp(Parameters.MCMCType, 'Inde')
    StarPar.Par = random(Parameters.DensityModel , 1)';  
    for i = 1:length(Names)
        StarPar.(Names{i}).TransfValue =  Parameters.(Names{i}).TransfValue;
        ParametersStar.(Names{i}).TransfValue = StarPar.(Names{i}).TransfValue;
    end
    qStarTemp = pdf(Parameters.PosteriorDensityModel, Scale(StarPar.Par,Parameters)');
end
ParametersStar.Pars = StarPar.Par;
ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);

% Step 2: test it
if strcmp(Parameters.ModelType,'SMC')   
    ResStar = EstimationSMCsmoothGen(Data,Model,ParametersStar);
elseif strcmp(Parameters.ModelType,'MCMC')
    ResStar = EstimationLogLikDeterm(Data,Model,ParametersStar);
end
% disp(ResStar.LogLik)
GradStar = ResStar.Grad;
GStar = ResStar.G;
StarPar.LogLik = -Inf;
StarPar.LogPost = -Inf;
StarPar.LogCorr = -Inf;


RandInd = ceil(rand(1,1)*Parameters.NbParticules);
if not(Parameters.NoPaths)
    try
        StarPar.Paths = squeeze(ResStar.CompletePaths(RandInd,:,:));
    catch
        StarPar.Paths = [];
    end
end


% Acceptance Ratios
if strcmp(Parameters.MCMCType, 'Lang')
    StarPar.Mu = StarPar.Par + epsil^2/2*GStar^-1*GradStar' ;%  - epsil^2/2*GStar^-2*GParDerStar ;
    StarPar.Cov = epsil^2*GStar^-1;
elseif strcmp(Parameters.MCMCType, 'Inde')
    StarPar.Ratio = exp(StarPar.LogPost)/qStarTemp;
elseif strcmp(Parameters.MCMCType, 'Rand')
    StarPar.Cov = Cov;
end