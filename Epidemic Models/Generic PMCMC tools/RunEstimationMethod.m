function Res = RunEstimationMethod(Data, Model, Parameters, TempPar, NbIterations)

TransfThetas = [];
PropTransfThetas = [];
Thetas = [];
if not(Parameters.NoPaths)
    Paths = zeros(NbIterations,length(Parameters.PathsToKeep),sum(Data.NbComputingSteps)+1);
else
    Paths = [];
end
LogLiks = [];
PropLogLiks = [];
LogPosts = [];
PropLogPosts = [];
Accepted = [];
Ratios = [];
Grads = [];
ProposedPars = [];
SamplCovs = [];
Mus = [];
Inds = [];
Pars = [];
MargLogLiks = [];
MargLogLikStars = [];
Renewal = [];
EstimatedSigs = [];
Coalescence = [];
InitCov = Parameters.G^(-1);
try 
    if not(Parameters.AdMet)
        Parameters.AdMet = 0;
    end
catch
    Parameters.AdMet = 0;
end

try 
    disp(Parameters.EpsAdaptRate)
catch
    Parameters.EpsAdaptRate = 0.1;
end


try
    if not(Parameters.SaveForMarcMCMC)
        Parameters.SaveForMarcMCMC = 0;
    end
catch
    Parameters.SaveForMarcMCMC = 0;
end


if Parameters.SaveForMarcMCMC
    try
        load(Parameters.NameToSave);
        MargLogLiks = Res.MargLogLiks;
        MargLogLiksStars = Res.MargLogLiksStars;
        Ratios = Res.Ratios;
        TransfThetas  = Res.TransfThetas;
        Thetas  = Res.Thetas;
        LogLiks = Res.LogLiks;
        LogPosts = Res.LogPosts;
        Accepted = Res.Accepted;
        Pars = Res.Pars;  
        startind = size(TransfThetas,2);
        TempPar = Res.TempPar;
        if not(Parameters.NoPaths)
            Paths = Res.Paths;
        end
    catch
        startind = 1;
    end
else
    startind = 1;
end

if not(Parameters.NoPaths)
   Names = Parameters.Names.Estimated;
   for i = 1:length(Names)
      Parameters.(Names{i}).TransfValue = TempPar.(Names{i}).TransfValue;
      Parameters = UpdateParsTransfToNoTransf(Parameters);
   end   
   RandInd = ceil(rand(1,1)*Parameters.NbParticules);
   ResTemp = EstimationSMCsmoothGen(Data,Model,Parameters); 
   TempPar.Paths = ResTemp.CompletePaths(RandInd,:,:);
end


Names = Parameters.Names.Estimated;
if not(isfield(Parameters,'PMCMC'))
    Parameters.PMCMC = 'PMMH';
end

% [TempPar.LogLik TempPar.LogPost]
% 'go'



for IndIt = startind:NbIterations
    if strcmp(Parameters.PMCMC,'PMMH')
       TempPar.Cov = InitCov;
       if Parameters.AdMet
            if IndIt > 100
                if rand(1,1)>Parameters.AdMetBeta
                    try
                        chol(cov(TransfThetas'));
                        TempPar.Cov = cov(TransfThetas');
%                         disp('changed cov')
%                         disp(eig(TempPar.Cov))
                    end
                end
            end
        end
        
        StarPar = ProposeParameter(Data,  Model,Parameters, TempPar); 
        for i = 1:length(StarPar.Par)
            if not(isreal(StarPar.Par(i)))
                disp('complex argument, this shouldn''t happen')
            end
        end
        Pars(IndIt,:) = TempPar.Par;
        ProposedPars(IndIt,:) = StarPar.Par;
    
%         [StarPar.LogLik StarPar.LogPost]
        
        LogRand = log(rand(1,1));
        
        
        if LogRand <= StarPar.LogAccRate
            TempPar = StarPar;
            Accepted(IndIt) = 1;
        elseif and(StarPar.LogAccRate>-10,and(length(Accepted)>1000,mean(Accepted(end-min(1000,length(Accepted))+1:end))==0))
            TempPar = StarPar;
            Accepted(IndIt) = 1;
            disp('forced it')
        else
            Accepted(IndIt) = 0;
        end
        

        for i = 1:length(Names)
            TransfThetas(Parameters.(Names{i}).Index,IndIt) = TempPar.(Names{i}).TransfValue;
        end
        for i = 1:length(Names)
            PropTransfThetas(Parameters.(Names{i}).Index,IndIt) = StarPar.(Names{i}).TransfValue;
        end
        if not(Parameters.NoPaths)
            try
                Paths((IndIt-1)*1+1:(IndIt)*1,:,:) = TempPar.Paths;
            catch
                Paths(IndIt,:,:) = TempPar.Paths;
            end
        end
%         TempPar.LogPost
%         StarPar.LogPost
        LogLiks(IndIt) = TempPar.LogLik;
        LogPosts(IndIt) = TempPar.LogPost;
        PropLogLiks(IndIt) = StarPar.LogLik;
        PropLogPosts(IndIt) = StarPar.LogPost;
%         Coalescence(IndIt) = TempPar.Coalescence;
%         EstimatedSigs(IndIt) = TempPar.EstimatedSig;
        if rand(1,1)<0.1
            disp(['Acceptance Rate: ' num2str(sum(Accepted)/length(Accepted),2) ' (it' num2str(IndIt) ')'])
        end
%         if strcmp(Parameters.MCMCType, 'Inde')
%             Ratios(IndIt) = TempPar.Ratio ;
%         end        
        try 
            Grads(:,IndIt) = TempPar.Grad;
        end
        try
            Parameters.AccRate = sum(Accepted)/length(Accepted);        
        catch
            'AccRate Problem'
            die
            Parameters.AccRate = 100;
        end
        if IndIt>150
            if rand(1,1)<0.2
                Parameters.Epsil = exp(log(Parameters.Epsil) + Parameters.AdaptC^IndIt*(Parameters.AccRate-0.23));

            end
        end
    
    end
    if Parameters.SaveForMarcMCMC
        Res.MargLogLiks = MargLogLiks;
        Res.MargLogLiksStars = MargLogLikStars;
        Res.Epsil = Parameters.Epsil;
        Res.Ratios = Ratios;
        Res.TempPar = TempPar;
        Res.TransfThetas  = TransfThetas;
        Res.Thetas  = Thetas;
        Res.Paths = Paths;
        Res.LogLiks = LogLiks;
        Res.LogPosts = LogPosts;
        Res.AccRate = sum(Accepted)/length(Accepted);
        Res.Accepted = Accepted;
        Res.Grads = Grads;
        Res.ESS = [];
        Res.CompInds = Inds;
        Res.Pars = Pars;
        Res.Parameters = Parameters;
        Res.Data = Data;
        Res.Model = Model;
        save(Parameters.NameToSave,'Res');
    end
end


for IndIt = 1:size(TransfThetas,2)
    for i = 1:length(Names)
        Parameters.(Names{i}).TransfValue = TransfThetas(Parameters.(Names{i}).Index,IndIt);
    end
    Parameters = UpdateParsTransfToNoTransf(Parameters);
    for i = 1:length(Names)
        Thetas(Parameters.(Names{i}).Index,IndIt) = Parameters.(Names{i}).Value ;
    end
end


Res.MargLogLiks = MargLogLiks;
Res.MargLogLiksStars = MargLogLikStars;
Res.Epsil = Parameters.Epsil;
Res.Ratios = Ratios;
Res.TempPar = TempPar;
Res.TransfThetas  = TransfThetas;
%Res.PropTransfThetas  = PropTransfThetas;
Res.Thetas  = Thetas;
Res.Paths = Paths;
%Res.PropLogLiks = PropLogLiks;
Res.LogLiks = LogLiks;
%Res.PropLogPosts = PropLogPosts;
Res.LogPosts = LogPosts;
Res.AccRate = sum(Accepted)/length(Accepted);
Res.Accepted = Accepted;
Res.Grads = Grads;
Res.ESS = [];
%Res.ProposedPars = ProposedPars;
Res.CompInds = Inds;
Res.Pars = Pars;
Res.Parameters = Parameters;
Res.Data = Data;
Res.Model = Model;

tmp = median(Thetas');
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).Value = tmp(Parameters.(Names{i}).Index);
end
Parameters = UpdateParsNoTransfToTransf(Parameters);

if strcmp(Parameters.ModelType,'SMC')
    TempThetaMean = EstimationSMCsmoothGen(Data, Model, Parameters);
    p_d = mean(-2*LogLiks) - (-2)*TempThetaMean.LogLik;
    DIC = (-2)*TempThetaMean.LogLik + 2*p_d;

    Res.p_d = p_d;
    Res.DIC = DIC;
end

for i = 1:size(TransfThetas,1)
    temp = AutoCorrelation(TransfThetas(i,:));
    Res.ESSTransf(i) = NbIterations/(1+2*sum(temp(2:end)));
    Res.RelESSTransf(i) = Res.ESSTransf(i)/NbIterations*100;
end
for i = 1:size(TransfThetas,1)
    temp = AutoCorrelation(Thetas(i,:));
    Res.ESS(i) = NbIterations/(1+2*sum(temp(2:end)));
    Res.RelESS(i) = Res.ESS(i)/NbIterations*100;
end

if strcmp(Parameters.PMCMC,'Gibbs')                
    Res.GibbsAccRates = Parameters.GibbsAccRates;
    Res.EstimatedSigs = EstimatedSigs; 
end
    
    

try
    if Parameters.SaveSpace
        NbWeKeep = min(5000,round(NbIterations/10));
        inds = sort(randsample(NbIterations,NbWeKeep));
        Res.Thetas = Res.Thetas(:,inds);
        Res.TransfThetas = Res.TransfThetas(:,inds);
        Res.Paths = Res.Paths(inds,:,:);
        Res.LogLiks = Res.LogLiks(inds);
        Res.LogPosts = Res.LogPosts(inds);
    end
end
    