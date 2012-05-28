function Result = RunEstimationMethod(Data, Model, Parameters, TempPar, NbIterations)

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

for IndIt = 1:NbIterations
    if strcmp(Parameters.PMCMC,'PMMH')
       TempPar.Cov = InitCov;
       if Parameters.AdMet
            if IndIt > 100
                if rand(1,1)>Parameters.AdMetBeta
                    try
                        if not(sum(eig(cov(TransfThetas'))<=eps))
                            TempPar.Cov = cov(TransfThetas');
%                             eig(TempPar.Cov)
                        end
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
%         elseif and(StarPar.LogAccRate>-50,mean(Accepted(end-min(100,length(Accepted))+1:end))==0)
%             TempPar = StarPar;
%             Accepted(IndIt) = 1;
%             disp('forced it')
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
        LogLiks(IndIt) = TempPar.LogLik;
        LogPosts(IndIt) = TempPar.LogPost;
        PropLogLiks(IndIt) = StarPar.LogLik;
        PropLogPosts(IndIt) = StarPar.LogPost;
%         Coalescence(IndIt) = TempPar.Coalescence;
%         EstimatedSigs(IndIt) = TempPar.EstimatedSig;
        if rand(1,1)<0.1
            disp(['Acceptance Rate: ' num2str(sum(Accepted)/length(Accepted),2) ' (it' num2str(IndIt) ')'])
        end
        if strcmp(Parameters.MCMCType, 'Inde')
            Ratios(IndIt) = TempPar.Ratio ;
        end        
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
        if IndIt>100
%             if rand(1,1)<0.2
                Parameters.Epsil = exp(log(Parameters.Epsil) + Parameters.AdaptC^IndIt*(Parameters.AccRate-0.23));

%             end
        end

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


Result.MargLogLiks = MargLogLiks;
Result.MargLogLiksStars = MargLogLikStars;
Result.Epsil = Parameters.Epsil;
Result.Ratios = Ratios;
Result.TempPar = TempPar;
Result.TransfThetas  = TransfThetas;
%Result.PropTransfThetas  = PropTransfThetas;
Result.Thetas  = Thetas;
Result.Paths = Paths;
%Result.PropLogLiks = PropLogLiks;
Result.LogLiks = LogLiks;
%Result.PropLogPosts = PropLogPosts;
Result.LogPosts = LogPosts;
Result.AccRate = sum(Accepted)/length(Accepted);
Result.Accepted = Accepted;
Result.Grads = Grads;
Result.ESS = [];
%Result.ProposedPars = ProposedPars;
Result.CompInds = Inds;
Result.Pars = Pars;
Result.Parameters = Parameters;
Result.Data = Data;
Result.Model = Model;

tmp = median(Thetas');
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).Value = tmp(Parameters.(Names{i}).Index);
end
Parameters = UpdateParsNoTransfToTransf(Parameters);
TempThetaMean = EstimationSMCsmoothGen(Data, Model, Parameters);

p_d = mean(-2*LogLiks) - (-2)*TempThetaMean.LogLik;
DIC = (-2)*TempThetaMean.LogLik + 2*p_d;

Result.p_d = p_d;
Result.DIC = DIC;

for i = 1:size(TransfThetas,1)
    temp = AutoCorrelation(TransfThetas(i,:));
    Result.ESSTransf(i) = NbIterations/(1+2*sum(temp(2:end)));
    Result.RelESSTransf(i) = Result.ESSTransf(i)/NbIterations*100;
end
for i = 1:size(TransfThetas,1)
    temp = AutoCorrelation(Thetas(i,:));
    Result.ESS(i) = NbIterations/(1+2*sum(temp(2:end)));
    Result.RelESS(i) = Result.ESS(i)/NbIterations*100;
end

if strcmp(Parameters.PMCMC,'Gibbs')                
    Result.GibbsAccRates = Parameters.GibbsAccRates;
    Result.EstimatedSigs = EstimatedSigs; 
end
    
    

try
    if Parameters.SaveSpace
        miness = round(min(Result.ESSTransf));
        NbWeKeep = min(NbIterations,miness*5);
        inds = (randsample(NbIterations,NbWeKeep));
        Result.Thetas = Result.Thetas(:,inds);
        Result.TransfThetas = Result.TransfThetas(:,inds);
        Result.Paths = Result.Paths(inds,:,:);
        Result.LogLiks = Result.LogLiks(inds);
        Result.LogPosts = Result.LogPosts(inds);
    end
end
    