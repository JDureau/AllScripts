function  TempPar = GibbsReparStep(Data, Model,Parameters, TempPar)

% Theta given beta
Names = Parameters.Names.Estimated;

for i = 1:length(Names)
    Parameters.CurrentPaths = TempPar.Paths;
    Parameters.(Names{i}).TransfValue = TempPar.(Names{i}).TransfValue;
    Parameters = UpdateParsTransfToNoTransf(Parameters);
end

IndBetaVar = Parameters.IndBetaVar;
for i = 1:length(Names)
    Parameters.CurrentPaths = TempPar.Paths;
    
    ParametersStar = Parameters;
    StarPar = TempPar;      
    StarPar.(Names{i}).TransfValue = TempPar.(Names{i}).TransfValue + randn(1,1)*Parameters.GibbsEpsils(Parameters.(Names{i}).Index)*Parameters.(Names{i}).GibbsSigma;
    ParametersStar.(Names{i}).TransfValue = StarPar.(Names{i}).TransfValue;
    ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);
    TellParsValues(ParametersStar);
    if or(strcmp(Names{i},'betainit'),strcmp(Names{i},'SigmaRW'))
        if strcmp(Parameters.Model,'SEIR')
            % take out trqnsf when working on DOBS
            ParametersStar.CurrentPaths(IndBetaVar,:) = ParametersStar.betainit.TransfValue + ((Parameters.CurrentPaths(IndBetaVar,:) - Parameters.betainit.TransfValue)/Parameters.SigmaRW.Value)*ParametersStar.SigmaRW.Value;
            StarPar.Paths = ParametersStar.CurrentPaths;
        elseif strcmp(Parameters.Model,'DObs')
            BetaStar = ParametersStar.betainit.Value + ((Parameters.CurrentPaths(IndBetaVar,:) - Parameters.betainit.Value)/Parameters.SigmaRW.Value)*ParametersStar.SigmaRW.Value;
            ParametersStar.CurrentPaths(IndBetaVar,:) = BetaStar;
            StarPar.Paths = ParametersStar.CurrentPaths;
        end
    end

    Temp = Model.PathLik(Data, Model,Parameters);
    TempStar =  Model.PathLik(Data, Model,ParametersStar);

%         LogLikStar = sum(log(normpdf(BetaStar,Data.Observations,Parameters.SigmaObs)));
%         LogLik = sum(log(normpdf(Beta,Data.Observations,Parameters.SigmaObs)));

    LogPrior = log(Parameters.(Names{i}).Prior(Names{i},Parameters));
    LogPriorStar = log(ParametersStar.(Names{i}).Prior(Names{i},ParametersStar));
    LogCorr = log(Parameters.(Names{i}).CorrFunct(Names{i},Parameters));
    LogCorrStar = log(ParametersStar.(Names{i}).CorrFunct(Names{i},ParametersStar));
    AccRate = TempStar.LogLik + LogPriorStar - LogCorrStar - (Temp.LogLik + LogPrior - LogCorr);
%         AccRate = LogLikStar + LogPriorStar - LogCorrStar - (LogLik + LogPrior - LogCorr);
    if log(rand(1,1))<AccRate
        TempPar = StarPar;
        Parameters = ParametersStar;
    end

end

% beta given Theta

% % beta given sigma
% Beta = Parameters.CurrentPaths(IndBetaVar,:);
% n = size(TempPar.Paths,2);
% BetasSamples = ParametersStar.betainit.Value * ones(Parameters.NbParticules,n);
% ind = 1;
% for IndStep = 2:n
%     for IndDiscr = 1:Data.NbComputingSteps(2)
%         ind = ind + 1;
%         BetasSamples(:,ind) = BetasSamples(:,ind-1)+randn(Parameters.NbParticules,1)*Parameters.SigmaRW.Value;
%         BetasSamples(1,ind) = Beta(ind);
%     end
%     Weigths = normpdf(BetasSamples(:,ind),Data.Observations(IndStep),Parameters.SigmaObs);
% 
% %     BetasSamples(:,IndStep) = BetasSamples(:,IndStep-1)+randn(Parameters.NbParticules,1)*Parameters.SigmaRW.Value;
% %     BetasSamples(1,IndStep) = Beta(IndStep);
% %     Weigths = normpdf(BetasSamples(:,IndStep),Data.Observations(IndStep),Parameters.SigmaObs);
% 
%     Weigths = Weigths/sum(Weigths);
%     u = rand(1,1)/Parameters.NbParticules;
%     s = 0;
%     KeptInds = [];
%     resind = 1;
%     for ipart = 1:Parameters.NbParticules
%         k = 0;
%         s = s+Weigths(ipart);
%         while s>u
%             k=k+1;
%             u = u+1/Parameters.NbParticules;
%             KeptInds(resind) = ipart;
%             resind = resind+1;
%         end
%     end
%     KeptInds(1) = 1; 
%     BetasSamples = BetasSamples(KeptInds,:);
% end
% RandInd = ceil(rand(1,1)*Parameters.NbParticules);
% Beta = BetasSamples(RandInd,:);
% TempPar.Paths = [Beta;Beta];
% TempPar.LogLik = 10;
% TempPar.Coalescence = 0;


Parameters.ForceTraj = 1;
Parameters.ForcedTraj = squeeze(TempPar.Paths);
Temp = EstimationSMCsmoothGen(Data,Model,Parameters);
ind = ceil(rand(1,1)*Parameters.NbParticules);
TempPar.Paths = squeeze(Temp.CompletePaths(ind,:,:));
TempPar.LogLik = Temp.LogLik;
disp(TempPar.LogLik)
TempPar.Coalescence = Temp.Coalescence;


