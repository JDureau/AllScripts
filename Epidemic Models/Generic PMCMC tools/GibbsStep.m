function  TempPar = GibbsStep(Data, Model,Parameters, TempPar)
% Simple Gibbs, without repar


% Theta given beta
Names = Parameters.Names.Estimated;

for i = 1:length(Names)
    Parameters.CurrentPaths = TempPar.Paths;
    Parameters.(Names{i}).TransfValue = TempPar.(Names{i}).TransfValue;
    Parameters = UpdateParsTransfToNoTransf(Parameters);
end
IndBetaVar = Parameters.IndBetaVar;
for i = 1:length(Names)
    if strcmp(Names{i},'SigmaRW')
        Beta = Parameters.CurrentPaths(IndBetaVar,:);
        a = 0.1;
        b = 0.1;
        if strcmp(Parameters.Model,'DObs')
            n = length(Data.Instants)-1;
            Parameters.SigmaRW.Value = sqrt(gamrnd((n-1)/2+a,(b+sum(diff(Beta(Data.Instants+1)).^2)/2)^-1)^-1);
        else
            n = length(Beta);
            Parameters.SigmaRW.Value = sqrt(gamrnd((n-1)/2+a,(b+sum(diff(Beta).^2)/2)^-1)^-1);
        end
        Parameters = UpdateParsNoTransfToTransf(Parameters);
        TempPar.(Names{i}).TransfValue = Parameters.(Names{i}).TransfValue;
    else
        Parameters.CurrentPaths = TempPar.Paths;
        ParametersStar = Parameters;
        StarPar = TempPar;      
        StarPar.(Names{i}).TransfValue = TempPar.(Names{i}).TransfValue + randn(1,1)*Parameters.GibbsEpsils(Parameters.(Names{i}).Index)*Parameters.(Names{i}).GibbsSigma;
        ParametersStar.(Names{i}).TransfValue = StarPar.(Names{i}).TransfValue;
        ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);

        if or(strcmp(Names{i},'betainit'),strcmp(Names{i},'SigmaRW'))
            if strcmp(Parameters.Model,'SEIR')
                ParametersStar.CurrentPaths(IndBetaVar,:) = ParametersStar.CurrentPaths(IndBetaVar,:) + ParametersStar.betainit.TransfValue -Parameters.betainit.TransfValue ;
                StarPar.Paths = ParametersStar.CurrentPaths;
            elseif strcmp(Parameters.Model,'DObs')
                BetaStar = ParametersStar.betainit.Value + ((Parameters.CurrentPaths(IndBetaVar,:) - Parameters.betainit.Value)/Parameters.SigmaRW.Value)*ParametersStar.SigmaRW.Value;
                ParametersStar.CurrentPaths(IndBetaVar,:) = BetaStar;
                StarPar.Paths = ParametersStar.CurrentPaths;
            end
        end

        Temp = Model.PathLik(Data, Model,Parameters);
        TempStar =  Model.PathLik(Data, Model,ParametersStar);

        LogPrior = log(Parameters.(Names{i}).Prior(Names{i},Parameters));
        LogPriorStar = log(Parameters.(Names{i}).Prior(Names{i},ParametersStar));
        LogCorr = log(Parameters.(Names{i}).CorrFunct(Names{i},Parameters));
        LogCorrStar = log(Parameters.(Names{i}).CorrFunct(Names{i},ParametersStar));
        AccRate = TempStar.LogLik + LogPriorStar - LogCorrStar - (Temp.LogLik + LogPrior - LogCorr);
        if log(rand(1,1))<AccRate
            TempPar = StarPar;
            Parameters = ParametersStar;
        end
    end
end

% beta given Theta
Parameters.ForceTraj = 1;
Parameters.ForcedTraj = squeeze(TempPar.Paths);
Temp = EstimationSMCsmoothGen(Data,Model,Parameters);
ind = ceil(rand(1,1)*Parameters.NbParticules);
TempPar.Paths = squeeze(Temp.CompletePaths(ind,:,:));
TempPar.LogLik = Temp.LogLik;
TempPar.Coalescence = Temp.Coalescence;

