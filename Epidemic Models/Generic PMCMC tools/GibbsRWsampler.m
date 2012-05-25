function Parameters = GibbsRWsampler(Name,Data,Model,Parameters)

ParametersStar = Parameters;
if and(strcmp(Name,'SigmaRW'),not(Parameters.GibbsRepar))
    tmp = diff(Parameters.CurrentPaths(6,1:end-1)).^2/(Parameters.ComputationTStep);
    est = log(sqrt(1/(length(tmp)-1)*sum(tmp)));
    rd = randn(1,1);
    ParametersStar.(Name).TransfValue = est*(1+0.02*Parameters.GibbsEpsils(Parameters.(Name).Index)*rd);
    ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);
    LogQstar = log(normpdf(rd,0,1));
    tmp = (Parameters.(Name).TransfValue-est)/(0.02*Parameters.GibbsEpsils(Parameters.(Name).Index));
    LogQ = log(normpdf(tmp,0,1));
else
%     ParametersStar.(Name).TransfValue = randn(1,1)*Parameters.(Name).GibbsSigma*Parameters.GibbsEpsils(Parameters.(Name).Index) + Parameters.(Name).TransfValue;
%     ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);
%     ParametersStar.(Name).Value = max(0,randn(1,1)*Parameters.(Name).GibbsSigma*Parameters.GibbsEpsils(Parameters.(Name).Index) + Parameters.(Name).Value);
%     ParametersStar = UpdateParsNoTransfToTransf(ParametersStar);
    ParametersStar.(Name).TransfValue = randn(1,1)*Parameters.(Name).GibbsSigma*Parameters.GibbsEpsils(Parameters.(Name).Index) + Parameters.(Name).TransfValue;
    ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);

    if strcmp(Parameters.Model,'DObs')
        if strcmp(Parameters.GibbsReparType,'classic')
            ParametersStar.CurrentPaths = repmat(ParametersStar.betainit.Value + ParametersStar.SigmaRW.Value*ParametersStar.B,2,1);
            Parameters.CurrentPaths = repmat(Parameters.betainit.Value + Parameters.SigmaRW.Value*Parameters.B,2,1);
        elseif strcmp(Parameters.GibbsReparType,'new')
            ParametersStar.CurrentPaths = log(ParametersStar.CurrentMin + (1+ParametersStar.SigmaRW.Value)*ParametersStar.B);
            ParametersStar.betainit.TransfValue = ParametersStar.CurrentPaths(6,1);
            ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);
        end
    elseif strcmp(Parameters.Model,'SEIR')
        if strcmp(Parameters.GibbsReparType,'classic')
            ParametersStar.CurrentPaths(6,:) = ParametersStar.betainit.TransfValue + ParametersStar.SigmaRW.Value*ParametersStar.B;
        elseif strcmp(Parameters.GibbsReparType,'new')
            ParametersStar.CurrentPaths(6,:) = log(ParametersStar.CurrentMin + (1+ParametersStar.SigmaRW.Value)*ParametersStar.B);
            ParametersStar.betainit.TransfValue = ParametersStar.CurrentPaths(6,1);
            ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);
        end
    end
end

tmp =  Model.PathLik(Data,Model,Parameters);
% Res = EstimationSMCsmoothGen(Data,Model,Parameters);

tmp =  Model.PathLik(Data,Model,ParametersStar);
tmp1 = tmp;
LogLikStar = tmp.LogLik;
tmp =  Model.PathLik(Data,Model,Parameters);
LogLik     = tmp.LogLik;

% 
% if strcmp(Parameters.Model,'DObs')
%     plot(tmp1.RecordVariables(1,:),'g')
%     hold on
%     plot(tmp.RecordVariables(1,:))
%     hold off
%     % title(num2str(LogGivenSigmaStar-LogGivenSigma))
%     pause(0.01)   
% elseif strcmp(Parameters.Model,'SEIR')
%     subplot(2,1,1)
%     plot(tmp1.RecordVariables(6,:),'g')
%     hold on
%     plot(tmp.RecordVariables(6,:))
%     plot(ParametersStar.betainit.TransfValue*ones(size(tmp.RecordVariables(6,:))),'--k')
%     hold off
%     subplot(2,1,2)
%     plot(tmp1.RecordVariables(5,:),'g')
%     hold on
%     plot(tmp.RecordVariables(5,:))
%     hold off
%     % title(num2str(LogGivenSigmaStar-LogGivenSigma))
%     pause(0.01)     
% end
LogPriorStar = log(Parameters.(Name).Prior(Name,ParametersStar));
LogPrior     = log(Parameters.(Name).Prior(Name,Parameters));
LogCorrStar = log(Parameters.(Name).CorrFunct(Name,ParametersStar));
LogCorr     = log(Parameters.(Name).CorrFunct(Name,Parameters));

if  and(strcmp(Name,'SigmaRW'),not(Parameters.GibbsRepar))
    LogGivenSigma = sum(log(normpdf(diff(Parameters.CurrentPaths(6,1:end)),0,(Parameters.SigmaRW.Value*sqrt(Parameters.ComputationTStep)))));
    LogGivenSigmaStar = sum(log(normpdf(diff(ParametersStar.CurrentPaths(6,1:end)),0,(ParametersStar.SigmaRW.Value*sqrt(ParametersStar.ComputationTStep)))));
    AccRate = LogLikStar + LogGivenSigmaStar + LogPriorStar - LogCorrStar - LogQstar - (LogLik + LogGivenSigma + LogPrior - LogCorr - LogQ);
else
    AccRate = LogLikStar +  LogPriorStar - LogCorrStar - (LogLik + LogPrior - LogCorr);
end
    
% disp('DeltaLogLik')
% disp(LogLikStar-LogLik)
% disp('DeltaLogPrior')
% disp(LogPriorStar-LogPrior)
% disp('DeltaLogCorr')
% disp(LogCorr-LogCorrStar)

if log(rand(1,1))<AccRate
    disp([Name ': 1' num2str(LogLikStar)])
    Parameters = ParametersStar;
else
    disp([Name ': 0' num2str(LogLikStar)])
end
Parameters.EstimatedSig = sqrt(mean(diff(Parameters.CurrentPaths(1,:)).^2)/Parameters.ComputationTStep);
Parameters.MargLogLik = LogLik;
Parameters.LogPost = LogLik - LogCorr;
Parameters.MargLogLikStar = LogLikStar;
Parameters.ParametersStar = ParametersStar;
