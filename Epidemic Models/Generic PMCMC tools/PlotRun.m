function [] = PlotRun(Res,Parameters)

maxind = max(Res.CompInds);

figure(1)
if Parameters.EstimatedVariables.Beta
    subplot(Parameters.NbParsEstimated*5+1,1,(Parameters.BetaIndex-1)*5+1:(Parameters.BetaIndex)*5)
    plot(Res.Pars(:,Parameters.BetaIndex),'k')
    hold on
    plot(Res.ProposedPars(:,Parameters.BetaIndex),'g')
    plot(Res.Mus(:,Parameters.BetaIndex))
    plot(Res.Mus(:,Parameters.BetaIndex) + sqrt(Res.SamplCovs(:,Parameters.BetaIndex,Parameters.BetaIndex)),'r')
    plot(Res.Mus(:,Parameters.BetaIndex) - sqrt(Res.SamplCovs(:,Parameters.BetaIndex,Parameters.BetaIndex)),'r')
    hold off
    ylim([min(Res.Pars(:,Parameters.BetaIndex))-2*std(Res.Pars(:,Parameters.BetaIndex)) max(Res.Pars(:,Parameters.BetaIndex))+2*std(Res.Pars(:,Parameters.BetaIndex))])
    title('Log plots')
end

if Parameters.EstimatedVariables.Gamma
    subplot(Parameters.NbParsEstimated*5+1,1,(Parameters.GammaIndex-1)*5+1:(Parameters.GammaIndex)*5)
    plot(Res.Pars(:,Parameters.GammaIndex),'k')
    hold on
    plot(Res.ProposedPars(:,Parameters.GammaIndex),'g')
    plot(Res.Mus(:,Parameters.GammaIndex))
    plot(Res.Mus(:,Parameters.GammaIndex) + sqrt(Res.SamplCovs(:,Parameters.GammaIndex,Parameters.GammaIndex)),'r')
    plot(Res.Mus(:,Parameters.GammaIndex) - sqrt(Res.SamplCovs(:,Parameters.GammaIndex,Parameters.GammaIndex)),'r')
    hold off
    ylim([min(Res.Pars(:,Parameters.GammaIndex))-2*std(Res.Pars(:,Parameters.GammaIndex)) max(Res.Pars(:,Parameters.GammaIndex))+2*std(Res.Pars(:,Parameters.GammaIndex))])
end

if Parameters.EstimatedVariables.SigmaRW
    subplot(Parameters.NbParsEstimated*5+1,1,(Parameters.SigmaRWIndex-1)*5+1:(Parameters.SigmaRWIndex)*5)
    plot(Res.Pars(:,Parameters.SigmaRWIndex),'k')
    hold on
    plot(Res.ProposedPars(:,Parameters.SigmaRWIndex),'g')
    plot(Res.Mus(:,Parameters.SigmaRWIndex))
    plot(Res.Mus(:,Parameters.SigmaRWIndex) + sqrt(Res.SamplCovs(:,Parameters.SigmaRWIndex,Parameters.SigmaRWIndex)),'r')
    plot(Res.Mus(:,Parameters.SigmaRWIndex) - sqrt(Res.SamplCovs(:,Parameters.SigmaRWIndex,Parameters.SigmaRWIndex)),'r')
    hold off
    ylim([min(Res.Pars(:,Parameters.SigmaRWIndex))-2*std(Res.Pars(:,Parameters.SigmaRWIndex)) max(Res.Pars(:,Parameters.SigmaRWIndex))+2*std(Res.Pars(:,Parameters.SigmaRWIndex))])
end


subplot(Parameters.NbParsEstimated*5+1,1,Parameters.NbParsEstimated*5+1)
plot(1,0,'o','MarkerFaceColor',[Res.CompInds(1)/maxind,0,0])
hold on
for i = 1:size(Res.Mus,1)
    plot(i,0,'o','MarkerFaceColor',[Res.CompInds(i)/maxind,0,0],'MarkerEdgeColor',[Res.CompInds(i)/maxind,0,0])
end
hold off
set(gca,'YTickLabel',{''})
set(gca,'XTickLabel',{''})


figure(2)
if Parameters.EstimatedVariables.Beta
    subplot(Parameters.NbParsEstimated*5+1,1,(Parameters.BetaIndex-1)*5+1:(Parameters.BetaIndex)*5)
    plot(exp(Res.Pars(:,Parameters.BetaIndex)),'k')
    hold on
    plot(exp(Res.ProposedPars(:,Parameters.BetaIndex)),'g')
    plot(exp(Res.Mus(:,Parameters.BetaIndex)))
    plot(exp(Res.Mus(:,Parameters.BetaIndex) + sqrt(Res.SamplCovs(:,Parameters.BetaIndex,Parameters.BetaIndex))),'r')
    plot(exp(Res.Mus(:,Parameters.BetaIndex) - sqrt(Res.SamplCovs(:,Parameters.BetaIndex,Parameters.BetaIndex))),'r')
    hold off
    ylim([min(exp(Res.Pars(:,Parameters.BetaIndex)))-2*std(exp(Res.Pars(:,Parameters.BetaIndex))) max(exp(Res.Pars(:,Parameters.BetaIndex)))+2*std(exp(Res.Pars(:,Parameters.BetaIndex)))])
    title('Natural scale plots')
end
    
if Parameters.EstimatedVariables.Gamma
    subplot(Parameters.NbParsEstimated*5+1,1,(Parameters.GammaIndex-1)*5+1:(Parameters.GammaIndex)*5)
    plot(exp(Res.Pars(:,Parameters.GammaIndex)),'k')
    hold on
    plot(exp(Res.ProposedPars(:,Parameters.GammaIndex)),'g')
    plot(exp(Res.Mus(:,Parameters.GammaIndex)))
    plot(exp(Res.Mus(:,Parameters.GammaIndex) + sqrt(Res.SamplCovs(:,Parameters.GammaIndex,Parameters.GammaIndex))),'r')
    plot(exp(Res.Mus(:,Parameters.GammaIndex) - sqrt(Res.SamplCovs(:,Parameters.GammaIndex,Parameters.GammaIndex))),'r')
    hold off
    ylim([min(exp(Res.Pars(:,Parameters.GammaIndex)))-2*std(exp(Res.Pars(:,Parameters.GammaIndex))) max(exp(Res.Pars(:,Parameters.GammaIndex)))+2*std(exp(Res.Pars(:,Parameters.GammaIndex)))])
end

if Parameters.EstimatedVariables.SigmaRW
    subplot(Parameters.NbParsEstimated*5+1,1,(Parameters.SigmaRWIndex-1)*5+1:(Parameters.SigmaRWIndex)*5)
    plot(exp(Res.Pars(:,Parameters.SigmaRWIndex)),'k')
    hold on
    plot(exp(Res.ProposedPars(:,Parameters.SigmaRWIndex)),'g')
    plot(exp(Res.Mus(:,Parameters.SigmaRWIndex)))
    plot(exp(Res.Mus(:,Parameters.SigmaRWIndex) + sqrt(Res.SamplCovs(:,Parameters.SigmaRWIndex,Parameters.SigmaRWIndex))),'r')
    plot(exp(Res.Mus(:,Parameters.SigmaRWIndex) - sqrt(Res.SamplCovs(:,Parameters.SigmaRWIndex,Parameters.SigmaRWIndex))),'r')
    hold off
    ylim([min(exp(Res.Pars(:,Parameters.SigmaRWIndex)))-2*std(exp(Res.Pars(:,Parameters.SigmaRWIndex))) max(exp(Res.Pars(:,Parameters.SigmaRWIndex)))+2*std(exp(Res.Pars(:,Parameters.SigmaRWIndex)))])
end


subplot(Parameters.NbParsEstimated*5+1,1,Parameters.NbParsEstimated*5+1)
plot(1,0,'o','MarkerFaceColor',[Res.CompInds(1)/maxind,0,0])
hold on
for i = 1:size(Res.Mus,1)
    plot(i,0,'o','MarkerFaceColor',[Res.CompInds(i)/maxind,0,0],'MarkerEdgeColor',[Res.CompInds(i)/maxind,0,0])
end
hold off
set(gca,'YTickLabel',{''})
set(gca,'XTickLabel',{''})
