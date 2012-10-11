function Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,NbIts,IndModel)


mu = Parameters.mu;
Epsil = 1;
NamesEst = Parameters.NamesEst;
indpars = [];
dim = length(NamesEst);
% for i = 1:length(NamesEst)
%     indpars(i) = Parameters.(NamesEst{i}).Index;
% end
indpars = 1:length(NamesEst);
% die
x = [];
for i = 1:length(NamesEst)
    x(i) = Parameters.(NamesEst{i}).TransfValue;
end
xnostransf = zeros(size(x));
for i = 1:length(NamesEst)
    xnotransf(indpars(i)) =   Parameters.(NamesEst{i}).Value;
end
if IndModel == 1
    TempSim = SimulateBertallanfyDetermEpid(Data,Parameters,HIVModel);
elseif IndModel == 2
    TempSim = SimulateSigmoidDetermEpid(Data,Parameters,HIVModel);
end    
tmpLogLik = 0;
tmpLogPrior = 0;
tmpLogCorr = 0;
for i = 2:length(Data.ObservedVariables);
%     tmpLogLik = tmpLogLik + max(-700,not(TempSim.Crash)*log(normpdf(Data.Observations(Data.ObservedVariables(i),i),TempSim.Observations(Data.ObservedVariables(i),i),sqrt(Data.Observations(Data.ObservedVariables(i),i)*(100-Data.Observations(Data.ObservedVariables(i),i))/400))));
    tmpLogLik = tmpLogLik + max(-700,log(binopdf(round(Parameters.NbSamples(i-1)*Data.Observations(Data.ObservedVariables(:,i),i)/100),Parameters.NbSamples(i-1),TempSim.Observations(Data.ObservedVariables(i),i)/100)));

end
for i = 1:length(NamesEst)
    tmp = Parameters.(NamesEst{i}).Prior(NamesEst{i},Parameters);
    tmpLogPrior = tmpLogPrior + log(tmp);
    tmpLogCorr = tmpLogCorr + log(Parameters.(NamesEst{i}).CorrFunct(NamesEst{i},Parameters));
end
LogPost = -Inf;
LogLik = tmpLogLik;
TempSim.LogLik = tmpLogLik;
TempSim.LogPrior = tmpLogPrior;
TempSim.LogCorr = tmpLogCorr;


Accepted = [];
TransfThetas = [];
Thetas = [];
LogPosts = [];
LogLiks = [];
AdaptC = 0.99;
tis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
Parameters.NbTSteps = length(tis);
Paths = zeros(NbIts,Parameters.NbTSteps);

CovInit = Cov;
CptNotAccepted = 0;
for  j = 1:NbIts
    tmpCov = CovInit;
%     if Parameters.AdMet
%         if j > 500
%             if rand(1,1)>Parameters.AdMetBeta
%                 try
%                     if not(sum(eig(cov(TransfThetas'))<=eps))
%                         tmpCov = 2.38^2/dim*cov(TransfThetas');
% %                         disp(eig(cov(TransfThetas'))')
%                     end
%                 end
%             end
%         end
%     end

    probother = 0;
    covcoeff  = 10; 
%     if rand(1,1)>probother
%         xstar = mvnrnd(x' + Epsil^2/2*(mu-x)',Epsil^2*2.38^2/dim*tmpCov);
%     else
%         xstar = mvnrnd(x' + Epsil^2/2*(mu-x)',Epsil^2*covcoeff*2.38^2/dim*tmpCov);
%     end
%     LogqStarTemp = log(mvnpdf(x',xstar' + Epsil^2/2*(mu-xstar)',Epsil^2*2.38^2/dim*tmpCov));
%     LogqTempStar = log(mvnpdf(xstar',x' + Epsil^2/2*(mu-x)'    ,Epsil^2*2.38^2/dim*tmpCov));
%     LogqStarTemp = log((1-probother)*mvnpdf(x',xstar' + Epsil^2/2*(mu-xstar)',Epsil^2*2.38^2/dim*tmpCov)+probother*mvnpdf(x',xstar' + Epsil^2/2*(mu-xstar)',Epsil^2*covcoeff*2.38^2/dim*tmpCov));
%     LogqTempStar = log((1-probother)*mvnpdf(xstar',x' + Epsil^2/2*(mu-x)'    ,Epsil^2*2.38^2/dim*tmpCov)+probother*mvnpdf(xstar',x' + Epsil^2/2*(mu-x)'    ,Epsil^2*covcoeff*2.38^2/dim*tmpCov));
    xstar = mvnrnd(x' ,Epsil^2*2.38^2/dim*tmpCov);
%     LogqStarTemp = log(mvnpdf(x',xstar' ,Epsil^2*tmpCov));
%     LogqTempStar = log(mvnpdf(xstar',x'    ,Epsil^2*tmpCov));
    
%     if CptNotAccepted > 500
%         xstar = mu;
%         LogqStarTemp = Inf;
%         LogqTempStar = 0;
%         disp('Forced move')
%     end
    
    
    
    
%     xstar = mvnrnd(x,Epsil^2*tmpCov); 
    ParametersStar = Parameters;
    for i = 1:length(NamesEst)        
        ParametersStar.(NamesEst{i}).TransfValue = (xstar(indpars(i)));
    end
    ParametersStar = UpdateParsTransfToNoTransf(ParametersStar);
    if IndModel == 1
        TempSimStar = SimulateBertallanfyDetermEpid(Data,ParametersStar,HIVModel);
    elseif IndModel == 2
        TempSimStar = SimulateSigmoidDetermEpid(Data,ParametersStar,HIVModel);
    end
    
    
  
    
    
    
    tmpLogLik = 0;
    tmpLogPrior = 0;
    tmpLogCorr = 0;
    for i = 2:length(Data.ObservedVariables)
%         tmpLogLik = tmpLogLik + max(-700,log(normpdf(Data.Observations(Data.ObservedVariables(i),i),TempSimStar.Observations(Data.ObservedVariables(i),i),sqrt(Data.Observations(Data.ObservedVariables(i),i)*(100-Data.Observations(Data.ObservedVariables(i),i))/400))));
        tmpLogLik = tmpLogLik + max(-700,log(binopdf(round(Parameters.NbSamples(i-1)*Data.Observations(Data.ObservedVariables(:,i),i)/100),Parameters.NbSamples(i-1),TempSim.Observations(Data.ObservedVariables(i),i)/100)));

    end
    for i = 1:length(NamesEst)
        tmp = Parameters.(NamesEst{i}).Prior(NamesEst{i},Parameters);
        tmpLogPrior = tmpLogPrior + log(tmp);
        tmpLogCorr = tmpLogCorr + log(ParametersStar.(NamesEst{i}).CorrFunct(NamesEst{i},ParametersStar));
    end
    LogPostStar = tmpLogPrior + tmpLogLik - tmpLogCorr ;
    LogLikStar = tmpLogLik;
    TempSimStar.LogLik = tmpLogLik;
    TempSimStar.LogPrior = tmpLogPrior;
    TempSimStar.LogCorr = tmpLogCorr;
    
%     disp([TempSimStar.LogLik - TempSim.LogLik  TempSimStar.LogPrior - TempSim.LogPrior TempSimStar.LogCorr-TempSim.LogCorr])
%     disp([LogPostStar - LogPost  LogqStarTemp - LogqTempStar
%     LogPostStar+LogqStarTemp-LogPost-LogqTempStar])
   
%     if log(rand(1,1))<LogPostStar+LogqStarTemp-LogPost-LogqTempStar
    if log(rand(1,1))<LogPostStar-LogPost
        x = xstar;
        Parameters = ParametersStar;
        LogPost = LogPostStar;
        LogLik = LogLikStar;
        TempSim = TempSimStar;
        Accepted(j) = 1;
        xnostransf = zeros(size(x));
        for i = 1:length(NamesEst)
            xnotransf(indpars(i)) =   ParametersStar.(NamesEst{i}).Value;
        end
        if TempSim.Crash
            'pb'
        end
        CptNotAccepted = 0;
    else
        Accepted(j) = 0;
        CptNotAccepted = CptNotAccepted+1; 
    end    
    TransfThetas(:,j) = x;
    Thetas(:,j) = xnotransf;
    LogPosts(j) = LogPost;
    LogLiks(j) = LogLik;
    Paths(j,:) = TempSim.BuiltTraj(:,9)';
    if j >30
%         if rand(1,1)<0.1
            Epsil = exp(log(Epsil) + AdaptC^j*(mean(Accepted)-0.23));
%         end
    end
    if rand(1,1)<0.05
        disp([num2str(j) '    ' num2str(mean(Accepted)) ])
    end
    tin = Parameters.BRtinfl.Value;
%     clf
%     plot([tin tin]*10, [0 1])
%     hold on
%     plot(TempSim.Fts)
%     hold off
%     pause(0.01)
end
mean(Accepted)  




Parameters.TypeWork = 'Boston Examples';

Res.Data = Data;
Res.Thetas = Thetas;
Res.TransfThetas = TransfThetas;
Res.Accepted = Accepted;
Res.Parameters = Parameters;
Res.Paths = Paths;
Res.LogPosts = LogPosts;
Res.LogLiks = LogLiks;

   

for i = 1:size(TransfThetas,1)
    temp = AutoCorrelation(TransfThetas(i,:),100);
    Res.ESSTransf(i) = NbIts/(1+2*sum(temp(2:end)));
    Res.RelESSTransf(i) = Res.ESSTransf(i)/NbIts*100;
end
for i = 1:size(TransfThetas,1)
    temp = AutoCorrelation(Thetas(i,:),100);
    Res.ESS(i) = NbIts/(1+2*sum(temp(2:end)));
    Res.RelESS(i) = Res.ESS(i)/NbIts*100;
end


try
    if Parameters.SaveSpace
        miness = min(Res.ESSTransf);
        NbWeKeep = min(NbIts,round(miness*5));
        disp(['NbWeKeep:' num2str(NbWeKeep)])
        inds = randsample(NbIts,NbWeKeep);
        length(inds)
        Res.Thetas = Res.Thetas(:,inds);
        Res.TransfThetas = Res.TransfThetas(:,inds);
        Res.Paths = Res.Paths(inds,:,:);
        Res.LogLiks = Res.LogLiks(inds);
        Res.LogPosts = Res.LogPosts(inds);
        size(Res.Thetas)
        Res
    end
end