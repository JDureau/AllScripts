function Res = RunMCMC(InitValue,Parameters,NbIts)


% Initialize Parameters
Dens = Parameters.Dens;
Parameters.Mus = Dens.mu;
Parameters.Sigmas = Dens.Sigma;
Parameters.Nc = Dens.NComponents;

LogRatioFun = Parameters.LogRatioFun;
SampleFun = Parameters.SampleFun;

TempValue = InitValue';
TempValue = Parameters.SampleFun(TempValue,Parameters);

logliks = [];
vals = [];
Potvals = [];
NbAcc = 0;
for i = 1:NbIts
    % sample new
    PotValue = Parameters.SampleFun(TempValue,Parameters);
    LogRatio = Parameters.LogRatioFun(TempValue,PotValue,Parameters);
    rd = rand(1,1);
    if log(rd)<LogRatio
        TempValue = PotValue;
        NbAcc = NbAcc+1;
%         ezcontour(@(x,y)pdf(Parameters.OptDens,[x y]),[-20 20],[-40 20],500);
%         hold on
%         plot(PotValue(1),PotValue(2),'og', 'MarkerFaceColor','g')
%         hold off
%          pause(0.001)
    else
%         ezcontour(@(x,y)pdf(Parameters.OptDens,[x y]),[-20 20],[-40 20],500);
%         hold on
%         plot(PotValue(1),PotValue(2),'or', 'MarkerFaceColor','r')
%         hold off
%         pause(0.001)
    end
%     try
%         logliks(i) = log(Parameters.f(TempValue',Parameters));
%     catch
%         disp('pb');
%     end
    vals(:,i) = TempValue;
    Potvals(:,i) = PotValue;
end

% 
% try
%     subplot(2,1,1)
%     hist(vals(1,:))
%     subplot(2,1,2)
%     hist(vals(2,:))
% catch
%     disp('pb')
% end

for i = 1:Parameters.Dim
    subplot(Parameters.Dim,1,i)
    hist(vals(i,:))
    autocor = autocorrelation(vals(i,:),100);
    Res.autocor(i,:) = autocor;
    Res.ESS(i) = NbIts/(1+2*sum(autocor));
    Res.RelESS(i) = Res.ESS(i)/NbIts*100;
    title([num2str(Res.ESS(i)/NbIts*100,3) '%'])
end

Res.LogLiks = logliks;
Res.Vals = vals;
Res.PotVals = Potvals;
Res.AccRate = NbAcc/NbIts;
Res.Parameters = Parameters;