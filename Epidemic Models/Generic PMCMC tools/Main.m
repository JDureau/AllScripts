% Tests.
k = 1;
theta = sqrt(1/k);
epss = 0.1:0.5:6;
NbIts = 30000;

AccRates = [];
Samples = zeros(length(epss),NbIts);
for IndEps = 1:length(epss)
    IndEps
    epsil = epss(IndEps);
    Temp = k*theta;
    Samples(IndEps,1) = Temp;
    LogLikTemp = log(gampdf(Temp,k,theta));
    AccRate = 1;
    for i = 2:NbIts 
        Pot = Temp + epsil*randn(1,1);
        LogLikPot = log(gampdf(Pot,k,theta));
        AccRatio = LogLikPot - LogLikTemp;
        if AccRatio>log(rand(1,1))
            Temp = Pot;
            LogLikTemp = LogLikPot;
            AccRate = (AccRate*(i-1) + 1)/i;
        else
            AccRate = (AccRate*(i-1))/i;
        end
        Samples(IndEps,i) = Temp;
    end
    AccRates(IndEps) = AccRate;
end

ESSs = [];
for IndEps = 1:length(epss)
    temp = AutoCorrelation(Samples(IndEps,:),150);
    ESSs(IndEps) = NbIts/(1+2*sum(temp));
end
tests = [];
pvalues = [];
for IndEps = 1:length(epss)
        CDF = [Samples(IndEps,:)' gamcdf(Samples(IndEps,:),k,theta)'];
        [h,p] = kstest(Samples(IndEps,:),CDF);
        tests(IndEps) = 1 - h;
        pvalues(IndEps) = p;
end
subplot(3,1,1)
plot(epss,AccRates)
subplot(3,1,2)
plot(epss,ESSs)
subplot(3,1,3)
plot(epss,pvalues)

TrueSamples = gamrnd(k,theta,20000,1);
[fT,xT] = ksdensity(TrueSamples);
for i = 1:length(epss)
    subplot(2,1,1)
    plot(epss,pvalues)
%     plot(epss,tests)
    hold on
    plot(epss(i),pvalues(i),'og')
%     plot(epss(i),tests(i),'og')
    hold off
    subplot(2,1,2)
    plot(xT,fT,'g')
    hold on
    [f,x] = ksdensity(Samples(i,:));
    plot(x,f)
    hold off
    pause()
end

% Does the 0.24 always induce good sampling? 
Samples = [];
Pars = [];
for IndPars = 1:20
    k = rand(1,1)*10;
    theta = sqrt(1/k);
    
    EpsStd = 1;
    CoolingParameter = 0.98;
    TempEps = 1;
    TempAccRate = 1;
    NbIts = 10000;
    IndEps = 0;
    Condition = 0;
    TempSamples = [];
    while not(Condition)

        IndEps = IndEps + 1;

        % Sample PotEps
        pos = 0;
        while not(pos)
            PotEps = TempEps + EpsStd * CoolingParameter^IndEps*randn(1,1);
            if PotEps > 0
                pos = 1;
            end
        end

        % ComputeAccRate
        Temp = k*theta;
        TempSamples(IndEps,1) = Temp;
        LogLikTemp = log(gampdf(Temp,k,theta));
        AccRate = 1;
        for i = 2:NbIts 
            Pot = Temp + PotEps*randn(1,1);
            LogLikPot = log(gampdf(Pot,k,theta));
            AccRatio = LogLikPot - LogLikTemp;
            if AccRatio>log(rand(1,1))
                Temp = Pot;
                LogLikTemp = LogLikPot;
                AccRate = (AccRate*(i-1) + 1)/i;
            else
                AccRate = (AccRate*(i-1))/i;
            end
            TempSamples(IndEps,i) = Temp;
        end
        PotAccRate = AccRate;
        AccRates(IndEps) = AccRate;

        % Accept/Reject
        if abs(PotAccRate - 0.23)< abs(TempAccRate - 0.23)
            TempEps = PotEps;
            TempAccRate = PotAccRate;
        end
        if and(TempAccRate>0.20,TempAccRate<0.26) 
            Condition = 1;
        end
        disp(['AccRate : ' num2str(TempAccRate) ', Eps = ' num2str(TempEps)])
    end

    Pars(IndPars,:) = [k theta TempEps];
    Samples(IndPars,:) = TempSamples(end,:);
end

for i = 1:size(Samples,1)
    plot(Samples(i,:))
    pause()
end
    
    

