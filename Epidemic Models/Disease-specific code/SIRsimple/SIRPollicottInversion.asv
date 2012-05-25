function Res = SIRPollicottInversion(Data,Parameters)

h = Parameters.ComputationTStep;
gamma = Parameters.gamma.Value;

% Step 1: smooth data.
TSeries = Data.Observations(1,:);
% for i = 2:length(Data.Observations)
%     TSeries(i) = TSeries(i-1) + Data.Observations(5,i) - gamma*TSeries(i-1)*Parameters.ComputationTStep*Data.NbComputingSteps(i);
% end

Instants = Data.Instants;
ComputingInstants = Parameters.ComputationTStep:Parameters.ComputationTStep:sum(Data.NbComputingSteps)*Parameters.ComputationTStep;
InterpTSeries = spline(Instants,TSeries,ComputingInstants);
n = length(InterpTSeries);
% plot(Instants,TSeries,'o',ComputingInstants,InterpTSeries)



% compute derivatives
tempf = [InterpTSeries(1) InterpTSeries(1) InterpTSeries InterpTSeries(end) InterpTSeries(end)];
f = InterpTSeries;
tempdf = (tempf(3:end)-tempf(1:end-2))/(2*h);
tempdf(end-1) = tempdf(end-2);
tempdf(end) = tempdf(end-1);
tempdf(2) = tempdf(3);
tempdf(1) = tempdf(2);
df = tempdf(2:end-1);
tempd2f = (tempf(3:end)+tempf(1:end-2)-2*tempf(2:end-1))/(h^2);
tempd2f(end-1) = tempd2f(end-2);
tempd2f(end) = tempd2f(end-1);
tempd2f(2) = tempd2f(3);
tempd2f(1) = tempd2f(2);
d2f = tempd2f(2:end-1);
tempd3f = (tempf(5:end)-tempf(1:end-4)-2*tempf(4:end-1)+2*tempf(2:end-3))/(2*h^3);
tempd3f(end-1) = tempd3f(end-2);
tempd3f(end) = tempd3f(end-1);
tempd3f(2) = tempd3f(3);
tempd3f(1) = tempd3f(2);
d3f = tempd3f;

Test1 = sum((df+gamma*f) <= 0);
if or(Test1,0)
    disp('pb with Pollicott')
end

% Step 2: compute p(t)
p = zeros(1,n);
p = (d2f.*f-df.^2)./(f.*(df+gamma*f));


% condition 3
disp(Parameters.betainit.Value)
Beta0 = 1.9;
P = cumsum(p.*(h*ones(1,n)));
Test3 = (1/Beta0 - sum(exp(P).*f.*(h*ones(1,n)))<=0); 
if Test3
    disp('pb with Pollicott (beta0)')
end

% Step 3: compute beta 
BetaRebuilt = ones(1,n)./(exp(-P)/Beta0-exp(-P).*cumsum(exp(P).*f.*(h*ones(1,n))));

plot(BetaRebuilt)
pause(0.01)

Res.Beta = BetaRebuilt;
Res.LogLik = sum(log(normpdf(diff(BetaRebuilt),0,Parameters.SigmaRW.Value*sqrt(Parameters.ComputationTStep))));














