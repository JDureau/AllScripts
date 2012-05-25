function Res = SEIRPollicottInversion(Data,Parameters)

h = Parameters.ComputationTStep;
gamma = Parameters.gamma.Value;
k = Parameters.k.Value;

% Step 1: smooth data.
% TSeries = 0;
% for i = 2:length(Data.Observations)
%     TSeries(i) = TSeries(i-1) + Data.Observations(5,i) - gamma*TSeries(i-1)*Parameters.ComputationTStep*Data.NbComputingSteps(i);
% end
TSeries = Data.ObsForPollicott;

TSeriesSmooth = smoothts(TSeries);
Instants = Data.Instants;
ComputingInstants = Parameters.ComputationTStep:Parameters.ComputationTStep:sum(Data.NbComputingSteps)*Parameters.ComputationTStep;
InterpTSeries = spline(Instants,TSeriesSmooth);
n = length(Instants);
% plot(Instants,TSeries,'o',ComputingInstants,InterpTSeries)


% compute derivatives
tempf = [InterpTSeries(1) InterpTSeries(1) InterpTSeries InterpTSeries(end) InterpTSeries(end)];
f = InterpTSeries;
[breaks,coefs,l,k,d] = unmkpp(f);
df = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
[breaks,coefs,l,k,d] = unmkpp(df);
d2f = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
[breaks,coefs,l,k,d] = unmkpp(d2f);
d3f = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
f = ppval(f,Instants);
df = ppval(df,Instants);
d2f = ppval(d2f,Instants);
d3f = ppval(d3f,Instants);

% tempdf = (tempf(3:end)-tempf(1:end-2))/(2*h);
% tempdf(end-1) = tempdf(end-2);
% tempdf(end) = tempdf(end-1);
% tempdf(2) = tempdf(3);
% tempdf(1) = tempdf(2);
% df = tempdf(2:end-1);
% tempd2f = (tempf(3:end)+tempf(1:end-2)-2*tempf(2:end-1))/(h^2);
% tempd2f(end-1) = tempd2f(end-2);
% tempd2f(end) = tempd2f(end-1);
% tempd2f(2) = tempd2f(3);
% tempd2f(1) = tempd2f(2);
% d2f = tempd2f(2:end-1);
% tempd3f = (tempf(5:end)-tempf(1:end-4)-2*tempf(4:end-1)+2*tempf(2:end-3))/(2*h^3);
% tempd3f(end-1) = tempd3f(end-2);
% tempd3f(end) = tempd3f(end-1);
% tempd3f(2) = tempd3f(3);
% tempd3f(1) = tempd3f(2);
% d3f = tempd3f;


% conditions
Test1 = sum((df+gamma*f) <= 0);
Test2 = sum((d2f + (gamma+k)*df + k*gamma*f)<=0);
if or(Test1,Test2)
    disp('pb with Pollicott')
end

% Step 2: compute p(t)
p = zeros(1,n);
p = (-k*d3f.*f-k*(gamma+k)*d2f.*f-k*k*gamma*df.*f+k*d2f.*df+k*(k+gamma)*df.^2+k*k*gamma*df.*f)./(k*f.*(d2f+(gamma+k)*df+k*gamma*f));
q = (-k*d2f.*f.^2-k*(k+gamma)*df.*f.^2-k^2*gamma*f.^3)./(k*f.*(d2f+(k+gamma)*df+k*gamma*f));

% condition 3
disp(Parameters.betainit.Value)
Beta0 = Parameters.betainit.Value/Parameters.TotalPopulation;
P = cumsum(p.*(h*ones(1,n)));
Test3 = (1/Beta0 + sum(exp(-P).*q.*(h*ones(1,n)))<=0); 
if Test3
    disp('pb with Pollicott')
end

% Step 3: compute beta 
Beta = Parameters.TotalPopulation*ones(1,n)./(exp(P)/Beta0+exp(P).*cumsum(exp(-P).*q.*(h*ones(1,n))));
Beta2 = spline(Instants,Beta,ComputingInstants);

Res.Beta = Beta;
Res.LogLik = sum(log(normpdf(diff(Beta),0,Parameters.SigmaRW.Value*sqrt(Parameters.ComputationTStep))));

% 
% TempVariables = Parameters.InitialState;
% ComputationTStep = Parameters.ComputationTStep;
% TotPop = Parameters.TotalPopulation;
% Variables = TempVariables;
% Rebuiltf = Variables(5);
% Is = Parameters.InitialState(3) ;
% for IndTime = 2:length(Data.Instants)
%     for IndDiscr = 1:Data.NbComputingSteps(IndTime)
%         % Variables
%         beta = Beta2(IndDiscr + sum(Data.NbComputingSteps(1:IndTime-1)))
%         TempVariables(1) = TempVariables(1) + (-beta.*Variables(1).*Variables(3)/TotPop)*ComputationTStep ;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
%         TempVariables(2) = TempVariables(2) + ( beta.*Variables(1).*Variables(3)/TotPop-k*Variables(2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands(:,IndDiscr,2);
%         TempVariables(3) = TempVariables(3) + (-gamma*Variables(3) + k*Variables(2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
%         TempVariables(4) = TempVariables(4) + ( gamma*Variables(3))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
%         TempVariables(5) = ( k*Variables(2))*ComputationTStep ; 
%         Variables = TempVariables;
%         Rebuiltf(end+1) = Rebuiltf(end) + Variables(5);
%     end
%     Is(end+1) = TempVariables(3);
% end
% clf
% plot(Instants,f,'g')
% hold on
% plot(Rebuiltf)
% hold off
% 
% % create points we will use to construct the spline
% x = 0:10;
% y = (x-5).^3+3*3+x+5;
% % points the spline will be plotted at
% xx = linspace(0,10,20);
% % create the spline
% pp = spline(x,y);
% % plot the spline and the initial points used to create it
% figure
% plot(x,y,'o',xx,ppval(pp,xx))
% hold on
% %%
% % extract details from piece-wise polynomial by breaking it apart
% [breaks,coefs,l,k,d] = unmkpp(pp);
% % make the polynomial that describes the derivative
% pp2 = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
% % plot the derivative of the polynomial
% plot(xx,ppval(pp2,xx),'-r')
% %%
% % to calculate 2nd derivative differentiate the 1st derivative
% [breaks,coefs,l,k,d] = unmkpp(pp2);
% pp3 = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
% plot(xx,ppval(pp3,xx),'-g')
% 
% 
% 


