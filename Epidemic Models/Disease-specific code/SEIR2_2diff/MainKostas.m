% Main Kostas

%% 1: generate data

Tau=0.05; %noise st. deviation
T=30; %end time
npoints=30; %number of observed points apart from the initial
obstep=T/npoints;
M=9; %number of imputed points between observations
step=obstep/(M+1); %euler step

%epidemic parameters
gamma = 0.94963^(-1);
k = 1.5647^(-1);
TotPop = 100000;
I0=1;
%sigmoid for beta
RealAmpl = 0.9;
RealBaseline = 1.0;
RealInflpt = 0.4*T/step;
RealSteepness = T/20;
xis = 0:step:T;
beta =  exp(RealBaseline-RealAmpl*Sigmoid((xis-RealInflpt*step)/RealSteepness));

[S,E,I,R,X]=ODESEIR(beta,gamma,k,I0,TotPop,npoints,M,step);
%X denote the latent true value of the observations.
Y=zeros(1,npoints);
Rands=randn(1,npoints);
for i=1:npoints
    Y(i)= exp(log(X(i))+Tau*Rands(i));%normal in the log scale
end

Parameters = struct();
Parameters.ComputationTStep = step;
Data.Observations = zeros(7,31);
Data.Observations(5,2:31) = Y;

Parameters.DiffusionType  = 'Add';
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:7]';
Parameters.TotalPopulation = TotPop;
Data.Instants = [1:size(Data.Observations,2)]/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];
Parameters.SigmaObs.Value = Tau;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';

Parameters.NbVariables = 7;
Parameters.SigmaRW.Value = 0.1;
Parameters.gammam1.Value = gamma^(-1);
Parameters.km1.Value = k^(-1);
Parameters.betainit.Value = beta(1);
Parameters.SInitProp.Value = 0.7;
Parameters.EInitProp.Value = I0/TotPop;
Parameters.IInitProp.Value = I0/TotPop;
Parameters.RInitProp.Value = 1-(Parameters.SInitProp.Value+Parameters.EInitProp.Value+Parameters.IInitProp.Value);
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)

    

SEIRModel = struct();
SEIRModel.InitializeParameters = @SEIRInitialize;
% SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
SEIRModel.LikFunction = 'normpdf(log(Variables(:,5)),log(Data.Observations(5,IndTime)),Parameters.SigmaObs.Value)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
% Parameters.SigmaRW.Value = 0.1327;
% Parameters.SigmaRW.Min = -10^14;
% Parameters.SigmaRW.Max = 10^14;
% Parameters.SigmaRW.MinLim = 0;
% Parameters.SigmaRW.MaxLim = 50;
% Parameters.SigmaRW.Estimated = 1;
% Parameters.SigmaRW.TransfType = 'Log';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

%% Go

% 1: knowing everything

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
% Parameters.SigmaRW.Value = 0.1327;
% Parameters.SigmaRW.Min = -10^14;
% Parameters.SigmaRW.Max = 10^14;
% Parameters.SigmaRW.MinLim = 0;
% Parameters.SigmaRW.MaxLim = 50;
% Parameters.SigmaRW.Estimated = 1;
% Parameters.SigmaRW.TransfType = 'Log';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Parameters.PMCMC = 'PMMH';
Parameters.GMeth = 'cst given';
Parameters.G = 0.05^(-1);
Parameters.MCMCType = 'Rand';
Parameters.Epsil = 1;
Parameters.ModelType = 'SMC';
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
Parameters.AdaptC = 0.99;
ResPMMH = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,1000);


PlotMarc(ResPMMH,8)

ind = 30;
temp = AutoCorrelation(ResPMMH.Paths(:,5,ind),300);
ESS = length(ResPMMH.Coalescence)/(1+2*sum(temp))
    
    


% 2: not knowing sigma and noise

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.SigmaRW.Value = 0.1327;
Parameters.SigmaRW.Min = -10^14;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.MinLim = 0;
Parameters.SigmaRW.MaxLim = 50;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Log';
Parameters.SigmaObs.Value = Tau;
Parameters.SigmaObs.Min = -10^14;
Parameters.SigmaObs.Max = 10^14;
Parameters.SigmaObs.MinLim = 0;
Parameters.SigmaObs.MaxLim = 50;
Parameters.SigmaObs.Estimated = 1;
Parameters.SigmaObs.TransfType = 'Log';

Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Parameters.PMCMC = 'PMMH';
Parameters.GMeth = 'cst given';
Parameters.G = 0.05^(-1)*eye(2);
Parameters.MCMCType = 'Rand';
Parameters.Epsil = 1;
Parameters.ModelType = 'SMC';
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
Parameters.AdaptC = 0.99;
ResPMMH = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,100);

subplot(5,2,1)
i = 1;
toplot = 5;
ResPMMH.BetatPath = beta;
Res = ResPMMH;
Paths = Res.Paths;
PathsInstant = 0:Res.Parameters.ComputationTStep:sum(Data.NbComputingSteps);
ciplot(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.025),quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.975),[172,215,255]/255)
hold on
ciplot(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.25),quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.75),[100,153,251]/255)
plot(mean(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps))))),'k','LineWidth',2)
plot(Data.Observations(5,:),'g','LineWidth',2)
hold off
xlim([1 size(Data.Observations,2)])
try
    set(gca,'XTick',[delta:delta:length(Data.Dates)])
%         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
    set(gca,'XTickLabel',dates)
end
title('Estimated Total Influenza Incidence')
subplot(5,2,2)
ciplot(quantile(squeeze(exp(Paths(:,6,:))),0.025),quantile(squeeze(exp(Paths(:,6,:))),0.975),[172,215,255]/255)
hold on
ciplot(quantile(squeeze(exp(Paths(:,6,:))),0.25),quantile(squeeze(exp(Paths(:,6,:))),0.75),[100,153,251]/255)
plot(mean(squeeze(exp(Paths(:,6,:)))),'k','LineWidth',2)
t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
try
    plot(Res.BetatPath,'g','LineWidth',2) 
end
hold off
xlim([0 Data.Instants(end)])
title('\beta_t')
ResPMMH.BetatPath = beta;
subplot(5,2,3)
ESSs = [];
for i = 2:length(Data.Instants)-1
    tmp = autocorrelation(ResPMMH.Paths(:,6,Data.Instants(i)),10);
    ESSs(i) = 100/(1+2*sum(tmp(2:end)));
end
plot(ESSs(2:end));
xlabel('time')
ylabel('ESS')
subplot(5,2,4)
plot(ResPMMH.Thetas(2,:),ResPMMH.Thetas(1,:),'.')
xlabel('\sigma')
ylabel('noise')
subplot(5,2,5)
plot(ResPMMH.Thetas(2,:),ResPMMH.Paths(:,1,150),'.')
xlabel('\sigma')
ylabel('\beta(T/2)')
subplot(5,2,6)
plot(ResPMMH.Thetas(1,:),ResPMMH.Paths(:,1,150),'.')
xlabel('noise')
ylabel('\beta(T/2)')
subplot(5,2,7)
plot(ResPMMH.Thetas(2,:))
tmp = 100/(1+2*sum(autocorrelation(ResPMMH.Thetas(2,:),100)));
title(['\sigma ESS= ' num2str(tmp,3) '%'])
xlabel('iterations')
ylabel('\sigma')
subplot(5,2,8)
plot(ResPMMH.Thetas(1,:))
tmp = 100/(1+2*sum(autocorrelation(ResPMMH.Thetas(1,:),100)));
title(['noise ESS= ' num2str(tmp,3) '%'])
xlabel('iterations')
ylabel('noise')
subplot(5,2,9)
[F,XI]=KSDENSITY(ResPMMH.Thetas(2,:)) ;
plot(XI,F)
hold on
plot(0.1, 0,'og')
hold off
xlabel('\sigma')
subplot(5,2,10)
[F,XI]=KSDENSITY(ResPMMH.Thetas(1,:)) ;
plot(XI,F)
hold on
plot(0.05, 0,'og')
hold off
xlabel('noise')

plot(ResPMMH.Thetas(1,:))
title('\sigma traceplot')
subplot(3,2,4)
plot(ResPMMH.Thetas(2,:))
title('noise traceplot')
subplot(3,2,5)
hist(ResPMMH.Thetas(1,:))
title('\sigma posterior')
subplot(3,2,6)
hist(ResPMMH.Thetas(2,:))
title('noise posterior')

ind = 30;
temp = AutoCorrelation(ResPMMH.Paths(:,5,ind),300);
ESS = length(ResPMMH.Coalescence)/(1+2*sum(temp))
    
    

    
    
    


