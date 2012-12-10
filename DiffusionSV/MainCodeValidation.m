% validating the code for fractional BM


cd('/Users/dureaujoseph/AllScripts')
addpath([pwd '/DiffusionSV'])

%% sampling from BM - comparing with Wikipedia plots$

N = 1000;
step = 0.1;
Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative



H = 0.75;
sigma_X = 0.1;
Z = Sample_Z(N);
Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh,sigma_X);
Y = SampleObs(X,step,Vol);
subplot(2,1,1)
plot(X)
subplot(2,1,2)
plot(Y)


H = 0.95;
Z = Sample_Z(N);
Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh,sigma_X);
Y = SampleObs(X,step,Vol);
subplot(2,1,1)
plot(X)
subplot(2,1,2)
plot(Y)

H = 0.55;
Z = Sample_Z(N);
Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh,sigma_X);
Y = SampleObs(X,step,Vol);
subplot(2,1,1)
plot(X)
subplot(2,1,2)
plot(Y)

H = 0.25;
Z = Sample_Z(N);
Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh,sigma_X);
Y = SampleObs(X,step,Vol);
subplot(2,1,1)
plot(X)
subplot(2,1,2)
plot(Y)



%% validating score computation with numerical computation

Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative
nobs = 10;
step = 0.01;
N = nobs/step;
H = 0.55;
sigma_X = 0.1;

SimDatafBM('scoreTest.mat',N,step,Vol,H,sigma_X);


ratios = [];

for i = 1:50
    % sample data
    H = 0.55;
    sigma_X = 0.4;

    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,H);
    X = Bh_to_X(Bh,sigma_X);
    Y = SampleObs(X,step,Vol);


    % dL/dZ 
    ind = ceil(rand(1,1)*N); 
    Par.H.Value = H;
    Par.sigma_X.Value = sigma_X;
    LogLik1 = ComputeLogLikZ(Z,Y,Vol,Par);
    epsil = 0.000001;
    Z2 = Z;
    Z2(ind) = Z2(ind) + epsil;
    LogLik2 = ComputeLogLikZ(Z2,Y,Vol,Par);
    ScoreNumerical = (LogLik2 - LogLik1)/epsil;
    ScoreAnalytic = ComputeScore(Z,Y,Vol,VolDer,Par);
    disp(['ratio between the two gradient estimates (comp ' num2str(ind) '): ' num2str(ScoreAnalytic(ind)/ScoreNumerical,10)]);
    ratios(i)=ScoreAnalytic(ind)/ScoreNumerical;
end

[fi,xi] = ksdensity(ratios);
plot(xi,fi)
title(['Density of numerical/analytic score ratios (H=' num2str(H) ')'],'FontSize',20)


% sigma_X score
Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative
nobs = 100;
step = 0.01;
N = nobs/step;
H = 0.55;
sigma_X = 0.1;

Par = struct();
Par.Names.All = {'H','sigma_X'};
Par.H.MinLim = 0;
Par.H.MaxLim = 1;
Par.H.Transf = @logit;
Par.H.InvTransf = @invlogit;
Par.H.Corr = @logitCorr;

Par.sigma_X.MinLim = 0;
Par.sigma_X.MaxLim = 1;
Par.sigma_X.Transf = @mylog;
Par.sigma_X.InvTransf = @invlog;
Par.sigma_X.Corr = @logCorr;

Par.GradCorr = 0;

Par.theta_sampler='JointHMC';

ratios = [];

for i = 1:100
    H = 0.55;
    sigma_X = 0.05+0.2*rand;

    
    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,H);
    X = Bh_to_X(Bh,sigma_X);
    Y = SampleObs(X,step,Vol);
    
    Par.H.Value = H;
    Par.H.Estimated = 0;
    Par.sigma_X.Value = sigma_X;
    Par.sigma_X.Estimated = 1;
    Par = DefineIndexes(Par);
    Par = NoTransfToTransf(Par);
    ScoreAnalytic = ComputeScore(Z,Y,Vol,VolDer,Par);

    LogLik1 = ComputeLogLikZ(Z,Y,Vol,Par);
    epsil = 0.0000001;    
    Par2 = Par;
    Par2.sigma_X.Value = Par.sigma_X.Value + epsil;
    Par2 = NoTransfToTransf(Par2);
    LogLik2 = ComputeLogLikZ(Z,Y,Vol,Par2);
    ScoreNumerical = (LogLik2 - LogLik1)/epsil;
    disp(['ratio between the two gradient estimates: ' num2str(ScoreAnalytic(end)/ScoreNumerical,10)]);
    ratios(i)=ScoreAnalytic(end)/ScoreNumerical;
    
end


clf
[fi,xi] = ksdensity(ratios);
plot(xi,fi)
title(['Density of numerical/analytic \sigma_X - gradient ratios'],'FontSize',20)



% sigma_X score with correction
Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative
nobs = 100;
step = 0.01;
N = nobs/step;
H = 0.55;
sigma_X = 0.1;

Par = struct();
Par.Names.All = {'H','sigma_X'};
Par.H.MinLim = 0;
Par.H.MaxLim = 1;
Par.H.Transf = @logit;
Par.H.InvTransf = @invlogit;
Par.H.Corr = @logitCorr;

Par.sigma_X.MinLim = 0;
Par.sigma_X.MaxLim = 1;
Par.sigma_X.Transf = @mylog;
Par.sigma_X.InvTransf = @invlog;
Par.sigma_X.Corr = @logCorr;

theta_sampler='JointHMC' ;
Par.theta_sampler = theta_sampler;

ratios = [];
Par.GradCorr = 1;

for i = 1:50
    H = 0.65;
    sigma_X = 0.1;

    
    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,H);
    X = Bh_to_X(Bh,sigma_X);
    Y = SampleObs(X,step,Vol);
    
    Par.JointHMC = 1;
    Par.H.Value = H;
    Par.H.Estimated = 0;
    Par.sigma_X.Value = sigma_X;
    Par.sigma_X.Estimated = 1;
    Par = DefineIndexes(Par);
    Par = NoTransfToTransf(Par);
    
    ScoreAnalytic = ComputeScore(Z,Y,Vol,VolDer,Par);

    LogLik1 = ComputeLogLikZ(Z,Y,Vol,Par);
    epsil = 0.0000001;    
    Par2 = Par;
    Par2.sigma_X.TransfValue = Par.sigma_X.TransfValue + epsil;
    Par2 = TransfToNoTransf(Par2);
    LogLik2 = ComputeLogLikZ(Z,Y,Vol,Par2);
    ScoreNumerical = (LogLik2 - LogLik1)/epsil;
    disp(['ratio between the two gradient estimates: ' num2str(ScoreAnalytic(end)/ScoreNumerical,10)]);
    ratios(i) = ScoreAnalytic(end)/ScoreNumerical;    
end


[fi,xi] = ksdensity(ratios);
plot(xi,fi)
title(['Density of numerical/analytic log(sigma_X) - gradient ratios'],'FontSize',20)



% H score
Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative
nobs = 100;
step = 0.01;
N = nobs/step;
H = 0.55;
sigma_X = 0.1;

Par = struct();
Par.Names.All = {'H','sigma_X'};
Par.H.MinLim = 0;
Par.H.MaxLim = 1;
Par.H.Transf = @logit;
Par.H.InvTransf = @invlogit;
Par.H.Corr = @logitCorr;

Par.sigma_X.MinLim = 0;
Par.sigma_X.MaxLim = 1;
Par.sigma_X.Transf = @mylog;
Par.sigma_X.InvTransf = @invlog;
Par.sigma_X.Corr = @logCorr;

Par.theta_sampler='JointHMC' ;
Par.GradCorr = 0;

ratios = [];

for i = 1:300
    H = 0.65+0.1*randn;
    sigma_X = 0.1;

    
    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,H);
    X = Bh_to_X(Bh,sigma_X);
    Y = SampleObs(X,step,Vol);
    
    Par.JointHMC = 1;
    Par.H.Value = H;
    Par.H.Estimated = 1;
    Par.H.Index = 1;
    Par.sigma_X.Value = sigma_X;
    Par.sigma_X.Estimated = 0;
    Par.sigma_X.Index = 1;
    Par = DefineIndexes(Par);
    Par = NoTransfToTransf(Par);
    ScoreAnalytic = ComputeScore(Z,Y,Vol,VolDer,Par);

    LogLik1 = ComputeLogLikZ(Z,Y,Vol,Par);
    epsil = 0.0000001;    
    Par2 = Par;
    Par2.H.Value = Par.H.Value + epsil;
    LogLik2 = ComputeLogLikZ(Z,Y,Vol,Par2);
    ScoreNumerical = (LogLik2 - LogLik1)/epsil;
    disp(['ratio between the two gradient estimates: ' num2str(ScoreAnalytic(end)/ScoreNumerical,10)]);
    ratios(i) = ScoreAnalytic(end)/ScoreNumerical;
    
end


[fi,xi] = ksdensity(ratios);
plot(xi,fi)
title(['Density of numerical/analytic H - gradient ratios'],'FontSize',20)


% H score with correction
Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative
nobs = 100;
step = 0.01;
N = nobs/step;
H = 0.55;
sigma_X = 0.1;

Par = struct();
Par.Names.All = {'H','sigma_X'};
Par.H.MinLim = 0;
Par.H.MaxLim = 1;
Par.H.Transf = @logit;
Par.H.InvTransf = @invlogit;
Par.H.Corr = @logitCorr;

Par.sigma_X.MinLim = 0;
Par.sigma_X.MaxLim = 1;
Par.sigma_X.Transf = @mylog;
Par.sigma_X.InvTransf = @invlog;
Par.sigma_X.Corr = @logCorr;

theta_sampler='JointHMC' ;
Par.theta_sampler = theta_sampler;

ratios = [];
Par.GradCorr = 1;

for i = 1:50
    H = 0.65;
    sigma_X = 0.1;

    
    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,H);
    X = Bh_to_X(Bh,sigma_X);
    Y = SampleObs(X,step,Vol);
    
    Par.JointHMC = 1;
    Par.H.Value = H;
    Par.H.Estimated = 1;
    Par.sigma_X.Value = sigma_X;
    Par.sigma_X.Estimated = 0;
    Par = DefineIndexes(Par);
    Par = NoTransfToTransf(Par);
    
    ScoreAnalytic = ComputeScore(Z,Y,Vol,VolDer,Par);

    LogLik1 = ComputeLogLikZ(Z,Y,Vol,Par);
    epsil = 0.0000001;    
    Par2 = Par;
    Par2.H.TransfValue = Par.H.TransfValue + epsil;
    Par2 = TransfToNoTransf(Par2);
    LogLik2 = ComputeLogLikZ(Z,Y,Vol,Par2);
    ScoreNumerical = (LogLik2 - LogLik1)/epsil;
    disp(['ratio between the two gradient estimates: ' num2str(ScoreAnalytic(end)/ScoreNumerical,10)]);
    ratios(i) = ScoreAnalytic(end)/ScoreNumerical;    
end


[fi,xi] = ksdensity(ratios);
plot(xi,fi)
title(['Density of numerical/analytic logit(H) - gradient ratios'],'FontSize',20)




%% confirming that cost is O(nlogn)



Par.H.Value = Data.Htrue;
Par.sigma_X.Value = Data.sigma_Xtrue;
Par = NoTransfToTransf(Par);

Par.loop = 20000;
Par.hH = 0.00000000004;
Par.hZ = 0.1;
Par.hsig = 1;
Par.nsteps = 1;




Ns = (1:10)*10000;
step = 0.01;

ts_sample = [];
ts_deriv = [];

for i = 1:length(Ns)
    i
    N = Ns(i);
    
    tic;
    Z = Sample_Z(N);
    Bh = Z_to_Bh(Z,N,step,H);
    X = Bh_to_X(Bh,sigma_X);
    Y = SampleObs(X,step,Vol);
    ts_sample(i)=toc;
    
    tic;
    ScoreAnal = ComputeScore(Z,Y,Vol,VolDer,H,sigma_X);
    ts_deriv(i)=toc;
end

subplot(2,1,1)
plot(Ns.*log(Ns),ts_sample,'-og')
xlabel('Nlog(N)','FontSize',20)
ylabel('Time for Sample','FontSize',20)
subplot(2,1,2)
plot(Ns.*log(Ns),ts_deriv,'-or')
xlabel('Nlog(N)','FontSize',20)
ylabel('Time for Score','FontSize',20)

%% playing with HMC


%Create data and likelihood components objects etc
nobs = 100; % Y(0) = 0 is counted as an observation
step = 0.05;
N = (nobs-1)/step;
Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative


%%% MODES SELECTION

% Qsampler='MALA';
Qsampler='HybridMC';

theta_sampler='JointHMC'; % JointHMC or GibbsRW
theta_sampler='GibbsRW';
theta_sampler='GibbsHMC';
theta_sampler='Blocks';

% PARAMETERS
Htrue = 0.6;
sigma_X = 0.1;
obstep = 1;
loop=20000;


data_file=strcat('fBMDataNoFine100_H=',num2str(Htrue),'_sigma=',num2str(sigma_X),'.mat');

data_files = {'fBMData_H=0.6_sigma=0.1.mat'}

% DATA
SimDatafBM(data_file,N,step,Vol,Htrue,sigma_X);

SavePath = '/Users/dureaujoseph/Documents/PhD_Data/fBM/';


% % datafiles = {'fBMData_H=0.6_sigma=0.1.mat','fBMData_H=0.6_sigma=0.01.mat','fBMData_H=0.8_sigma=0.1.mat'};
% % method = {'MALA','HMC'};
% % 
% % indData = 1;
% load([SavePath '/' datafiles{indData}]);

load([SavePath '/' data_file]);



Par = struct();
Par.Names.All = {'H','sigma_X'};
Par.H.MinLim = 0.4;
Par.H.MaxLim = 1;
Par.H.Transf = @logit;
Par.H.InvTransf = @invlogit;
Par.H.Corr = @logitCorr;
Par.H.Estimated = 1;
Par.sigma_X.MinLim = 0;
Par.sigma_X.MaxLim = 100;
Par.sigma_X.Transf = @mylog;
Par.sigma_X.InvTransf = @invlog;
Par.sigma_X.Corr = @logCorr;
Par.sigma_X.Estimated = 1;
Par = DefineIndexes(Par);

Par.H.Value = Data.Htrue;
Par.sigma_X.Value = Data.sigma_Xtrue;
Par = NoTransfToTransf(Par);

Par.loop = loop;
Par.Qsampler = Qsampler;
if strcmp(Qsampler,'MALA')
    Par.nsteps = 1;
else
    Par.nsteps = 20;
end

Par.theta_sampler = theta_sampler;
if strcmp(theta_sampler,'JointHMC')
    Par.h=0.06;
elseif strcmp(theta_sampler,'GibbsHMC')
    Par.hZ=0.19;
    Par.htheta=0.13;
elseif strcmp(theta_sampler,'GibbsRW')
    Par.hH = 0.5;
    Par.hZ = 0.1;
    Par.hsig = 1;
elseif strcmp(theta_sampler,'Blocks')
    Par.d = 10;
end


Res = RunJointMCMC(Data,Par);
PlotfBMoutput(Res)


save([SavePath '/Res_JointHMC20_Hsigest_' data_file],'Res')


est = 'No';
load([SavePath '/Res_RW_' est 'est_' data_file])
ResRW = Res;
load([SavePath '/Res_MALA_' est 'est_' data_file])
ResMALA = Res;
load([SavePath '/Res_HMC5_' est 'est_' data_file])
ResHMC5 = Res;
load([SavePath '/Res_HMC10_' est 'est_' data_file])
ResHMC10 = Res;
load([SavePath '/Res_HMC20_' est 'est_' data_file])
ResHMC20 = Res;

clf
[fi,xi] = ksdensity(ResRW.out_Hs);
plot(xi,fi,'k')
hold on
[fi,xi] = ksdensity(ResMALA.out_Hs);
plot(xi,fi,'g')
[fi,xi] = ksdensity(ResHMC.out_Hs);
plot(xi,fi,'b')
hold off
ylabel('p(H|y)','FontSize',20)
xlabel('H','FontSize',20)
h = legend('Gibbs-RW: ESS=0.4%','Joint-MALA: ESS=0.5%','Joint-HMC: ESS=6.4%');
set(h,'FontSize',20)

load([SavePath 'Res_MALA_Hsigest_fBMData1000_H=0.6_sigma=0.1.mat'])

Ress = {ResRW,ResMALA,ResHMC5,ResHMC10,ResHMC20};
Ress = {Res,Res};
ESSs = {};

for k = 2:2
    k
%     Ress{k}.h
    %figure(11); for i=10:10:min(500,loop); plot(out_Q(i,:),'b');hold on; end; plot(out_Q(1,:),'r');hold off
    ESS=zeros(size(Ress{k}.out_Zs,2),1);
    for i=1:size(Ress{k}.out_Zs,2)
        r=sum(autocorr(Ress{k}.out_Zs(:,i),1000));
        ESS(i)=100/(1+2*r);
    end
    ESSs{k} = ESS;
%     figure(12);plot(ESS);title('ESS (%) for each v over time')   
%     disp([min(ESS),median(ESS),max(ESS)])   
    r=sum(autocorr(Ress{k}.out_Hs,1200));
    disp(100/(1+2*r));
    r=sum(autocorr(Ress{k}.out_sigs,1200));
    disp(100/(1+2*r));
end
    
plot(ESSs{1},'k')
hold on
plot(ESSs{2},'g')
plot(ESSs{3},'b')
hold off
ylabel('ESS(Z_t)','FontSize',20)
xlabel('t','FontSize',20)
h = legend('Gibbs-RW: ESS_{min}=0.1%','Joint-MALA: ESS_{min}=0.1%','Joint-HMC: ESS_{min}=1.1%');
set(h,'FontSize',20)




load([SavePath 'Res_MALA_Hsigest_fBMData1000_H=0.6_sigma=0.1.mat'])
Res1000 = Res;
load([SavePath 'Res_JointHMC_Hsigest_fBMData_H=0.6_sigma=0.1.mat'])

clf
subplot(3,1,3)
[fi,xi]=ksdensity(Res.out_Hs);
plot(xi,fi);
hold on
[fi,xi]=ksdensity(Res1000.out_Hs);
plot(xi,fi,'g');
hold off
title('H','FontSize',20)
subplot(3,1,3)
[fi,xi]=ksdensity(Res.out_sigs);
plot(xi,fi);
hold on
[fi,xi]=ksdensity(Res1000.out_sigs);
plot(xi,fi,'g');
hold off
title('\sigma_X','FontSize',20)
h = legend('250 observations','1000 observations');
set(h,'FontSize',20)
subplot(3,1,3)
plot(Res.out_Hs,Res.out_sigs,'.')
hold on
plot(Res1000.out_Hs,Res1000.out_sigs,'.g')
hold off
xlabel('H')
ylabel('\sigma_X')

% r=sum(autocorr(out_kappa,100));
% ESSk=100/(1+2*r);
% %disp(ESSk)
% disp([min(ESS),median(ESS),max(ESS),ESSk])
end 

