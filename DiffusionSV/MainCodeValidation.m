% validating the code for fractional BM


cd('/Users/dureaujoseph/AllScripts')
addpath([pwd '/DiffusionSV'])

%% sampling from BM - comparing with Wikipedia plots$

N = 1000;
step = 0.1;

H = 0.75;
Z = Sample_Z(N);
Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh);
plot(X)


H = 0.95;
Z = Sample_Z(N);
Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh);
plot(X)

H = 0.55;
Z = Sample_Z(N);
Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh);
plot(X)

H = 0.25;
Z = Sample_Z(N);
Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh);
plot(X)



%% validating score computation with numerical computation

sigma = @ClassicVol; % how the volatility X plays on the price
dersigma = @DerClassicVol; % its derivative

nobs = 10;
step = 0.1;
N = nobs/step;


% sample data
H = 0.55;
Z = Sample_Z(N);
Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh);
Y = SampleObs(X,step,sigma);


% dL/dX
LogLik1 = ComputeLogLikX(X,Y,H,sigma);
epsil = 0.0000000001;
X2 = X;
X2(55) = X2(55) + epsil;
LogLik2 = ComputeLogLikX(X2,Y,H,sigma);
ScoreNum = (LogLik2 - LogLik1)/epsil
ScoreAnal = ComputeScore(Z,Y,sigma,dersigma,H)



% dL/dB
ind = 34; 
LogLik1 = ComputeLogLikB(Bh,Y,H,sigma);
epsil = 0.000000001;
Bh2 = Bh;
Bh2(ind) = Bh2(ind) + epsil;
LogLik2 = ComputeLogLikB(Bh2,Y,H,sigma);
ScoreNum = (LogLik2 - LogLik1)/epsil
ScoreAnal = ComputeScore(Z,Y,sigma,dersigma,H,ind)


% dL/dZ 
ind = 90; 
LogLik1 = ComputeLogLikZ(Z,Y,H,sigma);
epsil = 0.00000001;
Z2 = Z;
Z2(ind) = Z2(ind) + epsil;
LogLik2 = ComputeLogLikZ(Z2,Y,H,sigma);
ScoreNum = (LogLik2 - LogLik1)/epsil
ScoreAnal = ComputeScore(Z,Y,sigma,dersigma,H,ind)



subplot(2,1,1)
plot(B)



X = Bh_to_X(B);

subplot(2,1,2)
plot(X)

