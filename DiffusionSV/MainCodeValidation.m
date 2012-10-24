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

N = 100;
step = 0.01;

% sample data
H = 0.55;
Z = Sample_Z(N);
Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh);
Y = SampleObs(X,step,sigma);




LogLik1 = ComputeLogLikZ(Z,Y,H,sigma);
epsil = 0.0000001;
Z2 = Z;
Z2(1) = Z2(1) + epsil;
LogLik2 = ComputeLogLikZ(Z2,Y,H,sigma);
ScoreNum = (LogLik2 - LogLik1)/epsil
ScoreAnal = ComputeScore(Z,Y,sigma,dersigma,H)



subplot(2,1,1)
plot(B)



X = Bh_to_X(B);

subplot(2,1,2)
plot(X)

