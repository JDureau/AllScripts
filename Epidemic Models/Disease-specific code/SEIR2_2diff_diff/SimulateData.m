tau=0.05; %noise st. deviation
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
    Y(i)=exp(log(X(i))+tau*Rands(i));%normal in the log scale
end