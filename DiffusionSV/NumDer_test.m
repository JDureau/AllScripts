clear all
close all

kappa=0.2;
sigma=1;
T=1; %end time
npoints=200; %number of observed points apart from the initial
m=99; %number of imputed points for simulating
h=1e-6;
data_file=strcat('SVData_T=',num2str(T));

load(data_file)

obstep=T/npoints;%must be integer
step=obstep/(m+1);
Q=V(1:10:end)/sigma;

AnDeriv=VarDer(X,Q,kappa,sigma,step,npoints,m);
NumDeriv=AnDeriv;
for i=2:npoints*(m+1);
    Qph=Q;
    Qph(i)=Q(i)+h;
    Qmh=Q;
    Qmh(i)=Q(i)-h;
    NumDeriv(i)=(logLikeEval(X,Qph,kappa,sigma,step,npoints,m)-logLikeEval(X,Qmh,kappa,sigma,step,npoints,m))/(2*h);
end
figure(1);plot([(AnDeriv)',NumDeriv'])
disp([sum(AnDeriv),sum(NumDeriv)])
