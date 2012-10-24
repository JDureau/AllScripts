clear all
close all

kappa=0.2;
sigma=1;
T=10; %end time
npoints=10; %number of observed points apart from the initial
m=999; %number of imputed points for simulating
data_file=strcat('SVData_T=',num2str(T));

load(data_file)


obstep=T/npoints;
step=obstep/(m+1);
U=V(1:1:end)/sigma;
kappas=0.02:0.01:1.6;
klogL=zeros(1,length(kappas));
iter=0;
for kappa=kappas
    iter=iter+1;
    MU1 = -kappa*U;
    dMU1 = -kappa;    
    lG1 = logGirs(MU1,dMU1,U(end),kappa,step);
    klogL(iter)=lG1;
end
figure(1);plot(kappas,klogL);title('profile log-likelihood for kappa')

% sigmas=0.1:0.01:0.6;
% slogL=zeros(1,length(sigmas));
% iter=0;
% for sigma=sigmas
%     iter=iter+1;
%     lF1=0;
%     for i=1:npoints
%         lF1=lF1+logF(X(i+1)-X(i),U((i-1)*(m+1)+1:i*(m+1)),sigma,step);
%     end
%     slogL(iter)=lF1;
% end
% figure(2);plot(sigmas,slogL);title('profile log-likelihood for sigma')