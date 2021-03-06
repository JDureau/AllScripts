%overlapping blocks
%100 points, volatility endpoints unknown.




cd('/Users/dureaujoseph/AllScripts')
addpath([pwd '/DiffusionSV/'])



clear all
close all

draw_data='off';
plot_data='off';
run_mcmc='on';
load_output='on';
%Qsampler='test';
%Qsampler='ind_sampler';
%Qsampler='rw_diff';
Qsampler='MALA';
%Qsampler='HybridMC';
theta_sampler='none';
%theta_sampler='rw';
%theta_sampler='Gibbs';
%theta_sampler='';

loop=2000;
m=9;

%Set volartility parameters etc
kappa=0.03;  
sigma=sqrt(0.03);
T=250; %end time
npoints=250; %number of observed points apart from the initial
block=npoints; %block size in terms of observed points
blocks=1:round(block/2):npoints+1-block;
nblocks=length(blocks);
M=999; %number of imputed points for simulating
data_file=strcat('SVData_T=',num2str(T));

%Simulate data or load them
if strcmp(draw_data,'on')
    [X,V]=SimData1(data_file,kappa,sigma,T,npoints,M);
else
    load(data_file)
end

if strcmp(plot_data,'on')
    obstep=T/npoints;%must be integer
    step=obstep/(M+1);
    subplot(2,1,1);plot(0:obstep:T,X,'*')
    subplot(2,1,2);plot(0:step:T,V)
end    

if strcmp(run_mcmc,'on')
    %Set MCMC parameters    
    r=0.175;
    %st d' for the parameter random walk proposals
    pstd_kappa=0.3; %for ind_sampler
    pstd_sigma=0.7;
    %MALA/HybridMC parameters
    %Thomas algorithm objects (for inverting BM cov etc)
    obstep=T/npoints;%must be integer
    step=obstep/(m+1);
    time_scale_adj=(m+1)/obstep;
    %for the Brownian bridge case
    a1=[0;-ones(block*(m+1)-2,1)]*time_scale_adj;
    b1=2*ones(block*(m+1)-1,1)*time_scale_adj;
    c1=[-ones(block*(m+1)-2,1);0]*time_scale_adj;
    %for the Brownian motion case
    a2=[0;-ones(block*(m+1)-1,1)]*time_scale_adj;
    b2=[2*ones(block*(m+1)-1,1);1]*time_scale_adj;
    c2=[-ones(block*(m+1)-1,1);0]*time_scale_adj;
    
    %MALA HMC tuning parameters
    if strcmp(Qsampler,'MALA')
        %h=2.0*ones(nblocks,1);h(nblocks)=0.75;h(1)=3.0;h(76)=1.3;h(77)=1.5;
        h=0.085*ones(nblocks,1);
        nsteps=1;
    elseif strcmp(Qsampler,'HybridMC')
        h=0.075*ones(nblocks,1);
        nsteps=20;
    end
    %initial values
    V=V(1:100:end); %for the path
    %kappa=1;
    %sigma=0.3;
    output_file=strcat(data_file,'MCMCoutput_m=',num2str(m),'_loop=',num2str(loop),Qsampler,theta_sampler);
    
    %Run MCMC
    %RunMCMC
    tic;RunMCMC;toc
    save(output_file,'out*')    
end

%MCMC output
if strcmp(load_output,'on')
    output_file=strcat(data_file,'MCMCoutput_m=',num2str(m),'_loop=',num2str(loop),Qsampler,theta_sampler);
    load(output_file)
    
    %Analyze output
    enditer=loop;
    %enditer=2000;
%    disp(acrate_ov(out_kappa(1:enditer)))%,out_sigma(1:enditer)]))%,
    disp(acrate_ov(out_Q(1:min(500,enditer),5)))  %acceptance rates  
    ACRPath=acrate_ov(out_Q(:,m+2:m+1:npoints*(m+1)+1));
    figure(1);plot(ACRPath);title('paths acceptance probabilities')   
%    figure(2);plot(out_loglike(1:enditer));title('postrerior draws of log likelihood')
%    figure(3);plot(out_kappa(1:enditer));title('postrerior draws of kappa')
%     figure(4);plot(out_sigma(1:enditer));title('postrerior draws of sigma')        
%    figure(5);autocorr(out_kappa(1:enditer),100);title('autocorrelation of postrerior draws of kappa')
%     figure(6);autocorr(out_sigma(1:enditer),100);title('autocorrelation of postrerior draws of sigma')
%    figure(7);autocorr(out_loglike(1:enditer),100);title('autocorrelation of postrerior draws of log-likelihood')
%     figure(8);[f,x]=ksdensity(out_kappa(201:enditer));plot(x,f);title('postrerior density of kappa')
% %    figure(8);[f,x]=ksdensity(out_kappa(201:enditer));plot(x,f,x,exppdf(x));title('postrerior density of kappa versus prior')
% %     figure(9);[f,x]=ksdensity(out_sigma(201:enditer));plot(x,f,x,exppdf(x));title('postrerior density of sigma versus prior')        
figure(10);subplot(3,1,1);plot(out_Q(:,m+2));subplot(3,1,2);plot(out_Q(:,round(npoints*m/2)));
subplot(3,1,3);plot(out_Q(:,end))
%figure(11); for i=10:10:min(500,loop); plot(out_Q(i,:),'b');hold on; end; plot(out_Q(1,:),'r');hold off
ESS=zeros(npoints(1)-1,1);
for i=1:npoints(1)-1
    r=sum(autocorr(out_Qm(:,i),1200));
    ESS(i)=100/(1+2*r);
end
figure(12);plot(ESS);title('ESS (%) for each v over time')   
disp([min(ESS),median(ESS),max(ESS)])
% r=sum(autocorr(out_kappa,100));
% ESSk=100/(1+2*r);
% %disp(ESSk)
% disp([min(ESS),median(ESS),max(ESS),ESSk])
end 

