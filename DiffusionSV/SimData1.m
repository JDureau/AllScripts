function [X,V]=SimData1(data_file,kappa,sigma,T,npoints,M)

% npoints: no. of total points except from the initial one.
% m: no. of imputed points between each observed point. 
% T: ending time
% kappa sigma: parameters
obstep=T/npoints;%must be integer
step=obstep/(M+1);
sstep=sqrt(step);
V=zeros(1,npoints*(M+1)+1);
v=V(1);
Rands=randn(1,npoints*(M+1));
for i=1:npoints*(M+1)
    v=v-step*kappa*v+sigma*sstep*Rands(i);
    V(i+1)=v;

    if (rem(i,10000) == 0)
        disp((100*i)/(T*(M+1)))
    end
end

X=zeros(1,npoints+1);
Rands=randn(1,npoints+1);
x=X(1);
for i=1:npoints
    vi=V((i-1)*(M+1)+1:i*(M+1));
    sd=sqrt(sum(exp(vi))*step);
    x=x+sd*Rands(i);
    X(i+1)=x;
end

save(data_file,'X','V','kappa','sigma','T','npoints');













    
    
    