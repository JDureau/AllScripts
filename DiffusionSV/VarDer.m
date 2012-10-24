function x=VarDer(X,Q,kappa,sigma,step,npoints,m)
%derivative of minus log likelihood
%The endpoint should be ignored in the bridge case. Initial point is 0. 
x=zeros(1,npoints*(m+1)+1);
for i=1:npoints
    q=Q((i-1)*(m+1)+1:i*(m+1)+1);
    Ii=sum(exp(q*sigma))*step;
    k1=sigma*exp(sigma*q)/(Ii);
    k2=(sigma*exp(sigma*q)*(X(i+1)-X(i))^2)/((Ii^2));
    x((i-1)*(m+1)+1:i*(m+1)+1)=-step*(0.5*(-k1+k2)-q*kappa^2);
end
x(npoints*(m+1)+1)=x(npoints*(m+1)+1)+kappa*q(end); %for the non-bridge case
x(1)=0;
%x(npoints*(m+1)+1)=0;    