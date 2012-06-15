function [S,E,I,R,X]=ODESEIR(beta,gamma,k,I0,Tot,N,m,step)
S=zeros(1,N*(m+1)+1);
E=zeros(1,N*(m+1)+1);
I=zeros(1,N*(m+1)+1);
R=zeros(1,N*(m+1)+1);
C=zeros(1,N*(m+1)+1);
X=zeros(1,N);
S(1)=0.7*Tot;E(1)=I0;I(1)=I0;R(1)=0.3*Tot-2*I0;

for i = 2:N*(m+1)+1
    b = beta(i-1);
    S(i) = S(i-1) + (-b*S(i-1)*I(i-1)/Tot)*step ;
    E(i) = E(i-1) + ( (b*S(i-1)*I(i-1)/Tot)-k*E(i-1))*step ;
    I(i) = I(i-1) + (-gamma*I(i-1) + k*E(i-1))*step ;   
    if min([S(i),E(i),I(i)])<0
        disp('negative values in ODE')
    end
    R(i) = R(i-1) + ( gamma*I(i-1))*step ;
%    Observations(IndIt) = Variables(5)*(1+sigmaobs*randn(1,1));
end 
for i=1:N
    X(i)=sum(k*E((i-1)*(m+1)+1:i*(m+1)))*step;
end