function [q,qz,u,Vderiv]=ImplicitLeapFrogStepBM(q0,qz0,u0,h,VDeriv0,a1,b1,c1,X,kappa,sigma,step,npoints,m,T)
% for the bridge case
Z=[0,TDMAsolver(a1,b1,c1,VDeriv0(2:end-1)),0];
qz=(u0-(0.5*h*Z)+((1/h)-0.25*h)*qz0)/((1/h)+0.25*h);
q=eta(qz,T,npoints*(m+1)-1,q0(1),q0(end));
Vderiv=VarDer(X,q,kappa,sigma,step,npoints,m);
Z=[0,TDMAsolver(a1,b1,c1,Vderiv(2:end-1)),0];
u=(qz-qz0)/h;
u=u-(0.5*h*Z)-0.25*h*(qz0+qz);