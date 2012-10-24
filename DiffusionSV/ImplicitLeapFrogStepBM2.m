function [q,qz,u,Vderiv]=ImplicitLeapFrogStepBM2(q0,qz0,u0,h,VDeriv0,a1,b1,c1,X,kappa,sigma,step,npoints,m,T)
% for the non bridge case
Z=[0,TDMAsolver(a1,b1,c1,VDeriv0(2:end))];
qz=(u0-(0.5*h*Z)+((1/h)-0.25*h)*qz0)/((1/h)+0.25*h);
q=qz+q0(1);
Vderiv=VarDer(X,q,kappa,sigma,step,npoints,m);
Z=[0,TDMAsolver(a1,b1,c1,Vderiv(2:end))];
u=(qz-qz0)/h;
u=u-(0.5*h*Z)-0.25*h*(qz0+qz);