function [q,qz,u,VDeriv]=ImplicitHybridMC(q0,qz0,u0,h,VDeriv0,a1,b1,c1,...
    X,kappa,sigma,step,npoints,m,T,nsteps)
q=q0;
qz=qz0;
VDeriv=VDeriv0;
u=u0;
for i=1:nsteps
    [q,qz,u,VDeriv]=ImplicitLeapFrogStepBM(q,qz,u,h,VDeriv,a1,b1,c1,X,kappa,sigma,step,npoints,m,T);
end

