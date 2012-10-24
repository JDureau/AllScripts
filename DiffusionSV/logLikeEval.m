function x=logLikeEval(X,U,kappa,sigma,step,npoints,m)
MU =-kappa*U;
dMU=-kappa;
lG = logGirs(MU,dMU,U(end),kappa,step);
lF=0;
for j=1:npoints
    lF=lF+logF(X(j+1)-X(j),U((j-1)*(m+1)+1:j*(m+1)),sigma,step);
end
x=lG+lF;
%x=lF;
%x=lG;

