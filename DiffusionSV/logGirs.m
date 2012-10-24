function x=logGirs(MU,dMU,Q_T,kappa,step)
k1=-0.5*kappa*Q_T*Q_T;
K=MU.^2+dMU;
x=k1-0.5*sum(K)*step;