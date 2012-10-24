function x=logF(dy,Q,sigma,step)
vd=sum(exp(sigma*Q))*step;
x=-0.5*(log(vd)+(dy*dy/vd));