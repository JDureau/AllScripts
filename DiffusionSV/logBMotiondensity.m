function x=logBMotiondensity(BB,step)
x=-0.5*(diff(BB)*diff(BB)')/step;