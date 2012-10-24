function x=logBridgedensity(BB,a1,b1,c1,step)
B=BB(2:end-1);
Z=TDMAmultiplier(a1,b1,c1,B');
x=-0.5*B*Z;