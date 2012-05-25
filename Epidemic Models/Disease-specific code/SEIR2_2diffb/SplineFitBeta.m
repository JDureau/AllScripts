function Spline = SplineFitBeta(Betas,N)

n = length(Betas);
delta = floor(n/N);
tsinit = delta:delta:delta*N;
xisinit = Betas(tsinit);
xbetas = 1:n;

xInit = [tsinit,xisinit];
[x,fval,exitflag,output] = fminsearch(@(x) SplineFitR2(x,Betas),xInit,optimset('MaxFunEvals',1000,'MaxIter',1000,'TolX',1e-7,'TolFun',1e-7));


Spline.tis = ceil(x(1:length(x)/2));
Spline.yis = x(length(x)/2+1:end);