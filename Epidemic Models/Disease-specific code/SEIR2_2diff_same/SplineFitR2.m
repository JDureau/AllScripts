function res = SplineFitR2(x,Betas)


n = length(Betas);
ts = ceil(x(1:length(x)/2));
[ts,I,J] = UNIQUE(ts);

yis = x(length(x)/2+1:end);
yis = yis(I);
xbetas = 1:n;

BetaFit = spline(ts,yis,xbetas);

res = sum((Betas-BetaFit).^2);