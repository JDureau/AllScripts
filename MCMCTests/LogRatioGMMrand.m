function LogRatio = LogRatioGMMrand(x,xstar,Parameters)

Sigmas = Parameters.Sigmas;
Nc = Parameters.Nc;
Dens = Parameters.Dens;

% numerator (xstar)
temp = 0;
posts = posterior(Dens,[x;xstar]);
if strcmp(Parameters.Mode,'Log')
    post = log(1+posts(1,:))/sum(log(1+posts(1,:)));
elseif strcmp(Parameters.Mode,'L')
    post = posts(1,:)/sum(posts(1,:));
end
    
for IndComp = 1:Nc
    temp = temp + (mvnpdf(xstar,x,Parameters.Epsil^2*(squeeze(Sigmas(:,:,IndComp))))*post(IndComp));
end
LogqTempStar = log(temp);
LogLikxstar = log(Parameters.f(xstar,Parameters));


% Star to Temp 
temp = 0;
if strcmp(Parameters.Mode,'Log')
    post = log(1+posts(2,:))/sum(log(1+posts(2,:)));
elseif strcmp(Parameters.Mode,'L')
    post = posts(2,:)/sum(posts(2,:));
end

for IndComp = 1:Nc
    temp = temp + (mvnpdf(x,xstar,Parameters.Epsil^2*(squeeze(Sigmas(:,:,IndComp))))*post(IndComp));
end
LogqStarTemp = log(temp);
LogLikx = log(Parameters.f(x,Parameters));

LogNum = LogLikxstar + LogqStarTemp;
LogDenom = LogLikx + LogqTempStar;
LogRatio = LogNum - LogDenom;

if isnan(LogRatio)
    disp('Ratio GMM rand problem')
end