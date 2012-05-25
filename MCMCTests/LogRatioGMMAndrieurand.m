function LogRatio = LogRatioGMMAndrieurand(x,xstar,Parameters)

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
    
load('TempInd.mat')
temp = post(ind);
LogqTempStar = log(temp);
LogLikxstar = log(max(eps^3,Parameters.f(xstar,Parameters)));


% Star to Temp 
temp = 0;
if strcmp(Parameters.Mode,'Log')
    post = log(1+posts(2,:))/sum(log(1+posts(2,:)));
elseif strcmp(Parameters.Mode,'L')
    post = posts(2,:)/sum(posts(2,:));
end


temp = post(ind);
LogqStarTemp = log(temp);
LogLikx = log(max(eps^3,Parameters.f(x,Parameters)));

LogNum = LogLikxstar + LogqStarTemp;
LogDenom = LogLikx + LogqTempStar;
LogRatio = LogNum - LogDenom;

if isnan(LogRatio)
    disp('Ratio GMM rand problem')
end