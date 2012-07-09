function LogRatio = LogRatioGMMLang(x,xstar,Parameters)

Epsil = Parameters.Epsil;
Dens = Parameters.Dens;
Nc = Parameters.Nc;
Mus = Parameters.Mus;
Sigmas = Parameters.Sigmas;

posts = posterior(Dens,[x ; xstar]);

% numerator (xstar)
temp = 0;
post = posts(1,:);
if strcmp(Parameters.Mode,'Log')
    post = log(1+posts(1,:))/sum(log(1+posts(1,:)));
elseif strcmp(Parameters.Mode,'L')
    post = posts(1,:)/sum(posts(1,:));
end
for IndComp = 1:Nc
    mu = Mus(IndComp,:)';
    Sigma = squeeze(Sigmas(:,:,IndComp));
    temp = temp + (mvnpdf(xstar',x'-Epsil^2*(x'-mu)/2,Parameters.Epsil^2*Sigma)*post(IndComp));
end
LogqTempStar = log(temp);
LogLikxstar = log(Parameters.f(xstar,Parameters));


% Star to Temp 
temp = 0;
post = posts(2,:);
if strcmp(Parameters.Mode,'Log')
    post = log(1+posts(2,:))/sum(log(1+posts(2,:)));
elseif strcmp(Parameters.Mode,'L')
    post = posts(2,:)/sum(posts(2,:));
end
for IndComp = 1:Nc
    mu = Mus(IndComp,:)';
    Sigma = squeeze(Sigmas(:,:,IndComp));
    temp = temp + (mvnpdf(x',xstar'-Epsil^2*(xstar'-mu)/2,Parameters.Epsil^2*Sigma)*post(IndComp));
end
LogqStarTemp = log(temp);
LogLikx = log(Parameters.f(x,Parameters));

LogNum = LogLikxstar + LogqStarTemp;
LogDenom = LogLikx + LogqTempStar;
LogRatio = LogNum - LogDenom;

if isnan(LogRatio)
    disp('Problem')
end