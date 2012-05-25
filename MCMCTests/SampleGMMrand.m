function xstar = SampleGMMrand(x,Parameters)

Mixture = Parameters.Dens;


liks = posterior(Mixture,x);
if strcmp(Parameters.Mode,'Log')
    ind = PickRandInd(log(1+liks));
elseif strcmp(Parameters.Mode,'L')
    ind = PickRandInd(liks);
end
    
    
mu = Parameters.Mus(ind,:)';
Epsil = Parameters.Epsil;


xstar = mvnrnd(x',squeeze(Epsil^2*Parameters.Sigmas(:,:,ind)));

if isnan(xstar)
    disp('Sample GMM rand problem')
end
