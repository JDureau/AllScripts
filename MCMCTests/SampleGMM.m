function xstar = SampleGMM(x,Parameters)

Mixture = Parameters.SamplDens;


liks = posterior(Mixture,x);
ind = PickRandInd(liks);

mu = Parameters.Mus(ind,:)';
Epsil = Parameters.Epsil;


xstar = mvnrnd(x'- Epsil^2/2*(x'-mu),Epsil^2*Parameters.Sigmas);
