function xstar = SampleGMMind(x,Parameters)

Mixture = Parameters.Dens;

xstar = random(Mixture);

