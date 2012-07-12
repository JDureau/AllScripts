function Grad = ComputeGMMGrad(x,Parameters)


Mixture = Parameters.Dens;
ws = Parameters.Dens.PComponents;

Grad = zeros(Parameters.Dim,1);
tmp = 0;
Sigmas = Mixture.Sigma;
Mus = Mixture.mu;
for i = 1:Mixture.NComponents
    Sigma = squeeze(Sigmas(:,:,i));
    Mu = squeeze(Mus(i,:))';
    Grad = Grad - ws(i)*(Sigma^-(1))*(x-Mu)*mvnpdf(x,Mu,Sigma);
end
Grad = Grad/pdf(Parameters.Dens,x');


    
    
    
    
    