function Grad = ComputeGMMGrad(x,Parameters)


Mixture = Parameters.Dens;
ws = Parameters.Dens.PComponents;

Grad = zeros(Parameters.Dim,1);
tmp = 0;
for i = 1:Mixture.NComponents
    Sigma = squeeze(Mixture.Sigma(:,:,i));
    Mu = squeeze(Mixture.mu(i,:))';
    Grad = Grad - ws(i)*(Sigma^-(1))*(x-Mu)*mvnpdf(x,Mu,Sigma);
end
Grad = Grad/pdf(Parameters.Dens,x');


    
    
    
    
    