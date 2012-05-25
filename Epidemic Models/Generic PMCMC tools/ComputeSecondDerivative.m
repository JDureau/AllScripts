function res = ComputeSecondDerivative(Parameters)

X = Parameters.Pars;

DensityModel = Parameters.RealScaleDensityModel;

NbComps = DensityModel.NComponents;

d2fdx2 = zeros(length(X),length(X));
for i = 1:NbComps
    sigma =  squeeze(DensityModel.Sigma(:,:,i));
    mu = DensityModel.mu(i,:)';
    d2fdx2 = d2fdx2 - 1/2*DensityModel.PComponents(i)*((sigma^(-1))'+(sigma^(-1)))*mvnpdf(X,mu,sigma) ...
        + 1/4*DensityModel.PComponents(i)*(((sigma^(-1))'+(sigma^(-1)))*(X-mu))*(((sigma^(-1))'+(sigma^(-1)))*(X-mu))'*mvnpdf(X,mu,sigma);
end

f = pdf(DensityModel,X');
dfdx = ComputeFirstDerivative(Parameters)*f;

res = (d2fdx2*f-dfdx*dfdx')/f^2;
    
