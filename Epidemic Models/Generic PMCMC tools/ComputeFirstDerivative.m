function res = ComputeFirstDerivative(Parameters)

X = Parameters.Pars;

DensityModel = Parameters.RealScaleDensityModel;

NbComps = DensityModel.NComponents;

dfdx = zeros(size(X));
for i = 1:NbComps
    sigma =  squeeze(DensityModel.Sigma(:,:,i));
    mu = DensityModel.mu(i,:)';
    dfdx = dfdx - 1/2*DensityModel.PComponents(i)*((sigma^(-1))'*(X-mu)+(sigma^(-1))*(X-mu))*mvnpdf(X,mu,sigma);
end

res = dfdx/pdf(DensityModel,X');
    
