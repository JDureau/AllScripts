function Y = SampleObs(X,step,sigma)

N = length(X);
nobs = N*step;
npoints = N/nobs;


Y = zeros(npoints+1,1);
for i = 1:nobs
    Y(i+1) = sqrt(sum(sigma(X((i-1)*npoints+1:i*npoints)).^2))*randn(1,1);   
end
    