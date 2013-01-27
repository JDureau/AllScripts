function Y = SampleObs_Full(X,B,step,sigma,Par)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


N = length(X);
nobs = N*step+1;
npoints = N/(nobs-1);
rho = Par.rho.Value;
mu_Y = Par.mu_Y.Value;

Y = zeros(nobs,1);
Y(1) = 0;
for i = 2:nobs
    Y(i) = Y(i-1) + sqrt(1-rho^2)*sqrt(sum(sigma(X((i-2)*npoints+1:(i-1)*npoints)).^2))*randn(1,1)*sqrt(step);  
    Y(i) = Y(i)   + sum((mu_Y-sigma(X((i-2)*npoints+1:(i-1)*npoints)).^2/2)*step) ;
    Y(i) = Y(i)   + rho*sum(sigma(X((i-2)*npoints+1:(i-1)*npoints)).*B((i-2)*npoints+1:(i-1)*npoints));
end

