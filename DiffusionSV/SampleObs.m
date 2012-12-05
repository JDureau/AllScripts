function Y = SampleObs(X,step,sigma)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


N = length(X)+1;
nobs = N*step+1;
npoints = N/(nobs-1);


Y = zeros(nobs,1);
Y(1) = 0;
Y(2) = Y(1) + sqrt(sum(sigma(X(1:npoints-1)).^2))*randn(1,1)*sqrt(step);  % accounting for the fact that X(0) = 0 so we only take first npoints-1 elements of X 
for i = 3:nobs
    Y(i) = Y(i-1) + sqrt(sum(sigma(X((i-2)*npoints:(i-1)*npoints-1)).^2))*randn(1,1)*sqrt(step);   
end
    