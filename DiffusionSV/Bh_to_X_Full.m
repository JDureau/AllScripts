function X = Bh_to_X_Full(B,step,Par)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i-1) with X(0) = 0;


H = Par.H.Value;
sigma_X = Par.sigma_X.Value;
kappa = Par.kappa.Value;
X0 = Par.X0.Value;
mu_X = Par.mu_X.Value;
    
X = zeros(length(B),1); 
X(1) = X0;
for i = 1:length(B)-1
    X(i+1,1) = X(i,1)*(1-kappa*step) + kappa*mu_X*step + sigma_X*B(i);
%     X(i+1) = X(i) + kappa * X(i) * step + sigma_X*(step^H)*B(i);
end
    
    
