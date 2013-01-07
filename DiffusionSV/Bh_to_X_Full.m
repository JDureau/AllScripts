function X = Bh_to_X_Full(B,H,step,sigma_X,kappa)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i-1) with X(0) = 0;


X = zeros(length(B),1); 
for i = 1:length(B)-1
    X(i+1,1) = X(i,1)*(1+kappa*step) + sigma_X*B(i);
%     X(i+1) = X(i) + kappa * X(i) * step + sigma_X*(step^H)*B(i);
end
    
    
    
    