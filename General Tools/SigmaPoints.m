function[Xi W] = SigmaPoints(xm,P,kappa)

n = numel(xm);
Xi = zeros(n,2*n+1);
W = zeros(2*n+1,1);

Xi(:,1) = xm;
W(1) = kappa / (n+kappa);
U = chol((n+kappa)*P);

for k = 1:n
    Xi(:,k+1) = xm +U(k,:)';
    W(k+1) = 1/(2*(n+kappa));
end

for k = 1:n
    Xi(:,k+n+1) = xm -U(k,:)';
    W(k+n+1) = 1/(2*(n+kappa));
end



