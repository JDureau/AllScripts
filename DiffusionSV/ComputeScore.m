function Score = ComputeScore(Z,Y,sigma,dersigma,H,ind)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;

N = length(Z)/2+1;
nobs = length(Y);
step = (nobs-1)/N;
npoints = N/(nobs-1);

Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh);

Sigmas_s = sigma(X);
Sigmas_s_prime = dersigma(X);

% a = dlog(Y|X)/dX
a = zeros(N-1,1);
for j = 1:N-1
    k = floor(j/npoints)+1;
    if k == 1
        denom = sum(Sigmas_s(1:k*npoints-1).^2)*step;
    else
        denom = sum(Sigmas_s((k-1)*npoints:k*npoints-1).^2)*step;
    end
    a(j) = -Sigmas_s_prime(j)*Sigmas_s(j)*step/denom;
    a(j) = a(j) + (Y(k+1)-Y(k-1+1))^2*Sigmas_s_prime(j)*Sigmas_s(j)*step/(denom^2);
end

% B = dX/d\DeltaB_H
B = zeros(N-1,N-1);
for j = 1:N-1
    for i = j:N-1
        B(i,j) = 1;
    end
end

% d\DeltaB_H/dZ 
M = ComputeM(N);
Lambda = ComputeLambda(N,step,H);
SqrtLambda = zeros(2*(N-1),2*(N-1));
for j = 1:2*(N-1)
    SqrtLambda(j,j) = sqrt(real(Lambda(j,j)));
end
Extr = zeros(N-1,2*(N-1));
for j = 1:(N-1)
    Extr(j,j) = 1;
end
P = zeros(2*(N-1),2*(N-1));
for j = 1:2*(N-1)
    for k = 1:2*(N-1)
        P(j,k) = 1/sqrt(2*(N-1))*exp(-2*pi*1i*(j-1)*(k-1)/(2*(N-1)));
    end
end

C = Extr*P*SqrtLambda*M;

% Score = a(ind);

% Score = 0;
% for k = 1:N-1
%     Score = Score + a(k)*B(k,ind);
% end


tmp = a'*B*Extr;
tmp = ifft(tmp);
tmp = real(tmp*SqrtLambda*M);

Score = tmp(ind);
% 
% Score = 0;
% for k = 1:N-1
%     for i = 1:N-1
%         Score = Score + a(k) * B(k,i) * real(C(i,ind));
%     end
% end

% Score = a'*(B');
% Score = Score*(C);
