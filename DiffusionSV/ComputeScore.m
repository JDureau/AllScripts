function Score = ComputeScore(Z,Y,sigma,dersigma,H)

N = length(Z)/2;
nobs = length(Y);
step = (nobs-1)/N;
npoints = N/(nobs-1);

Bh = Z_to_Bh(Z,N,step,H);
X = Bh_to_X(Bh);

Sigmas_s = sigma(X);
Sigmas_s_prime = dersigma(X);

% a = dlog(Y|X)/dX
a = zeros(N,1);
for j = 1:N
    k = floor((j-1)/npoints)+1;
    denom = sum(Sigmas_s((k-1)*npoints+1:k*npoints).^2)*step;
    a(j) = -Sigmas_s_prime(j)*Sigmas_s(j)*step/denom;
    a(j) = a(j) + (Y(k+1)-Y(k-1+1))^2*Sigmas_s_prime(j)*Sigmas_s(j)*step/(denom^2);
end

% B = dX/d\DeltaB_H
B = zeros(N,1);
for j = 2:N
    B(j,j) = 1;
end

% d\DeltaB_H/dZ 
M = ComputeM(N);
Lambda = ComputeLambda(N,step,H);
SqrtLambda = zeros(2*N,2*N);
for j = 1:2*N
    SqrtLambda(j,j) = sqrt(real(Lambda(j,j)));
end
Extr = zeros(N,2*N);
for j = 1:N
    Extr(j,j) = 1;
end
P = zeros(2*N,2*N);
for j = 1:2*N
    for k = 1:2*N
        P(j,k) = 1/sqrt(2*N)*exp(-2*pi*1i*(j-1)*(k-1)/(2*N));
    end
end

C = Extr*P*SqrtLambda*M;


% Score = a(1);

Score = 0;
j = 1;
for k = 1:N
    for i = 1:N
        Score = Score + a(k) * B(k,i) * real(C(i,j));
    end
end

% Score = a'*(B');
% Score = Score*(C);
