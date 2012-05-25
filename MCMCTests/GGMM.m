function G = GGMM(x,Parameters)

G = zeros(Parameters.Dim,Parameters.Dim);

if and(size(x,1) == 1,size(x,2) == Parameters.Dim)
    p = pdf(Parameters.RealDens,x);
    x = x';
elseif and(size(x,2) == 1,size(x,1) == Parameters.Dim)
    p = pdf(Parameters.RealDens,x');
else
    'Wrong argument in fGMM'
    die
end

gammas = Parameters.Dens.PComponents;
mus = Parameters.Dens.mu;
exps = [];
Sm1s = {};
Ns = [];
for i = 1:Parameters.Dens.NComponents
    Ns(i) = (2*pi)^(Parameters.Dim/2)*sqrt(det(squeeze(Parameters.Dens.Sigma(:,:,i))));
    Sm1s{i} = squeeze(Parameters.Dens.Sigma(:,:,i))^(-1);
    exps(i) = exp(-0.5*(x-mus(i,:)')'*Sm1s{i}*(x-mus(i,:)'));
end

A = zeros(Parameters.Dim,Parameters.Dim);
for i = 1:Parameters.Dens.NComponents
    A = A + 1/p*(-gammas(i)/(2*Ns(i))*(Sm1s{i}*exps(i)-0.5*Sm1s{i}*(x-mus(i,:)')*(Sm1s{i}*(x-mus(i,:)'))'*exps(i)));
end
B = zeros(Parameters.Dim,Parameters.Dim);
for i = 1:Parameters.Dens.NComponents
    for j = 1:Parameters.Dens.NComponents
        B = B + 1/(p^2)*(-gammas(i)/(2*Ns(i))*(Sm1s{i}*(x-mus(i,:)')*exps(i)))*(-gammas(j)/(2*Ns(j))*(Sm1s{j}*(x-mus(j,:)')*exps(j))');
    end
end

G = -A - B;

G
eig(G)
det(G)
% 
% 'pause'