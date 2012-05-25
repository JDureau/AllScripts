function error = TestNormalModel(x,xvalues,yvalues,MinVect,MaxVect,MinC,MaxC,ampl,inds,weigths)

dim = size(xvalues,2);
if nargin == 9
    weigths = ones(size(yvalues));
end

c = InvLogitTransf(x(1),MinC,MaxC);
mu = InvLogitTransf(x(2:1+dim)',MinVect',MaxVect');

a = x(dim+2:end); 
% a = InvLogitTransf(a,MinBounds,MaxBounds);
b = triu(ones(dim),1)+eye(dim); 
b(~~b)=a ;
sigma = b*b';

temp = min(1000,mvnpdf(xvalues(inds,:),mu',sigma));
% error = log((mean(weigths.*((c*yvalues-mvnpdf(xvalues,mu,sigma)')./(c*yvalues)).^2)));

error = (sum(weigths(inds).*((c*exp(yvalues(inds))-temp')./(max(eps,temp))').^2))*(1 + sum(max(0,sqrt(diag(sigma))-ampl')));
if isinf(error)
    disp('inf value')
else
%     disp((1 + sum(max(0,sqrt(diag(sigma))-ampl'))))
end
