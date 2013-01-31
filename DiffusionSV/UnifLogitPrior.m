function prior = UnifLogitPrior(Name,Parameters,TransfValues,Values)
x1 = Parameters.(Name).MinLim;
x2 = Parameters.(Name).MaxLim;

% prior = 1/(x2-x1);

delta = 0.10*(x2-x1);


alpha = 1/(x2-x1-2/3*delta);

a = -alpha/(delta^2);
b = -2*(x1+delta)*a;
c = alpha+a*(x1+delta)^2;

if nargin == 2
    x = Parameters.(Name).Value;
%     x = min(x,x2-0.0000001*delta);
%     x = max(x,x1+0.0000001*delta);
    if x<x1+delta
        tmp = a*x^2 + b*x + c;
    elseif x<x2-delta
        tmp = alpha;
    else
        tmp = a*(x-(x2-x1-2*delta))^2 + b*(x-(x2-x1-2*delta)) + c;
    end
    prior = tmp;     
    
else 
    Values = min(Values,x2-eps*delta);
    Values = max(Values,x1+eps*delta);
    inds1 = find(Values<a+delta);
    inds2 = find(and(Values>x1+delta,Values<x2-delta));
    inds3 = find(Values>=x2+delta);
    tmp = zeros(size(Values));
    tmp(inds1) = a*Values(inds1).^2  + b*Values(inds1) + c;
    tmp(inds2) = alpha*ones(size(inds2));
    tmp(inds3) = a*(Values(inds3)-(x2-x1-2*delta)*ones(size(inds3))).^2  + b*(Values(inds3)-(x2-x1-2*delta)*ones(size(inds3))) + c ;
    
    prior = tmp;     
end

prior = max(tmp,eps);
