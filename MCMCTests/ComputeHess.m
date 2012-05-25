function Parameters = ComputeHess(x,Parameters)

  
% Grad and Hess
xs = [];
hs = [];
Grad = [];
Hess = [];

for i = 1:length(x)
    Temp = log(Parameters.f(x,Parameters));
    fxs(i) = Temp;
    hs(i) = 0.01 ;
    xphs =  x;
    xmhs =  x;
    xphs(i) =  x(i) + hs(i);
    xmhs(i) =  x(i) - hs(i);
    Temp = log(Parameters.f(xphs,Parameters));
    fxphs(i) = Temp;
    Temp = log(Parameters.f(xmhs,Parameters));
    fxmhs(i) = Temp;
    
    Grad(i) =   ( fxphs(i) - fxmhs(i) )/ (2*hs(i));
    Hess(i,i) = ( fxphs(i) - 2*fxs(i) + fxmhs(i) )/ (hs(i)^2);

end



for i = 1:length(x)-1
    for j = i+1:length(x)
        % x+hx z+hy
        xpphs = x;
        xpphs(i) = x(i) + hs(i);
        xpphs(j) = x(j) + hs(j);
        Temp = log(Parameters.f(xpphs,Parameters));
        fxphzph = Temp;
        % x+h z-h
        xpmhs = x;
        xpmhs(i) = x(i) + hs(i);
        xpmhs(j) = x(j) - hs(j);
        Temp = log(Parameters.f(xpmhs,Parameters));
        fxphzmh = Temp;
        % x-h z+h
        xmphs = x;
        xmphs(i) = x(i) - hs(i);
        xmphs(j) = x(j) + hs(j);
        Temp = log(Parameters.f(xmphs,Parameters));
        fxmhzph = Temp;
        % x-h z-h
        xmmhs = x;
        xmmhs(i) = x(i) - hs(i);
        xmmhs(j) = x(j) - hs(j);
        Temp = log(Parameters.f(xmmhs,Parameters));
        fxmhzmh = Temp;

        Hess(i,j) = (fxphzph - fxphzmh - fxmhzph + fxmhzmh)/(4*hs(i)*hs(j));
        Hess(j,i) = (fxphzph - fxphzmh - fxmhzph + fxmhzmh)/(4*hs(i)*hs(j));    
    end
end

Parameters.Hess = Hess;