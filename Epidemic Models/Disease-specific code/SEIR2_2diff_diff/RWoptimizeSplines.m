function x = RWoptimizeSplines(x,Data,Parameters,N,Sigma)

TempLogLik = SplinePrevLik(x,Data,Parameters);


acc = [];
LogLiks = [];
for i = 1:N
    xstar = x+Sigma*rand(size(x));
    StarLogLik = SplinePrevLik(xstar,Data,Parameters);
    
    acc(i)=0;
    if StarLogLik<TempLogLik
        x = xstar;
        TempLogLik = StarLogLik;
        acc(i)=1;
    end
    LogLiks(i) = TempLogLik;
    disp(mean(acc))
end
clf
plot(LogLiks)
disp(['AccRate: ' num2str(mean(acc))]);
