function [inf,sup] = Quantiles95(Values,Weigths)

    if nargin == 1
        Weigths = ones(length(Values),1);
    end

    temp = [];
    for i = 1:length(Values)
        temp = [temp; Values(i)*ones(max(round(Weigths(i)*1000),1),1)];
    end
    sup = quantile(temp,0.975);
    inf = quantile(temp,0.225);
    