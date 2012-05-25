function [AutoCor] = AutoCorrelation( X, MaxLag )

if nargin == 1
    MaxLag = length(X);
end

X = (X-mean(X))/std(X);

AutoCor = zeros(1,MaxLag);
test = 0;
k = 0;
while not(test)
    k = k+1;
    n = length(X)-k+1;
    AutoCor(1,k)=1/n*(X(k:length(X))*(X(1:length(X)-k+1)'));
    test = 1;
    if k < MaxLag
        test = 0;
    elseif and(isempty(find(AutoCor<0.05)),k < length(X)/2)
        test = 0;
    end
end