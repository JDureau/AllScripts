function [AutoCor] = AutoCorrelation( X, MaxLag )

if nargin == 1
    MaxLag = length(X)/20;
else
    MaxLag = min(MaxLag,length(X)/20);
end

X = (X-mean(X))/std(X);

AutoCor = [];
test = 0;
k = 0;
while not(test)
    k = k+1;
    n = length(X)-k+1;
    AutoCor(1,end+1)=1/n*(X(k:length(X))*(X(1:length(X)-k+1)'));
    test = 0;
    if k > MaxLag
        test = 1;
    end
    if k>300
        if sum(AutoCor(min(end-300,round(end/2)):end))<0.05
            test = 1;
        end
    end
end