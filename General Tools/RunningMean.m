function res = RunningMean(x,k)

res = zeros(size(x));
for i = 1:length(x)
    res(i) = mean(x(max(1,i-k):min(length(x),i+k)));
end