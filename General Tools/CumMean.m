function res = CumMean(v)

res = [];
for i = 1:length(v)
    res(i) = mean(v(1:i));
end