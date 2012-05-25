function res = rank(v) 

[bof,inds] = sort(-v);
res = [];
for i = 1:length(v)
    res(inds(i)) = i;
end