function res = MyBIC(GMMobj,samples)

dim = GMMobj.NDimensions;
ncomps = GMMobj.NComponents;
k = ncomps-1 + dim*ncomps + ncomps*dim*(dim+1)/2;
res = -2*sum(log(pdf(GMMobj,samples))) + k*log(max(size(samples)));