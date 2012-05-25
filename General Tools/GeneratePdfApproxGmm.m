function GmmPdfApprox = GeneratePdfApproxGmm(x)

% x should be 1-d data.
Liks = [];
BICs = [];
objs = {};

for NbComp = 1:4
    obj = gmdistribution.fit(x',NbComp);
    Liks(NbComp) = obj.NlogL ;
    BICs(NbComp) = obj.BIC;
    objs{NbComp} = obj;
end

[bof,Ind] = min(BICs);
disp([num2str(Ind) ' components'])
obj = objs{Ind};

mult = 100;
res = random(obj,length(x)*mult);
[ord,abs] = hist(res);
hist(x)
hold on
plot(abs,ord/mult,'r','LineWidth',2)
hold off

GmmPdfApprox = obj;
