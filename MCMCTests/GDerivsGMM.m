function GDerivs = GDerivsGMM(x,Parameters)

GDerivs = {};

for i = 1:length(x)
    GDerivs{i} = zeros(length(x),length(x));
end
