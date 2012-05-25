function SamplesOut = GenerateCopulaSamples(SamplesIn)

U = zeros(size(SamplesIn));
for i = 1:size(SamplesIn,2)
    i
    U(:,i) = ksdensity(SamplesIn(:,i),SamplesIn(:,i),'function','cdf');
end

[Rho,nu] = copulafit('t',U,'Method','ApproximateML');

U = copularnd('t',Rho,nu,size(SampelsIn,1));

SamplesOut = zeros(size(SamplesIn));
for i = 1:size(SampelsIn,2)
    SamplesOut(:,i) = ksdensity(U(:,i),U(:,i),'function','cdf');
end
