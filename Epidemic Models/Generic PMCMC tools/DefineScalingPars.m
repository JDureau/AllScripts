function Parameters = DefineScalingPars(ThetasSamples, Parameters)

Parameters.ScalingMeans = mean(ThetasSamples');
Stds = [];
for i = 1:length(Parameters.ScalingMeans)
    Stds(i) = std(ThetasSamples(i,:) - Parameters.ScalingMeans(i));
end
Parameters.ScalingStds = Stds;
