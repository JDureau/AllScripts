function ScaledThetas = Scale(Thetas,Parameters)


ScaledThetas = [];

for i = 1:length(Parameters.ScalingMeans)
    ScaledThetas(i,:) = (Thetas(i,:) - Parameters.ScalingMeans(i))./Parameters.ScalingStds(i);
end