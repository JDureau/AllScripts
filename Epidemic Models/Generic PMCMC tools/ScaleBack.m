function Thetas = ScaleBack(ScaledThetas,Parameters)

Thetas = [];

for i = 1:length(Parameters.ScalingMeans)
    Thetas(i,:) = ScaledThetas(i,:)*Parameters.ScalingStds(i) + Parameters.ScalingMeans(i);
end