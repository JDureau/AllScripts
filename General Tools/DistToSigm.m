function d = DistToSigm(x,data,Parameters)

xis = (1:length(data))*Parameters.ComputationTStep;
sigm = x(1)+x(2)*Sigmoid((xis(1:length(data))-x(3)*Parameters.ComputationTStep)/x(4));

d = mean((sigm-data).^2);