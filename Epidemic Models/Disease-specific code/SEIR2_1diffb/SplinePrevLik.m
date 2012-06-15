function res = SplinePrevLik(x,Data,Parameters)


ts = ceil(x(1:length(x)/2));
[ts,I,J] = UNIQUE(ts);

yis = x(length(x)/2+1:end);
yis = yis(I);

n = sum(Data.NbComputingSteps);
xbetas = 1:n;

BetaFit = spline(ts,yis,xbetas);


SEIRModel = struct();
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);

Parameters.ObsNoise = 0;
SimData = SEIR_CreateData(BetaFit,Parameters,Data,SEIRModel);

LogLik = 0;
for i = 2:size(Data.Observations,2)
    LogLik = LogLik + log(normpdf(SimData.Observations(5,i),Data.Observations(5,i),Data.Observations(5,i)*Parameters.SigmaObs.Value));
end
res = -LogLik;
% disp(res)
