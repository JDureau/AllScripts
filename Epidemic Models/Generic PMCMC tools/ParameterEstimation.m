function Result = ParameterEstimation(Data,Parameters)




% Initialise TempPar
Parameters.NoLambdas = 1;
if strcmp(Parameters.MCMCType,'Rand')
    Parameters.NoDerivatives = 1;
else
    Parameters.NoDerivatives = 0;
end
TempPar = ProposeInitialParameter(Data, Parameters);

% Calibrate the method
if Parameters.NbParsEstimated>1
    Parameters.ComputeRWsamplCov = 1;
else
    Parameters.ComputeRWsamplCov = 0;
end
Parameters.Calibrate = {};
[Parameters, TempPar] = CalibrateMethod( Data, Parameters, TempPar);

% Burn-In
Res = RunEstimationMethod(Data,Parameters,TempPar,500);


TempPar = Res.TempPar;

% Calibrate the method after BurnIn
[Parameters, TempPar] = CalibrateMethod( Data, Parameters, TempPar);


% Run The algorithm
Parameters.NoLambdas = 0;
Res = RunEstimationMethod(Data,Parameters,TempPar,Parameters.NbIterations);

Result = Res;
Result.Parameters = Parameters;


