function LogLik = SEIR_ComputeLogLik(Data,Model,Parameters)

Parameters.LengthPaths = size(Parameters.Paths,2);
Parameters.Pars = Parameters.Paths(6,:);
Temp = SEIR_PathLik(Data, Model,Parameters);
LogLik = Temp.LogLik;
