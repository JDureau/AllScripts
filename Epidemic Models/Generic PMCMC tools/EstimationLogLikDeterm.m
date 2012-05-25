function Res = EstimationLogLikDeterm(Data,Model,ParametersStar)

ParametersStar = Model.BuildPath(ParametersStar,Data); % build the path from pars and put it into Parameters.Path
Res = Model.PathLik(Data, Model,ParametersStar); % call a function like SEIR_PathLik

Res.Grad = [];
Res.Cov = ParametersStar.G^-1;
Res.G = ParametersStar.G;
Res.CompInd = [];