function ParametersStar = SEIR_BuildPathConstantBeta(ParametersStar,Data)

% should give the path of log(beta)
indbeg = ceil(ParametersStar.Begining.Value);
Path = ones(1,sum(Data.NbComputingSteps));
Path(1:sum(Data.NbComputingSteps(1:indbeg))) = ParametersStar.FirstBeta.TransfValue*ones(size(sum(Data.NbComputingSteps(1:indbeg))));
Path(sum(Data.NbComputingSteps(1:indbeg))+1:sum(Data.NbComputingSteps(1:9))) = ParametersStar.SecondBeta.TransfValue*ones(size(sum(Data.NbComputingSteps(1:9))));
Path(sum(Data.NbComputingSteps(1:9))+1:sum(Data.NbComputingSteps(1:15))) = ParametersStar.ThirdBeta.TransfValue*ones(size(sum(Data.NbComputingSteps(1:9)+1:sum(Data.NbComputingSteps(1:15)))));
Path(sum(Data.NbComputingSteps(1:15))+1:sum(Data.NbComputingSteps(1:end))) = ParametersStar.FourthBeta.TransfValue*ones(size(sum(Data.NbComputingSteps(1:15)+1:sum(Data.NbComputingSteps(1:end)))));

try
    ParametersStar.Pars = [Path'; ParametersStar.Pars];
catch
    ParametersStar.Pars = [Path'; ParametersStar.Pars'];
end
ParametersStar.LengthPaths = sum(Data.NbComputingSteps);
