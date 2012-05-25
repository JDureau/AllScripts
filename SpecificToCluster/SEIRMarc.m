function [] = SEIRMarc(ind)


cd('/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllScripts/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR'])


SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc';
Name = [SavePath '/MarcData_Add_LogNorm_Tau2.mat'];

load(Name)

Data = Res3.Data;
Parameters = Res3.Parameters;
SEIRModel = Res3.Model;


TempRes = Res3;

Cov = 2.38^2/dim*cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
TempPar = TempRes.TempPar;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);

% if strcmp(Parameters.DiffusionType,'IBM')
%     Parameters.NoPaths = 0;
% else
%     Parameters.NoPaths = 0;
% end
% Parameters.PathsToKeep = [1:7]';
if ind < 20
    Res3 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,5000);
    Res3.Description = 'This one has unif prior on days for latent and inf periods. It also uses the actual averaged data and not only >65 as FirstEst';
    save([Name num2str(ind)],'Res3')
end
