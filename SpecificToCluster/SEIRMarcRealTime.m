function [] = SEIRMarcRealTime(ind,length,inddiff)
s = RandStream('mcg16807','Seed',1000*ind)
RandStream.setDefaultStream(s)

ind = ind-floor(ind/6)*6;


if inddiff == 1
  typediff = 'Add'
else
  typediff = 'IBM'
end

%rng('shuffle')

cd('/users/ecologie/dureau/src/AllScripts/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR'])


SavePath = '/users/ecologie/dureau/src/AllData/ResultsMarc/';
Name = [SavePath '/TestsRealTime_' typediff '_' num2str(length) 'weeks_NoPaths.mat'];

load(Name)

dim = 8;
Data = Res3.Data;
Parameters = Res3.Parameters;
SEIRModel = Res3.Model;

%Parameters.ComputationTStep = 0.1;
%Data.Instants = [0:size(Data.Observations,2)-1]*7/Parameters.ComputationTStep;
%Data.ObservedVariables = 5*ones(1,length(Data.Instants));
%Data.NbComputingSteps = [0 diff(Data.Instants)];

TempRes = Res3;

Cov = 2.38^2/dim*cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
TempPar = TempRes.TempPar;
Parameters.NbParticules = 2000;
Parameters.AdaptC = 0.98;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);

% if strcmp(Parameters.DiffusionType,'IBM')
%     Parameters.NoPaths = 0;
% else
%     Parameters.NoPaths = 0;
% end
% Parameters.PathsToKeep = [1:7]';
n = 15000;
%Name = [SavePath 'Sim0p1_Alts' num2str(ceil(ind/10)) 'Cluster'];
Name = [SavePath 'RealTime_' typediff '_' num2str(length) '_cluster' num2str(ind) '.mat'];

if ind > 3
  Parameters.NoPaths = 0;
  Parameters.PathsToKeep = [1:7]';
  n = 1000;
end
Res3 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,n);
Res3.Description = 'This one has unif prior on days for latent and inf periods. It also uses the actual averaged data and not only >65 as FirstEst';
save([Name],'Res3')

