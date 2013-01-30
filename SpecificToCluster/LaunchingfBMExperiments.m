function [] = LaunchingfBMExperiments(ind,DataSet)

ind = ind+1;

s = RandStream('mcg16807','Seed',ind);
RandStream.setDefaultStream(s)


cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/DiffusionSV'])

SavePath = '/users/ecologie/dureau/src/AllData/fBM/';

loop = 100;

if DataSet == 1
    load([SavePath '/DataSet1.mat'])
else
    load([SavePath '/DataSet2.mat'])
end

Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative


Par = Data.ParTrue;

switch ind

    case 1
        Par.theta_sampler = 'GibbsHMC';
        Par.nsteps = 1;
        if DataSet == 1
            Par.hZ = 1;
            Par.hP = 0.11;
        else
            Par.hZ = 0.1;
            Par.hP = 0.04;
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(ind)],'Res')
    case 2
        Par.theta_sampler = 'GibbsHMC';
        Par.nsteps = 10;
        if DataSet == 1
            Par.hZ = 1;
            Par.hP = 0.11;
        else
            Par.hZ = 0.1;
            Par.hP = 0.04;
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(ind)],'Res')
    case 3
        Par.theta_sampler = 'GibbsHMC';
        Par.nsteps = 20;
        if DataSet == 1
            Par.hZ = 1;
            Par.hP = 0.11;
        else
            Par.hZ = 0.1;
            Par.hP = 0.04;
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(ind)],'Res')
    case 4
        Par.theta_sampler = 'JointHMC';
        Par.nsteps = 1;
        if DataSet == 1
            Par.h = 0.1;
        else
            Par.h = 0.035;
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(ind)],'Res')
    case 5
        Par.theta_sampler = 'JointHMC';
        Par.nsteps = 10;
        if DataSet == 1
            Par.h = 0.1;
        else
            Par.h = 0.035;
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(ind)],'Res')
    case 6
        Par.theta_sampler = 'JointHMC';
        Par.nsteps = 20;
        if DataSet == 1
            Par.h = 0.1;
        else
            Par.h = 0.035;
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(ind)],'Res')
    case 7
        if DataSet == 1
            Par.Epsil = 1;
            Par.MCMCType = 'Rand';
            Par.G = Data.Cov^(-1);
            Par.ModelType='SMC';
            Par.NbVariables = 3;
            Par.NbParticules = 100;
            Par.NoPaths = 0;
            Par.PathsToKeep = [1];
            Par.NbParsEstimated  =length(Par.Names.Estimated);
            Par.ComputationTStep = Data.step;
            Par.Vol = Vol;
            Par.Problem = 'vol';
            Data.ObservedVariables = 1;
            Par.AdaptC = 0.999;
            Par.GMeth =  'cst given';
            Data.NbComputingSteps = [0 Data.obsstep*ones(1,Data.nobs)] ;
            fullvolModel.InitializeParameters = @fullvolInitialize;
            fullvolModel.SMC_projection = @fullvol_SMC_projection;
            fullvolModel.LikFunction = 'normpdf(Data.Y(IndTime)-Data.Y(IndTime-1),Variables(:,2),Variables(:,3))';
            TempPar = ProposeInitialParameter(Data, fullvolModel, Par);
            Res = RunEstimationMethod(Data, fullvolModel, Par, TempPar, loop);
            save([SavePath '/D' num2str(DataSet) '_Exp' num2str(ind)],'Res')
        end
end
    

