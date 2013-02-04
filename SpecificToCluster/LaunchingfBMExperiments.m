function [] = LaunchingfBMExperiments(ind,DataSet,OptionCov,OptionS,indExp)

ind = ind+1;

s = RandStream('mcg16807','Seed',indExp+1);
RandStream.setDefaultStream(s)


cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/DiffusionSV'])

SavePath = '/users/ecologie/dureau/src/AllData/fBM/';

loop = 20000;

if DataSet == 3
    load([SavePath '/DataSet3.mat'])
else
    load([SavePath '/DataSet4.mat'])
end

Vol = @ClassicVol; % how the volatility X plays on the price
VolDer = @DerClassicVol; % its derivative


Par = Data.ParTrue;


if not(OptionCov)
    Data.Cov = eye(length(Par.Names.Estimated));
    Data.CovS = eye(length(Par.Names.Estimated)-1);
end
if OptionS
    Par.H.Estimated = 0;
    Par = DefineIndexes(Par);
    Par = NoTransfToTransf(Par);
    Data.Cov = Data.CovS;
end

switch ind

    case 1
        Par.theta_sampler = 'GibbsHMC';
        Par.nsteps = 1;
        if not(OptionCov)
            if DataSet == 3
                Par.hZ=0.6;
                Par.hP=0.025;
            else
                Par.hZ=0.18;
                Par.hP=0.025;%
            end
        else
            if DataSet == 3
                Par.hZ=0.6;
                Par.hP=0.18;
            else
                Par.hZ=0.2;
                Par.hP=0.12;%
            end
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(indExp) '_' num2str(ind) '_' num2str(OptionS) '_' num2str(OptionCov)],'Res')
    case 2
        Par.theta_sampler = 'GibbsHMC';
        Par.nsteps = 10;
        if not(OptionCov)
            if DataSet == 3
                Par.hZ=0.6;
                Par.hP=0.025;
            else
                Par.hZ=0.18;
                Par.hP=0.025;%
            end
        else
            if DataSet == 3
                Par.hZ=0.6;
                Par.hP=0.18;
            else
                Par.hZ=0.2;
                Par.hP=0.12;%
            end
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(indExp) '_' num2str(ind) '_' num2str(OptionS) '_' num2str(OptionCov)],'Res')
    case 3
        Par.theta_sampler = 'GibbsHMC';
        Par.nsteps = 20;
        if not(OptionCov)
            if DataSet == 3
                Par.hZ=0.6;
                Par.hP=0.025;
            else
                Par.hZ=0.18;
                Par.hP=0.025;%
            end
        else
            if DataSet == 3
                Par.hZ=0.6;
                Par.hP=0.18;
            else
                Par.hZ=0.2;
                Par.hP=0.12;%
            end
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(indExp) '_' num2str(ind) '_' num2str(OptionS) '_' num2str(OptionCov)],'Res')
    case 4
        Par.theta_sampler = 'JointHMC';
        Par.nsteps = 1;
        if not(OptionCov)
            if DataSet == 3
                Par.h = 0.025;
            else
                Par.h = 0.025;
            end
        else
            if DataSet == 3
                Par.h=0.18;
            else
                Par.h=0.09;%
            end
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(indExp) '_' num2str(ind) '_' num2str(OptionS) '_' num2str(OptionCov)],'Res')
    case 5
        Par.theta_sampler = 'JointHMC';
        Par.nsteps = 10;
        if not(OptionCov)
            if DataSet == 3
                Par.h = 0.025;
            else
                Par.h = 0.025;
            end
        else
            if DataSet == 3
                Par.h=0.18;
            else
                Par.h=0.09;%
            end
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(indExp) '_' num2str(ind) '_' num2str(OptionS) '_' num2str(OptionCov)],'Res')
    case 6
        Par.theta_sampler = 'JointHMC';
        Par.nsteps = 20;
        if not(OptionCov)
            if DataSet == 3
                Par.h = 0.025;
            else
                Par.h = 0.025;
            end
        else
            if DataSet == 3
                Par.h=0.18;
            else
                Par.h=0.09;%
            end
        end
        Par.loop = loop;
        Res = RunJointMCMC_Full(Data,Par);
        save([SavePath '/D' num2str(DataSet) '_Exp' num2str(indExp) '_' num2str(ind) '_' num2str(OptionS) '_' num2str(OptionCov)],'Res')
    case 7
        if DataSet == 3
%             Par = Data.ParTrue;
            Par.Epsil = 1;
            Par.MCMCType = 'Rand';
            % Par.G = eye(length(Par.Names.Estimated));
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
            fullvolModel.LikFunction = 'normpdf(Data.Y(IndTime)-Data.Y(IndTime-1),Variables(:,2),sqrt(Variables(:,3)))';
            TempPar = ProposeInitialParameter(Data, fullvolModel, Par);
            Res = RunEstimationMethod(Data, fullvolModel, Par, TempPar, loop);
            save([SavePath '/D' num2str(DataSet) '_Exp' num2str(indExp) '_' num2str(ind) '_' num2str(OptionS) '_' num2str(OptionCov)],'Res')
        end
end
    

