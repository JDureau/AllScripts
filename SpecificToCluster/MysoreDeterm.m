function [] = MysoreDeterm()


cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])


SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';

temp = load([SavePath 'ParametersMysore10.mat']);
Parameters = temp.Parameters;

ObsYearsMysore3rds = [2004.667  2006.917    2008.834    2009.26];
CIMysore3rds = [];
ObsMysore3rds = [0.2611 0.2424 0.054 0.111];
ObsMysore3rdsMin = [0.2193 0.1911 0.032663 0.06975];
ObsMysore3rdsMax = [0.3028 0.2945 0.07557 0.14818];
ObsVarsMysore3rds = [7 7 8 7];



HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIV_Initialize;
HIVModel.SMC_projection = @HIV_SMC_projection;
Parameters.NbVariables = 9;
Parameters.Problem = 'ImperialHIV';
Parameters.ObservationLength = 25*12;
Parameters.ComputationTStep = 0.5;
Parameters.TypeWork = 'Normal';

ObsVars = ObsVarsMysore3rds;
ObsMin = ObsMysore3rdsMin;
ObsMax = ObsMysore3rdsMax;
ObsYears = ObsYearsMysore3rds;

Data.Observations = zeros(10,length(ObsVars)+1);
for i = 1:length(ObsVars)
    Data.Observations(ObsVars(i),1+i) = (ObsMin(i)+ObsMax(i))/2*100;
    Data.ObsSigmas(i+1) = ((ObsMax(i)-ObsMin(i))*100/4);
end
Instants = round((ObsYears-1985)*12);
Data.Instants = round([0 Instants]/(Parameters.ComputationTStep));
Data.ObservedVariables = [ 0 ObsVars];
Data.NbComputingSteps = [0 diff(Data.Instants)];
Parameters.NbTSteps = sum(Data.NbComputingSteps);

temp7 = zeros(1,9);
temp8 = zeros(1,9);
temp7(1,7) = 1;
temp8(1,8) = 1;
temps{7} = temp7;
temps{8} = temp8;
HIVModel.ObservationJacobian = {};
for i = 1:length(ObsVars)
    HIVModel.ObservationJacobian{i+1} = temps{ObsVars(i)};
end
HIVModel.ObservationMeasurementNoise = {};
for i = 1:length(ObsVars)
    HIVModel.ObservationMeasurementNoise{i+1} = ((ObsMax(i)-ObsMin(i))*100/4)^2;%(Data.Observations(ObsVars(i),i+1)*(100-Data.Observations(ObsVars(i),i+1))/400);
end




Parameters.InitialFt.Estimated = 0;
Parameters.SigmaRW.Estimated = 0;
Parameters.BRmm1.Estimated = 1;
Parameters.BRmu.Estimated = 1;
Parameters.BRbase.Estimated = 1;
Parameters.BRtinfl.Estimated = 1;
Parameters.SigmaObs = 0.1;

Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


Parameters.DiffusionType = 'Bertallanfy';

disp('Optim')
Names = Parameters.Names.Estimated;
NamesEst = {};
for i = 1:length(Names)
    if not(or(strcmp(Names{i},'SigmaRW'),strcmp(Names{i},'InitialFt')))
        NamesEst{end+1} = Names{i};
    end
end
Initialization = [];
for i = 1:length(NamesEst)
    Initialization(i) = Parameters.(NamesEst{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) SquarredErrorBRtoMinimize(x,NamesEst,Data,Parameters,HIVModel),Initialization,optimset('MaxIter',2000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',4000));
for i = 1:length(NamesEst)
    Parameters.(NamesEst{i}).TransfValue = (x(i));
end
finalx = x;
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)
InitParameters = Parameters;


disp('Hessian')
fun =  @(x)SquarredErrorBRtoMinimize(x,NamesEst,Data,Parameters,HIVModel);
% [HD,err] = hessian(fun,Parameters.Finalx);
% Hess = -HD;
[HD,err] = hessdiag(fun,finalx); 
% [Hess,err] = hessian(fun,finalx); 

disp(HD)
Hess = diag(-HD);
% disp(eig(Hess))

NbIts = 5000;
Parameters = InitParameters;
Cov = 0.01*2.38^2/length(finalx)*(-Hess)^-1;
test =  (mean(eig(-Hess)>0)==1)*isreal(eig(-Hess));

% test = 1;
cpt = 1;
while not(test)
    disp('Optim')
    Names = Parameters.Names.Estimated;
    NamesEst = {};
    for i = 1:length(Names)
        if not(or(strcmp(Names{i},'SigmaRW'),strcmp(Names{i},'InitialFt')))
            NamesEst{end+1} = Names{i};
        end
    end
    Initialization = [];
    for i = 1:length(NamesEst)
        Initialization(i) = Parameters.(NamesEst{i}).TransfValue ;
    end
    [x,fval,exitflag,output] = fminsearch(@(x) SquarredErrorBRtoMinimize(x,NamesEst,Data,Parameters,HIVModel),Initialization,optimset('MaxIter',5000,'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',10000));
    for i = 1:length(NamesEst)
        Parameters.(NamesEst{i}).TransfValue = (x(i));
    end
    finalx = x;
    Parameters = UpdateParsTransfToNoTransf(Parameters);
    TellParsValues(Parameters)
    InitParameters = Parameters;
    disp('Hessian')
    fun =  @(x)SquarredErrorBRtoMinimize(x,NamesEst,Data,Parameters,HIVModel);
    % [HD,err] = hessian(fun,Parameters.Finalx);
    % Hess = -HD;
    [HD,err] = hessdiag(fun,finalx);       
%     [Hess,err] = hessian(fun,finalx); 
    Hess = diag(-HD);
    disp(eig(Hess))
    test =  (mean(eig(-Hess)>0)==1)*isreal(eig(-Hess));
    cpt = cpt+1;
    if cpt>5
        test = 1;
    end
end

% 
% fvals = fval;
% Pars{1} = Parameters;
% for k = 2:3
%     rdind = ceil(rand(1,1)*length(ResGens));
%     TempPars = ResGens{rdind}.Parameters;
%     TempPars.DiffusionType = 'Bertallanfy';
%     TempPars.BRm.Min = -10^14;
%     TempPars.BRm.Max = 10^14;
%     TempPars.BRm.MaxLim = 5000;
%     Parameters.BRm.MinLim = 1.01;
%     if strcmp(Type,'Logistic')
%         Parameters.BRm.Value = 2;
%         Parameters.BRm.Estimated = 0;
%     else
%         Parameters.BRm.Estimated = 1;
%     end
%     TempPars.BRm.TransfType = 'Logit';
%     TempPars.BRm.Init = 1;
%     TempPars.BRmu.Min = -10^14;
%     TempPars.BRmu.Max = 10^14;
%     TempPars.BRmu.MaxLim = 0.99;
%     TempPars.BRmu.MinLim = 0.01;
%     TempPars.BRmu.Estimated = 1;
%     TempPars.BRmu.TransfType = 'Logit';
%     TempPars.BRmu.Init = 1;
%     TempPars.BRbase.Min = -10^14;
%     TempPars.BRbase.Max = 10^14;
%     TempPars.BRbase.MaxLim = 0.9;
%     TempPars.BRbase.MinLim = 0.01;
%     TempPars.BRbase.Estimated = 1;
%     TempPars.BRbase.TransfType = 'Logit';
%     TempPars.BRbase.Init = 1;
%     TempPars.BRtinfl.Min = -10^14;
%     TempPars.BRtinfl.Max = 10^14;
%     TempPars.BRtinfl.MinLim = 1;
%     TempPars.BRtinfl.MaxLim = 600;
%     TempPars.BRtinfl.Estimated = 1;
%     TempPars.BRtinfl.TransfType = 'Logit';
%     TempPars.BRtinfl.Init = 1;
%     TempPars = DefineEstimatedParametersIndexes(TempPars);
%     TempPars = DefineTransfFunctions(TempPars);
%     TempPars = DefinePriors(TempPars);
%     TempPars = UpdateParsNoTransfToTransf(TempPars);
%     Initialization = [];
%     for i = 1:length(NamesEst)
%         Initialization(i) = TempPars.(NamesEst{i}).TransfValue ;
%     end
%     [x,fval,exitflag,output] = fminsearch(@(x) SquarredErrorBRtoMinimize(x,NamesEst,Data,TempPars,HIVModel),Initialization,optimset('MaxIter',3000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',6000));
%     for i = 1:length(NamesEst)
%         TempPars.(NamesEst{i}).TransfValue = (x(i));
%     end
%     finalx = x;
%     TempPars = UpdateParsTransfToNoTransf(TempPars);
%     TellParsValues(TempPars)
%     fvals(k) = fval;
%     Pars{k} = TempPars;
% end
% fvals
% [b,minind] = min(fvals)
% Parameters = Pars{minind};


InitParameters = Parameters;



Parameters.NamesEst = NamesEst;
Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,NbIts);

for i = 1:5
    
    disp(['Burn-in ' num2str(i)])
    Parameters = Res.Parameters;
    Cov = 2.38^2/length(finalx)*cov(Res.TransfThetas');
    Parameters.NamesEst = NamesEst;
    Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,NbIts);
    Parameters = Res.Parameters;
end

disp('MCMC')
Parameters = Res.Parameters;
Cov = 2.38^2/length(finalx)*cov(Res.TransfThetas');
Parameters.NamesEst = NamesEst;
Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,40000);

save([SavePath '/HIV_Mysore_3roundsDetBertallanfy10.mat'],'Res')






