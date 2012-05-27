function [] = HIVRealData_Determ(ind,IndModel)

ind = ind+1;

cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])


SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';

HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIV_Initialize;
HIVModel.SMC_projection = @HIV_SMC_projection;

NbIts = 20000;

if IndModel == 1
    s = '_dBR';
elseif IndModel == 2
    s = '_Sigm';
end


if or(ind == 1,ind == 2)
    temp = load([SavePath 'ParametersMysore.mat']);
    Parameters = temp.Parameters;
    
    if ind == 1
        ObsYears = [2004.667  2006.917    2008.834    2009.26];
        ObsMin = [0.2193 0.1911 0.032663 0.06975];
        ObsMax = [0.3028 0.2945 0.07557 0.14818];
        ObsVars = [7 7 8 7];
        Name = ['/HIV_Mysore_3rounds' s '.mat'];
    else
        ObsYears = [2004.667	2006.917    2008.834	];
        ObsMin = [0.2193 0.1911 0.032663 ];
        ObsMax = [0.3028 0.2945 0.07557 ];
        ObsVars = [7 7 8];
        Name = ['/HIV_Mysore_2rounds' s '.mat'];
    end
    
elseif or(ind == 3,ind == 4)
    temp = load([SavePath 'ParametersBelgaum.mat']);
    Parameters = temp.Parameters;
    
    if ind == 3
        ObsYears = [2005.834	2007.824    2008.584	2010.71];
        ObsMin = [0.2762 0.0363 0.2217  0.17552];
        ObsMax = [0.4018 0.0877 0.3251  0.2695];
        ObsVars = [7 8 7 7];
        Name = ['/HIV_Belgaum_3rounds' s '.mat'];
    else
        ObsYears = [2005.834	2007.824    2008.584 ];
        ObsMin = [0.2762  0.0363 0.2217];
        ObsMax = [0.4018  0.0877 0.3251];
        ObsVars = [7 8 7];
        Name = ['/HIV_Belgaum_2rounds' s '.mat'];
    end
    
elseif ind == 5
    temp = load([SavePath 'ParametersBellary.mat']);  
    Parameters = temp.Parameters;

    ObsYears = [2005.908	2007.873  2008.642	2010.865];
    ObsMin = [0.1106 0.0258 0.1048 0.0377];
    ObsMax = [0.2003 0.0946 0.1776 0.0892];
    ObsVars = [7 8 7 7];
    Name = ['/HIV_Bellary' s '.mat'];
    
elseif ind == 6
    temp = load([SavePath 'ParametersEastGodavry.mat']);
    Parameters = temp.Parameters;
    
    ObsYears = [2006.25  2006.84 2009.25 2009.35];
    ObsMin = [0.2004 0.0483 0.1514 0.041215];
    ObsMax = [0.3247 0.1184 0.3139 0.150556];
    ObsVars = [7 8 7 8];
    Name = ['/HIV_EastGodavry' s '.mat'];

elseif ind == 7
    temp = load([SavePath 'ParametersGuntur.mat']);  
    Parameters = temp.Parameters;
    
    ObsYears = [2006.38	2006.905   2009.53	2009.59];
    ObsMin = [0.1639 0.0369 0.0429 0.021891];
    ObsMax = [0.262 0.0956 0.125 0.1206];
    ObsVars = [7 8 7 8];

    Name = ['/HIV_Guntur' s '.mat'];


elseif ind == 8
    temp = load([SavePath 'ParametersHyderabad.mat']); 
    Parameters = temp.Parameters;
    
    ObsYears = [2006.16	2006.96   2009.47	2009.52];
    ObsMin = [0.0906 0.007 0.045 0];
    ObsMax = [0.1954 0.0405 0.147 0.0841];
    ObsVars = [7 8 7 8];
    Name = ['/HIV_Hyderabad' s '.mat'];
    

elseif ind == 9
    temp = load([SavePath 'ParametersYevatmal.mat']);  
    Parameters = temp.Parameters;
    
    ObsYears = [2006.337	2006.902    2009.728	2009.9];
    ObsMin = [0.2391 0.075527 0.1805 0.065032];
    ObsMax = [0.506 0.142355 0.3546 0.168657];
    ObsVars = [7 8 7 8];
    Name = ['/HIV_Yevatmal' s '.mat'];

elseif ind == 10
    temp = load([SavePath 'ParametersShimoga.mat']); 
    
    Parameters = temp.Parameters;
    ObsYears = [2005.688	2007.943  2008.718];
    Obs = [0.0968	 0.023 0.0896];
    ObsMin = [0.0632 0.0085 0.0566];
    ObsMax = [0.1305 0.0514 0.1226];
    ObsVars = [7 8 7];
    Name = ['/HIV_Shimoga' s '.mat'];

end






Parameters.NbVariables = 9;
Parameters.Problem = 'ImperialHIV';
Parameters.ObservationLength = 26*12;
Parameters.ComputationTStep = 0.5;
Parameters.TypeWork = 'Normal';



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
    HIVModel.ObservationMeasurementNoise{i+1} = 0.05*(ObsMax(i)+ObsMin(i))*100/2;%((ObsMax(i)-ObsMin(i))*100/4)^2;%(Data.Observations(ObsVars(i),i+1)*(100-Data.Observations(ObsVars(i),i+1))/400);
end




Parameters.InitialFt.Estimated = 0;
Parameters.SigmaRW.Estimated = 0;
Parameters.BRsigma.Estimated = 0;
Parameters.Sigmsigma.Estimated = 0;
if IndModel == 1
    Parameters.Sigmrate.Estimated = 0;
    Parameters.Sigmmu.Estimated = 0;
    Parameters.Sigmbase.Estimated = 0;
    Parameters.Sigmtinfl.Estimated = 0;
    Parameters.BRmm1.Estimated = 1;
    Parameters.BRmu.Estimated = 1;
    Parameters.BRbase.Estimated = 1;
    Parameters.BRtinfl.Estimated = 1;
elseif IndModel == 2
    Parameters.BRmm1.Estimated = 0;
    Parameters.BRmu.Estimated = 0;
    Parameters.BRbase.Estimated = 0;
    Parameters.BRtinfl.Estimated = 0;
    Parameters.Sigmrate.Estimated = 1;
    Parameters.Sigmmu.Estimated = 1;
    Parameters.Sigmbase.Estimated = 1;
    Parameters.Sigmtinfl.Estimated = 1;
end
Parameters.SigmaObs = 0.05;

Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

if IndModel == 1
    Parameters.DiffusionType = 'Bertallanfy';
elseif IndModel == 2
    Parameters.DiffusionType = 'Sigmoid';
end    
    
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
if IndModel == 1
    [x,fval,exitflag,output] = fminsearch(@(x) SquarredErrorBRtoMinimize(x,NamesEst,Data,Parameters,HIVModel),Initialization,optimset('MaxIter',5000,'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',10000));
elseif IndModel == 2
    [x,fval,exitflag,output] = fminsearch(@(x) SquarredErrorSigmtoMinimize(x,NamesEst,Data,Parameters,HIVModel),Initialization,optimset('MaxIter',5000,'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',10000));
end
for i = 1:length(NamesEst)
    Parameters.(NamesEst{i}).TransfValue = (x(i));
end
finalx = x;
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)
InitParameters = Parameters;

if IndModel == 1
    fun =  @(x)SquarredErrorBRtoMinimize(x,NamesEst,Data,Parameters,HIVModel);
elseif IndModel == 2
    fun =  @(x)SquarredErrorSigmtoMinimize(x,NamesEst,Data,Parameters,HIVModel);
end
[HD,err] = hessdiag(fun,finalx);       
Hess = diag(-HD);
try
    test =  (mean(eig(-Hess)>0)==1)*isreal(eig(-Hess));
end
Cov = (-Hess)^-1;
    
    

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
    if IndModel == 1
        [x,fval,exitflag,output] = fminsearch(@(x) SquarredErrorBRtoMinimize(x,NamesEst,Data,Parameters,HIVModel),Initialization,optimset('MaxIter',5000,'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',10000));
    elseif IndModel == 2
        [x,fval,exitflag,output] = fminsearch(@(x) SquarredErrorSigmtoMinimize(x,NamesEst,Data,Parameters,HIVModel),Initialization,optimset('MaxIter',5000,'TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',10000));
    end
    for i = 1:length(NamesEst)
        Parameters.(NamesEst{i}).TransfValue = (x(i));
    end
    finalx = x;
    Parameters = UpdateParsTransfToNoTransf(Parameters);
    TellParsValues(Parameters)
    InitParameters = Parameters;
    disp('Hessian')
    if IndModel == 1
        fun =  @(x)SquarredErrorBRtoMinimize(x,NamesEst,Data,Parameters,HIVModel);
    elseif IndModel == 2
        fun =  @(x)SquarredErrorSigmtoMinimize(x,NamesEst,Data,Parameters,HIVModel);
    end
    % [HD,err] = hessian(fun,Parameters.Finalx);
    % Hess = -HD;
    [HD,err] = hessdiag(fun,finalx);       
%     [Hess,err] = hessian(fun,finalx); 
    Hess = diag(-HD);
%     disp(eig(Hess))
    try
        test =  (mean(eig(-Hess)>0)==1)*isreal(eig(-Hess));
    end
    cpt = cpt+1;
    if cpt>5
        test = 1;
    end
end
Cov = (-Hess)^-1;

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


Parameters.mu = finalx;
Parameters.NamesEst = NamesEst;
Parameters.SaveSpace = 0;
Parameters.AdMet = 1;
Parameters.AdMetBeta = 0.05;
Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,NbIts,IndModel);
Res.Parameters.TypeWork='Normal';
PlotResHIV(Res,Res.Parameters)



for i = 1:5
    
    disp(['Burn-in ' num2str(i)])
    Parameters = Res.Parameters;
    if not(sum(eig(cov(Res.TransfThetas'))<=0))
        Cov = 2.38^2/length(finalx)*cov(Res.TransfThetas');
    end
    Parameters.NamesEst = NamesEst;
    Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,NbIts,IndModel);
    Parameters = Res.Parameters;
end

disp('MCMC')
Parameters = Res.Parameters;
Cov = 2.38^2/length(finalx)*cov(Res.TransfThetas');
Parameters.NamesEst = NamesEst;
Parameters.SaveSpace = 1;
clear Res;
Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,200000,IndModel);


save([SavePath Name],'Res')





