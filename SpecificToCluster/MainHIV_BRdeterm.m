function [] =MainHIV_BRdeterm(ind,IndModel)

 ModelType = 'Bert';
%ModelType = 'Bert10';
%ModelType = 'BertmInf';
% ModelType = 'Step';
% ModelType = 'Affine';
% ModelType = 'Logist';
%ModelType = 'Sigm';


Constr = 0;

ind = ind+1;

%% load paths                                                                                                                                                                                              


cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])


%% Go                                                                                                                                                                                                      


SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';


if strcmp(ModelType,'Bert')
    load([SavePath '/ResGenBR.mat'])    
elseif strcmp(ModelType,'Bert10')
    load([SavePath '/ResGenBR10.mat'])
elseif strcmp(ModelType,'BertmInf')
    load([SavePath '/ResGenBRmInf.mat'])
elseif strcmp(ModelType,'Logist')
    load([SavePath '/ResGenLogist.mat'])
elseif strcmp(ModelType,'Sigm')
    load([SavePath '/ResGenSigm.mat'])
elseif strcmp(ModelType,'Step')
    load([SavePath '/ResGenStep.mat'])
elseif strcmp(ModelType,'Step10')
    load([SavePath '/ResGenStep10.mat'])
elseif strcmp(ModelType,'Affine')
    load([SavePath '/ResGenAffine.mat'])
elseif strcmp(ModelType,'Affine10')
    load([SavePath '/ResGenAffine10.mat'])
end

HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIV_Initialize;
HIVModel.SMC_projection = @HIV_SMC_projection;
% % HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';



Data = ResGens{ind}.Data;
Parameters = ResGens{ind}.Parameters;
if IndModel == 1
    Parameters.DiffusionType = 'Bertallanfy';
elseif IndModel == 2
    Parameters.DiffusionType = 'Sigmoid';
end
    
if Constr
    Names = Parameters.Names.Estimated;
    for i = 1:length(Names)
        Parameters.(Names{i}).Estimated = 0;
    end
    Parameters = DefineEstimatedParametersIndexes(Parameters);
    Parameters = DefineTransfFunctions(Parameters);
    Parameters = DefinePriors(Parameters);
    Parameters = UpdateParsNoTransfToTransf(Parameters);
end


Parameters.InitialFt.Estimated = 0;
Parameters.SigmaRW.Estimated = 0;
Parameters.BRsigma.Estimated = 0;

if IndModel == 1
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
    Parameters.Sigmrate.Value = 6;
    Parameters.Sigmrate.Min = -10^14;
    Parameters.Sigmrate.Max = 10^14;
    Parameters.Sigmrate.MaxLim = 1000;
    Parameters.Sigmrate.MinLim = 0;
    % Parameters.Sigmmm1.MaxLim = 10;
    Parameters.Sigmrate.Estimated = 1;
    Parameters.Sigmrate.TransfType = 'Logit';
    Parameters.Sigmrate.Init = 1;
    Parameters.Sigmmu.Value = 0.8;
    Parameters.Sigmmu.Min = -10^14;
    Parameters.Sigmmu.Max = 10^14;
    Parameters.Sigmmu.MaxLim = 0.99;
    Parameters.Sigmmu.MinLim = 0.01;
    Parameters.Sigmmu.Estimated = 1;
    Parameters.Sigmmu.TransfType = 'Logit';
    Parameters.Sigmmu.Init = 1;
    Parameters.Sigmbase.Value = 0.05;
    Parameters.Sigmbase.Min = -10^14;
    Parameters.Sigmbase.Max = 10^14;
    Parameters.Sigmbase.MaxLim = 0.9;
    % Parameters.Sigmbase.MaxLim = 0.1;
    Parameters.Sigmbase.MinLim = 0;
    Parameters.Sigmbase.Estimated = 1;
    Parameters.Sigmbase.TransfType = 'Logit';
    Parameters.Sigmbase.Init = 1;
    Parameters.Sigmtinfl.Value = 100;
    Parameters.Sigmtinfl.Min = -10^14;
    Parameters.Sigmtinfl.Max = 10^14;
    Parameters.Sigmtinfl.MinLim = 1;
    Parameters.Sigmtinfl.MaxLim = 300;
    Parameters.Sigmtinfl.Estimated = 1;
    Parameters.Sigmtinfl.TransfType = 'Logit';
    Parameters.Sigmtinfl.Init = 1;
end

Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);



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
    [x,fval,exitflag,output] = fminsearch(@(x) SquarredErrorBRtoMinimize(x,NamesEst,Data,Parameters,HIVModel),Initialization,optimset('MaxIter',2000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',4000));
elseif IndModel == 2
    [x,fval,exitflag,output] = fminsearch(@(x) SquarredErrorSigmtoMinimize(x,NamesEst,Data,Parameters,HIVModel),Initialization,optimset('MaxIter',2000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',4000));
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
% [Hess,err] = hessian(fun,finalx); 

disp(HD)
Hess = diag(-HD);
% disp(eig(Hess))

NbIts = 20000;
Parameters = InitParameters;
Cov = (-Hess)^-1;
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
    end      % [HD,err] = hessian(fun,Parameters.Finalx);
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


for i = 1:5
    
    disp(['Burn-in ' num2str(i)])
    Parameters = Res.Parameters;
    Cov = cov(Res.TransfThetas');
    Parameters.NamesEst = NamesEst;
    Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,NbIts,IndModel);
    Parameters = Res.Parameters;
end

disp('MCMC')
Parameters = Res.Parameters;
Cov = 2.38^2/length(finalx)*cov(Res.TransfThetas');
Parameters.NamesEst = NamesEst;
Parameters.SaveSpace = 1;
Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,130000,IndModel);

if IndModel == 1
    save([SavePath '/ResGenBR' ModelType '_BRdeterm_' num2str(ind) '.mat'],'Res')
elseif IndModel == 2
    save([SavePath '/ResGenBR' ModelType '_Sigmdeterm_' num2str(ind) '.mat'],'Res')
end

