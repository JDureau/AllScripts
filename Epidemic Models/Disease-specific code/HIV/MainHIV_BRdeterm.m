function [] =MainHIV_BRdeterm(ind)

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
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';



Data = ResGens{ind}.Data;
Parameters = ResGens{ind}.Parameters;
Parameters.DiffusionType = 'Bertallanfy';

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


Parameters.BRmm1.Estimated = 1;



Parameters.BRmu.Estimated = 1;
Parameters.BRbase.Estimated = 1;
Parameters.BRtinfl.Estimated = 1;


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

NbIts = 10000;
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
Parameters.SaveSpace = 0;
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
Parameters.SaveSpace = 1;
Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,150000);

if Constr
    save([SavePath '/ResGenBR' ModelType '_BRdetermConstr_' num2str(ind) '.mat'],'Res')
else
    save([SavePath '/ResGenBR' ModelType '_BRdeterm_' num2str(ind) '.mat'],'Res')
end

