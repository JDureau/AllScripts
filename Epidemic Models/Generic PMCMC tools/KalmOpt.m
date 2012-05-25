function Parameters = KalmOpt(Parameters,Data,Model,n)

if nargin == 3
    n = 15000;
end



Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
% Parameters.Correction = 1;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,Model,Parameters),Initialization,optimset('MaxFunEvals',Parameters.NbParsEstimated*n,'MaxIter',n,'TolX',1e-7,'TolFun',1e-7));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
Parameters.Finalx = x;
TellParsValues(Parameters)  

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).SamplStd = 0.00001;
end
Test = 0;
try
    ResKal = KalmanNumericDerivativesWithPrior(Data,Model,Parameters);
    Cov = (-ResKal.Hess)^-1;
    Test = IsDefPos(Cov);
end



disp(['Maximum reached:' num2str(Test)])

Parameters.KalmMaxReached = Test;
try
    Parameters.KalHess = ResKal.Hess;
end