function Res = KalmanNumericDerivativesWithPrior(Data,Model,Parameters)

CheckParametersGen(Parameters)
  
% Grad and Hess
xs = [];
hs = [];
Grad = [];
Hess = [];
dxdTransfx = [];


% fun = @(x,y) x.^2 + y.^2;
% xy = [2 3];
% gradvec = [derivest(@(x) fun(x,xy(2)),xy(1),'d',1), ...
%            derivest(@(y) fun(xy(1),y),xy(2),'d',1)]
% [HD,err,finaldelta] = hessdiag(fun,xy')
fun =  @(x)KalmanToOptimizeWithPrior(x,Data,Model,Parameters);
% [HD,err] = hessian(fun,Parameters.Finalx);
% Hess = -HD;

[HD,err] = hessdiag(fun,Parameters.Finalx);       
Hess = diag(-HD);

% Names = Parameters.Names.Estimated;
% for i = 1:length(Names)
%     xs(i) = Parameters.(Names{i}).TransfValue;
%     Temp = EstimationEKFGen(Data, Model,Parameters);
%     fxs(i) = Temp.LogLik + Temp.LogPrior ;%  - Temp.LogCorr ;
%     hs(i) = Parameters.(Names{i}).SamplStd ;
%     xphs(i) =  xs(i) + hs(i);
%     xmhs(i) =  xs(i) - hs(i);
%     TempParameters = Parameters;
%     TempParameters.(Names{i}).TransfValue = xphs(i);
% %     TempParameters = UpdateParsNoTransfToTransf(TempParameters);
%     TempParameters = UpdateParsTransfToNoTransf(TempParameters);
%     Temp = EstimationEKFGen(Data, Model,TempParameters);
%     fxphs(i) = Temp.LogLik+ Temp.LogPrior ;% - Temp.LogCorr;
%     TempParameters.(Names{i}).TransfValue = xmhs(i);
% %     TempParameters = UpdateParsNoTransfToTransf(TempParameters);
%     TempParameters = UpdateParsTransfToNoTransf(TempParameters);
%     Temp = EstimationEKFGen(Data,Model, TempParameters);
%     fxmhs(i) = Temp.LogLik+ Temp.LogPrior ;% - Temp.LogCorr;
%     if strcmp(Parameters.(Names{i}).TransfType , 'Log')
%         dxdTransfx(i) = xs(i);
%     elseif strcmp(Parameters.(Names{i}).TransfType , 'Logit')
%         Transfx = Parameters.(Names{i}).TransfValue;
%         try
%             m = Parameters.(Names{i}).MinLim;
%             M = Parameters.(Names{i}).MaxLim;
%         catch
%             m = Parameters.(Names{i}).Min;
%             M = Parameters.(Names{i}).Max;
%         end
%         dxdTransfx(i) = (exp(Transfx)*(M-m))/((1+exp(Transfx))^2);
%     end
% %     Grad(Parameters.(Names{i}).Index) =  dxdTransfx(i)*( fxphs(i) - fxmhs(i) )/ (2*hs(i));
% %     Hess(Parameters.(Names{i}).Index,Parameters.(Names{i}).Index) =  dxdTransfx(i)*dxdTransfx(i)*( fxphs(i) - 2*fxs(i) + fxmhs(i) )/ (hs(i)^2);
%     Grad(Parameters.(Names{i}).Index) =  ( fxphs(i) - fxmhs(i) )/ (2*hs(i));
%     Hess(Parameters.(Names{i}).Index,Parameters.(Names{i}).Index) =  ( fxphs(i) - 2*fxs(i) + fxmhs(i) )/ (hs(i)^2);
% 
%     disp(Names{i})
%     disp([fxs(i) fxphs(i)-fxs(i) fxmhs(i)-fxs(i)])
% 
% end
% 
% 
% 
% for i = 1:length(Names)-1
%     for j = i+1:length(Names)
%         % x+hx z+hy
%         TempParameters = Parameters;
%         TempParameters.(Names{i}).TransfValue = xphs(i);
%         TempParameters.(Names{j}).TransfValue = xphs(j);
%         %     TempParameters = UpdateParsNoTransfToTransf(TempParameters);
%         TempParameters = UpdateParsTransfToNoTransf(TempParameters);
%         Temp = EstimationEKFGen(Data,Model, TempParameters);
%         fxphzph = Temp.LogLik+ Temp.LogPrior ;% - Temp.LogCorr;
%         % x+h z-h
%         TempParameters = Parameters;
%         TempParameters.(Names{i}).TransfValue = xphs(i);
%         TempParameters.(Names{j}).TransfValue = xmhs(j);
%        %     TempParameters = UpdateParsNoTransfToTransf(TempParameters);
%         TempParameters = UpdateParsTransfToNoTransf(TempParameters);
%         Temp = EstimationEKFGen(Data,Model, TempParameters);
%         fxphzmh = Temp.LogLik+ Temp.LogPrior ;% - Temp.LogCorr;
%         % x-h z+h
%         TempParameters = Parameters;
%         TempParameters.(Names{i}).TransfValue = xmhs(i);
%         TempParameters.(Names{j}).TransfValue = xphs(j);
%   %     TempParameters = UpdateParsNoTransfToTransf(TempParameters);
%         TempParameters = UpdateParsTransfToNoTransf(TempParameters);
%         Temp = EstimationEKFGen(Data,Model, TempParameters);
%         fxmhzph = Temp.LogLik+ Temp.LogPrior ;% - Temp.LogCorr;
%         % x-h z-h
%         TempParameters = Parameters;
%         TempParameters.(Names{i}).TransfValue = xmhs(i);
%         TempParameters.(Names{j}).TransfValue = xmhs(j);
%  %     TempParameters = UpdateParsNoTransfToTransf(TempParameters);
%         TempParameters = UpdateParsTransfToNoTransf(TempParameters);
%         Temp = EstimationEKFGen(Data,Model, TempParameters);
%         fxmhzmh = Temp.LogLik+ Temp.LogPrior ;% - Temp.LogCorr;
% 
% %         Hess(Parameters.(Names{i}).Index,Parameters.(Names{j}).Index) = dxdTransfx(i)*dxdTransfx(j)*(fxphzph - fxphzmh - fxmhzph + fxmhzmh)/(4*hs(i)*hs(j));
% %         Hess(Parameters.(Names{j}).Index,Parameters.(Names{i}).Index) = dxdTransfx(i)*dxdTransfx(j)*(fxphzph - fxphzmh - fxmhzph + fxmhzmh)/(4*hs(i)*hs(j));
%         Hess(Parameters.(Names{i}).Index,Parameters.(Names{j}).Index) = (fxphzph - fxphzmh - fxmhzph + fxmhzmh)/(4*hs(i)*hs(j));
%         Hess(Parameters.(Names{j}).Index,Parameters.(Names{i}).Index) = (fxphzph - fxphzmh - fxmhzph + fxmhzmh)/(4*hs(i)*hs(j));
% 
%         disp(Names{i})
%         disp(Names{j})
%         disp([fxs(i) fxphzph-fxs(i) fxphzmh-fxs(i) fxmhzph-fxs(i) fxmhzmh-fxs(i)])
%     end
% end

    
Res.Grad = Grad;
Res.Hess = Hess;
Cov = (-Hess)^-1;
Test = IsDefPos(Cov);
disp(Test)    