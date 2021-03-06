 function Result = EstimationEKFGen(Data, Model, Parameters)

Parameters = Model.InitializeParameters(Parameters);
 
CheckParametersGen(Parameters);

if not(isfield(Parameters,'RWinEKF'))
    Parameters.RWinEKF = 0;
end

if not(isfield(Parameters,'RunningMif'))
    Parameters.RunningMif = 0;
end


tic
NbVarsTot = Parameters.NbVariables;
Observations = Data.Observations;
ObservationInstants = Data.Instants;
NbComputingSteps = Data.NbComputingSteps;
Parameters.NbObs = length(ObservationInstants) ;

TStep     = Parameters.ComputationTStep;

InitialState = Parameters.InitialState;

if Parameters.RunningMif 
    names = Parameters.Names.Estimated;
    NbVarsTot = NbVarsTot + length(names);
    for i = 1:length(names)
        InitialState(end+1) = Parameters.(names{i}).TransfValue;
    end
end
    


m = InitialState;
Cov = zeros(NbVarsTot,NbVarsTot);
try  
    Ps = 0.00*diag(m);
%     Ps = Parameters.InitialCov;
catch
    Ps = Cov;
end

Cov = Ps;

Liks = 1;
LogLik = 0;
LiksMeanTraj = 1;
LogLikMeanTraj = 0;
LiksTheoret = 1;
LogLikTheoret = 0;
ms = repmat(m,1,length(ObservationInstants));
Covs = [];
Covs(1,:,:)=Cov;
endit = 0;
% disp(m(3))

for IndTime = 2:length(ObservationInstants) 
%     IndTime
    if not(mean(Cov == Cov')==1)
        disp('pb')
    end
    Res = Model.EKF_projection(Data,Model,m,Cov,NbComputingSteps(IndTime),IndTime,Parameters);
    Model = Res.Model;
    mpred = Res.m;
    
%     mpred
    
%     diag(Res.Cov)
%     die
    
%     mpred(1)
%     mpred(2)
%     mpred(3)
%     mpred(4)
%     mpred(5)
%     mpred(6)
%     mpred(7)
%     mpred(8)
%     mpred(9)
    
    Cov1 = Res.Cov;
%     'before'
%     for i = 1:4
%         disp(mpred(i));
%     end
%     disp(Cov1)
%     for i = 1:9
% %         disp(Cov1(i,i))
%         disp(mpred(i))
%     end
%     die
    
    if not(mean(Cov1 == Cov1')==1)
        disp('pb')
    end
    
    if Res.Crash
        endit = 1;
    end

    try if endit
            LogLik = -Inf;
        else
            die
        end    
    catch

% 
%         Res.m
% %         Res.Cov
%         die


        Ck = Model.ObservationJacobian{IndTime};

        % new attempt
        predCov = Cov1;

        try
            coeff = Parameters.MultCoeff.Value;
        catch
            coeff = 1;
        end
        
        inds = Data.ObservedVariables{IndTime};
        if not(isempty(inds))
            IndObservedVar = Data.ObservedVariables{IndTime};
            if and(strcmp(Parameters.Problem,'ImperialHIV2'),length(IndObservedVar) == 2)
                ypred = [];
                ypred(1,1) = mpred(IndObservedVar(1));
                ypred(2,1) = Parameters.Rho.Value*(InvLogitTransf(mpred(IndObservedVar(2)),0,1))*100;
%                 disp([ Observations(IndObservedVar,IndTime)  ypred])
            else
                ypred = mpred(IndObservedVar);  
            end
            vk = Observations(IndObservedVar,IndTime) - coeff*ypred; 

           
%             'x'
%             ypred
%             
%             'modx'
%             coeff*ypred
%             
%             'obs'
%             Observations(IndObservedVar,IndTime)
            
%             'pred_eror:'
%             vk
%             
         

            
            Rk = diag(Model.ObservationMeasurementNoise{IndTime}');
           
            try
                Sk = Ck*Cov1*Ck' + Rk;
                
            catch
                'yo'
            end
            
%             'rk'
%             Rk
%             
%             'sk'
%             Sk
            
        %         
        %     Sk
        %     die
            if Sk<0
                disp('stop')
            end
            Kk = Cov1*Ck'*(Sk^-1);

%             'Ck'
%             Ck
%             
%             'Kk'
%             Kk
            
            m = mpred + (Kk*vk);
            
%             [mpred m]
%             die
        %      Kk
        %     m(4)
        %     ypred
        %     Sk
        %     die
            temp = (Kk*Sk*Kk' + (Kk*Sk*Kk')')/2; 
            Cov = Cov1 - temp;
        %     if not(mean(Cov1 == Cov1')==1)
        %         IndTime
        %         disp('pb')
        %     end
        %     if sum(not(isreal(Cov)))
        %         IndTime
        %         disp('pb')
        %     end
            Cov = (Cov + Cov')/2;
            try
                [V,D] = eig(Cov);
                temp = real(eig(Cov));
                temp = max(temp,0);
                D = diag(temp);
                temp = ( V*D*V' + (V*D*V')')/2;
                Cov = temp;
            catch
                'pb'
            end
            
      
            
%             Cov1 = Res.Cov;
%             'after'
%             for i = 1:4
%                 disp(m(i));
%             end
%             disp(Cov)
            

             if sum(not(isreal(Cov)))
                 disp('pb')
             end
            if sum((isnan(Cov)))
                 disp('pb')
             end


            ms(:,IndTime) = m';% transposed this for simf, may have to take it off
            Covs(IndTime,:,:) = Cov;

            if isnan(Sk)
                disp('oups')
            end

            try
        %         if strcmp(Parameters.DiffusionType,'AffineAdd')
        %             if Parameters.RWinEKF
        %                 t = sum(log(normpdf(Res.deltabetas,0,sqrt(Parameters.ComputationTStep)*Parameters.SigmaRW.Value)));
        %             else 
        %                 t = 0;
        %             end
        %             tempLogLik = t + log(mvnpdf(vk,zeros(size(vk)),Sk));
        %         elseif strcmp(Parameters.DiffusionType,'AffineInt')
        %             t = 0;%sum(log(normpdf(diff(Res.deltabetas)/Parameters.ComputationTStep,0,sqrt(Parameters.ComputationTStep)*max(eps,Parameters.SigmaRW.Value))));
        %             tempLogLik = t + log(mvnpdf(vk,zeros(size(vk)),Sk));
        %         else
        Sk = max(0.1,Sk); %Just commented it june 6th
                    tempLogLik = max(-700,log(not(Res.Crash)*mvnpdf(vk,zeros(size(vk)),Sk)));
        %         end
            catch
                disp('pb')
                tempLogLik = -700;
            end
            
%             'loglik'
%             tempLogLik


            Liks(IndTime) = exp(tempLogLik);
            LogLik = LogLik + tempLogLik;
            if isinf(LogLik)
                disp('inf loglik')
            else if not(isreal(LogLik))
                   'pb'
                end
            end
        end
    end
%      'after'
%     disp(m)
%     disp(Cov)
%    die
end

LogPrior = 0;
LogCorr = 0;
Names = Parameters.Names.Estimated;

for i = 1:length(Names)
    temp = Parameters.(Names{i}).Prior(Names{i},Parameters);
    LogPrior = LogPrior +log(temp);
    temp = Parameters.(Names{i}).CorrFunct(Names{i},Parameters);
    LogCorr = LogCorr + log(temp);
%     if strcmp(Names{i},'BRmm1')
%         disp(['m' ' ' num2str(Parameters.(Names{i}).Value) ' ' num2str(log(temp))])
%     end
%     if strcmp(Names{i},'SigmaRW')
%         disp(['sigma' ' ' num2str(Parameters.(Names{i}).Value) ' ' num2str(log(temp))])
%     end
end
Result.LogPrior = LogPrior;
Result.LogCorr = LogCorr;

Result.Data = Data;
Result.Parameters = Parameters;
Result.PosteriorMeans = ms;
Result.PosteriorCovs = Covs;
Sigmas = [];
for i = 1:NbVarsTot
    Sigmas(i,:) = sqrt(Covs(:,i,i));
end
try
    Result.Posterior975 = ms + 2*Sigmas;
    Result.Posterior025 = ms - 2*Sigmas;
end
Result.Likelihood = prod(Liks);
Result.Liks = Liks;
if not(isreal(LogLik))
    LogLik = 0;
end
% LogLik = LogLik + log(normpdf(exp(ms(6,4))-exp(ms(6,17)),0,0.1));
Result.LogLik = LogLik;
% disp(Parameters.SigmaRW.Value)

