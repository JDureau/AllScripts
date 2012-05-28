function [] = SEIRMarcTestsKal(ind,indeps)
ind = ind+1;
indeps = indeps + 1;

cd('/users/ecologie/dureau/src/AllScripts/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR'])


SavePath = '/users/ecologie/dureau/src/AllData/ResultsMarc/';

load([SavePath '/RealTime_Add_28_cluster5.mat']);

PostCov = cov(Res3.TransfThetas');
SEIRModel = Res3.Model;
SEIRModel.LikFunction = 'normpdf(log(Variables(:,5)),transpose(log(coeff*Data.Observations(5,IndTime))-log(Parameters.SigmaObs.Value^2+1)/2),sqrt(log(Parameters.SigmaObs.Value^2+1)))';%Parameters.SigmaObs.Value)';
Parameters = Res3.Parameters;
Data = Res3.Data;
TempPar = Res3.TempPar;
Parameters.NbParticules = 3000;
Parameters.Correction = 1;
Parameters.Problem = 'MarcFluPlusObs';
Parameters = KalmOpt(Parameters,Data,SEIRModel,2000);
Test = 0;
NbIts = 0;
while not(Test)
    Parameters = KalmOpt(Parameters,Data,SEIRModel,1500);
    try
        KalHess = Parameters.KalHess;
        Test = Parameters.KalmMaxReached;
    end
    NbIts = NbIts + 1;
    disp(NbIts)
    if NbIts>50
        return
    end
end


TempRess = {};

NbIters = 150000;
dim = 7;

load([SavePath '/ResTrueCov.mat']);


switch ind
    case 1
        % 1: from posterior cov, w-o eps 
        Cov = cov(Res.TransfThetas(:,100000:end)');
        Parameters.G = Cov^-1;
        Parameters.NoPaths = 1;
        Parameters.NbVariables = 7;
        Parameters.aim = 0.23;
        Parameters.Epsil = 1;
        Parameters.AdaptC = 0;
        TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
        Parameters.ModelType='SMC';
        Parameters.AdMet = 0;
        if indeps == 1
            % no eps, no AM
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 2
            % eps, no AM
            Parameters.AdaptC = 0.99;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 3
            % AM, no eps
            Parameters.AdaptC = 0;
            Parameters.AdMet = 1;
            Parameters.AdMetBeta = 0.05;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 4
            % eps and AM
            Parameters.AdaptC = 0.99;
            Parameters.AdMet = 1;
            Parameters.AdMetBeta = 0.05;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        end

    case 2
        % 2: from EKF mode cov, w-o eps
        Cov = (-KalHess)^-1;
        Parameters.G = Cov^-1;
        Parameters.NoPaths = 1;
        Parameters.NbVariables = 7;
        Parameters.aim = 0.23;
        Parameters.Epsil = 1;
        Parameters.AdaptC = 0;
        TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
        Parameters.ModelType='SMC';
        
        if indeps == 1
            % no eps, no AM
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 2
            % eps, no AM
            Parameters.AdaptC = 0.99;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 3
            % AM, no eps
            Parameters.AdaptC = 0;
            Parameters.AdMet = 1;
            Parameters.AdMetBeta = 0.05;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 4
            % eps and AM
            Parameters.AdaptC = 0.99;
            Parameters.AdMet = 1;
            Parameters.AdMetBeta = 0.05;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        end

    case 3
        % 3: from EKF post cov, w-o eps
        %%%EKF - MCMC
        Cov = (-KalHess)^-1;
        Parameters.G = Cov^-1;
        Parameters.NoPaths = 1;
        Parameters.ModelType='SMC';
        Parameters.AdaptC = 0.99;
        Parameters.NbVariables = 7;
        Parameters.aim = 0.23;
        Parameters.Epsil = 1;
        TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
        Parameters.ModelType='Kalman';
        Parameters.AdaptC = 0.99;
        ResKal = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,20000);
        Cov = cov(ResKal.TransfThetas');
        Parameters.G = Cov^-1;
        Parameters.ModelType='Kalman';
        Parameters.AdaptC = 0.99;
        ResKal = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,20000);
        save([SavePath 'TestingDifferentCovs_ResKal_' num2str(ind) '_' num2str(indeps) '.mat'], 'ResKal');

        
        %%%Go
        Cov = cov(ResKal.TransfThetas');
        Parameters.G = Cov^-1;
        Parameters.NoPaths = 1;
        Parameters.NbVariables = 7;
        Parameters.aim = 0.23;
        Parameters.Epsil = 1;
        Parameters.AdaptC = 0;
        Parameters.ModelType='SMC';
        TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
        
        if indeps == 1
            % no eps, no AM
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 2
            % eps, no AM
            Parameters.AdaptC = 0.99;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 3
            % AM, no eps
            Parameters.AdaptC = 0;
            Parameters.AdMet = 1;
            Parameters.AdMetBeta = 0.05;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 4
            % eps and AM
            Parameters.AdaptC = 0.99;
            Parameters.AdMet = 1;
            Parameters.AdMetBeta = 0.05;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        end

    case 4
        % 4: from Identity, w-o eps
        Cov = eye(dim);
        Parameters.G = Cov^-1;
        Parameters.NoPaths = 1;
        Parameters.NbVariables = 7;
        Parameters.aim = 0.23;
        Parameters.Epsil = 1;
        Parameters.AdaptC = 0;
        Parameters.ModelType='SMC';
        TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
        
        if indeps == 1
            % no eps, no AM
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 2
            % eps, no AM
            Parameters.AdaptC = 0.99;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 3
            % AM, no eps
            Parameters.AdaptC = 0;
            Parameters.AdMet = 1;
            Parameters.AdMetBeta = 0.05;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 4
            % eps and AM
            Parameters.AdaptC = 0.99;
            Parameters.AdMet = 1;
            Parameters.AdMetBeta = 0.05;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        end
    case 5
        % 1: from diagonal posterior cov, 
        Cov = diag(diag(cov(Res.TransfThetas')));
        Parameters.G = Cov^-1;
        Parameters.NoPaths = 1;
        Parameters.NbVariables = 7;
        Parameters.aim = 0.23;
        Parameters.Epsil = 1;
        Parameters.AdaptC = 0;
        TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
        Parameters.ModelType='SMC';
        Parameters.AdMet = 0;

        if indeps == 1
            % no eps, no AM
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 2
            % eps, no AM
            Parameters.AdaptC = 0.99;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 3
            % AM, no eps
            Parameters.AdaptC = 0;
            Parameters.AdMet = 1;
            Parameters.AdMetBeta = 0.05;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        elseif indeps == 4
            % eps and AM
            Parameters.AdaptC = 0.99;
            Parameters.AdMet = 1;
            Parameters.AdMetBeta = 0.05;
            Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,NbIters);
            save([SavePath 'TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat'], 'Res');
        end
    
end


Ress = {};
for ind = 1:5
    for indeps = 1:4
        load([SavePath '/TestingDifferentCovs_' num2str(ind) '_' num2str(indeps) '.mat']);
        Ress{ind,indeps} = Res;
    end
end
% % %         
for i = 1:5
    disp(num2str(i))
    disp('acc rate')
    disp([Ress{i,1}.AccRate   Ress{i,2}.AccRate  Ress{i,3}.AccRate   Ress{i,4}.AccRate]) 
    disp('min ESS transf')
    disp([min(Ress{i,1}.RelESSTransf)   min(Ress{i,2}.RelESSTransf)  min(Ress{i,3}.RelESSTransf)   min(Ress{i,4}.RelESSTransf)])
%     disp('min ESS ')
%     disp([min(Ress{i,1}.RelESS)   min(Ress{i,2}.RelESS)  min(Ress{i,3}.RelESS)   min(Ress{i,4}.RelESS)])
end
% 
% inds = [1 500 1000 5000 10000 15000 20000 25000];
% vect = [];
% for i = 1:length(inds)
%     temp = AutoCorrelation(Ress{1,3}.TransfThetas(4,inds(i):size(Res.TransfThetas,2)-max(inds)+inds(i)),1000);
%     vect(i) = 1/(1+2*sum(temp(2:end)));
% end
% vect
% 
% % 
for i = 1:5
    for j = 1:4
        [i j]
        Ress{i,j}.AutoCorr = {};
        for k = 1:size(Ress{i,j}.TransfThetas,1)
            temp = AutoCorrelation(Ress{i,j}.TransfThetas(k,50000:end));
            Ress{i,j}.ESSTransf(k) = 100000/(1+2*sum(temp(2:end)));
            Ress{i,j}.RelESSTransf(k) = 1/(1+2*sum(temp(2:end)))*100;
            Ress{i,j}.AutoCorr{k} = temp;
        end
    end
end
% 
% 
% i = 5
% j = 1
% plot(Ress{i,j}.AutoCorr{1})
% hold on
% for k = 2:size(Ress{i,j}.TransfThetas,1)
%     plot(Ress{i,j}.AutoCorr{k})
% end
% hold off




% 
% clf
% cols = rand(7,3);
% for i = 1:2
% %     subplot(2,1,i)
%     for k = 1:size(Ress{1,i}.TransfThetas,1)
%         temp = AutoCorrelation(Ress{1,i}.TransfThetas(k,:),500);
%         if i == 1
%             plot(temp,'col',cols(k,:))
%         else
%             plot(temp,'--','col',cols(k,:))
%         end
%         hold on
%     end
% end
% hold off
%         

