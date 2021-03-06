function [] = PostTreating()

cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])
SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';

Resss = {};

MeanErrors = [];
Errors = [];
amplin = [];
amplout = [];
slope = [];
inflptin = [];
inflptout = [];
corrs = [];
baselinein = [];
baselineout = [];
endlinein = [];
endlineout = [];
sigmsin = [];
sigmsout = [];
%propaftintin = [];
%propaftintout = [];
partafterin = [];
partafterout = [];
PrevObs = {};
InflPts = {};
Ampls = {};
valsinbeg = [];
valsinend = [];
valsoutbeg = [];
valsoutend = [];
for BigInd = 1:7
    for SmallInd = 1:10
      disp(['SmallInd=' num2str(SmallInd)])
      if not(and(SmallInd==2,BigInd==100))
           load([SavePath 'VIH_PlayingWithSigm_ResultsLogistic__Biased_' num2str(BigInd) '_' num2str(SmallInd) '.mat'])          
         %   load([SavePath '/VIH_PlayingWithSigm_ResultsLogistic__' num2str(BigInd) '_' num2str(SmallInd) '.mat'])
            disp(['SmallInd=' num2str(SmallInd)])
	    for j = 1:length(Ress)
                Res = Ress{j};
                Parameters = Res.Parameters;
                j
                disp(Res.AccRate)
                try
                    if and(Res.AccRate > 0.18, Res.AccRate < 0.30)
                        tmp = diff(Res.Data.Fts);
                        [b,ind] = max(tmp);
                        if b<0.1*Parameters.ComputationTStep
                            indampl = Res.Parameters.CUdelta.Index;
                            indbase = Res.Parameters.CUinit.Index;

                            mu = Res.Thetas(indbase,:) + Res.Thetas(indampl,:); 
                            tmp = squeeze(Res.Paths(50:end,3,:));
                            tmp = diag(mu)*exp(tmp)./(1+exp(tmp));
                            if max(mean(tmp))-min(mean(tmp))>0.9
                                meurt
                            end
%                             Pars = FitSigmoids((tmp),Res.Parameters);
                            
                            %Pars = FitSigmoid(mean(tmp),Res.Parameters);
% 			    mu = Res.Parameters.CUinit.Value + Res.Parameters.CUdelta.Value;                            
%                             Betas = min(1,max(0,mu*exp(Res.Paths(:,3,:))./(1+exp(Res.Paths(:,3,:)))));
%                             
%                             indinfl = Res.Parameters.CUinfl.Index;
%                             indampl = Res.Parameters.CUdelta.Index;
%                             indbase = Res.Parameters.CUinit.Index;
%                             
%                             
%                             
%                             amplin(end+1) = Res.Data.Fts(585)-Res.Data.Fts(1)   ;
%                             amplout(end+1) = mean(Res.Thetas(indampl,:));           
%                             inflptin(end+1) = ind;
%                             baselineout(end+1) = mean(Res.Thetas(indbase,:));
%                             baselinein(end+1) = Res.Data.Fts(1);
%                             endlinein(end+1) = Res.Data.Fts(585);
%                             endlineout(end+1) = mean(Res.Thetas(indbase,:)) + mean(Res.Thetas(indampl,:));
%                             inflptout(end+1) = mean(Res.Thetas(indinfl,:));
%                             
%                             partafterin(end+1) = Res.Data.Fts(585)-Res.Data.Fts(456);
% %                              propaftintin(end+1) = partfter/(Res.Data.Fts(585)-Res.Data.Fts(1));
%                             valsinbeg(end+1) = Res.Data.Fts(384);
%                             valsinend(end+1) = Res.Data.Fts(585);
%                             valsoutbeg(end+1) = mean(Betas(:,384));
%                             valsoutend(end+1) = mean(Betas(:,384));
%                             partafterout(end+1) = mean(Betas(:,585))-mean(Betas(:,456));
% %                             propaftintout(end+1) = partfter/(Pars.Sigm(585)-Pars.Sigm(1));
%                             Res.Parameters.TypeWork='Boston Examples HIV2';
% %                             Res.Parameters.Sigm = Pars.Sigm; 
%                             PrevObs{end+1} = Res.Data.Observations; 
                            
                            %% BEFORE LOGISTIC:
                            amplin(end+1) = Res.Data.Fts(585)-Res.Data.Fts(1)   ;
                            amplout(end+1) = Pars.Ampl;           
                            inflptin(end+1) = ind;
                            baselineout(end+1) = Pars.Baseline;
                            baselinein(end+1) = Res.Data.Fts(1);
                            endlinein(end+1) = Res.Data.Fts(585);
                            endlineout(end+1) = Pars.Endline;
                            [b,indmed] = min(abs(Pars.Sigm-(Pars.Sigm(585)+Pars.Sigm(1))/2));
                            inflptout(end+1) = indmed;
                            sigmsin(end+1,:) = Res.Data.Fts;
                            sigmsout(end+1,:) = Pars.Sigm;
                            partfter = Res.Data.Fts(585)-Res.Data.Fts(456);
                            propaftintin(end+1) = partfter/(Res.Data.Fts(585)-Res.Data.Fts(1));
                            valsinbeg(end+1) = Res.Data.Fts(384);
			    valsinend(end+1) = Res.Data.Fts(585);
			    valsoutbeg(end+1) = Pars.Sigm(384);
			    valsoutend(end+1) = Pars.Sigm(585);
			    partfter = Pars.Sigm(585)-Pars.Sigm(456);
                            propaftintout(end+1) = partfter/(Pars.Sigm(585)-Pars.Sigm(1));
                            Res.Parameters.TypeWork='Boston Examples HIV2';
                            Res.Parameters.Sigm = Pars.Sigm; 
			    PrevObs{end+1} = Res.Data.Observations; 
                            
                        end
                    end
                end
            end
	end
	disp('clear')
	    clear('Ress')
        
        length(inflptout)
    end
end
% 
% BigInd = 1;
% SmallInd = 4;
% SavePath = '/Users/dureaujoseph/Documents/Taf/These/Results/Avahan/';
% load([SavePath 'VIH_PlayingWithSigm_RightPriors' num2str(BigInd) '_' num2str(SmallInd) '.mat'])
% j = 2;     
% Res = Ress{j};
% Res.Parameters.TypeWork='Boston Examples HIV2';
% PlotResHIV(Res,Res.Parameters)


Res = struct();
Res.amplin = amplin;
Res.amplout = amplout;
Res.inflptin = inflptin;
Res.inflptout = inflptout;
Res.baselinein = baselinein;
Res.baselineout = baselineout;
Res.endlinein = endlinein;
Res.endlineout = endlineout;
Res.sigmsin = sigmsin;
Res.sigmsout = sigmsout;
% Res.partafterin = partafterin;
% Res.partafterout = partafterout;
Res.propaftintin = propaftintin;
Res.propaftintout = propaftintout;
Res.PrevObs = PrevObs;
Res.valsinbeg = valsinbeg;
Res.valsinend = valsinend;
Res.valsoutbeg = valsoutbeg;
Res.valsoutend = valsoutend;
%Res.Ampls = Ampls;
%Res.InflPts = InflPts;
%save([SavePath '/TreatedRessRightPriors.mat'],'Res')
save([SavePath '/TreatedRessRightPriors_Biased.mat'],'Res')