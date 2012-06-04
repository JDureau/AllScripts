function [] = PostThreeMethods(Constr)

Med = 0; % use median trajectory instead of mean

%ModelType = 'mInf';
ModelType = '';%Bert
%ModelType = '10';
% ModelType = 'Affine'; 
%  ModelType = 'Step';
% ModelType = 'Logist10';
%ModelType = 'Sigm';


ESSthr = 100;

SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';
cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])




if strcmp(ModelType,'')
    load([SavePath '/ResGenBR_v2.mat'])    
elseif strcmp(ModelType,'10')
    load([SavePath '/ResGenBR10.mat'])
elseif strcmp(ModelType,'mInf')
    load([SavePath '/ResGenBRmInf.mat'])
elseif strcmp(ModelType,'Logist')
    load([SavePath '/ResGenLogist.mat'])
elseif strcmp(ModelType,'Sigm')
    load([SavePath '/ResGenSigm.mat'])
elseif strcmp(ModelType,'Step10')
    load([SavePath '/ResGenStep10.mat'])
elseif strcmp(ModelType,'Step')
    load([SavePath '/ResGenStep.mat'])
elseif strcmp(ModelType,'Affine')
    load([SavePath '/ResGenAffine.mat'])
elseif strcmp(ModelType,'Affine10')
    load([SavePath '/ResGenAffine10.mat'])
end




ampls = [];
partampls2003 = [];
partampls1995 = [];
asympts = [];
baselines = [];
amplsALL = [];
partampls2003ALL = [];
partampls1995ALL = [];
asymptsALL = [];
baselinesALL = [];
sigms = [];
amplsROC = [];
partampls2003ROC = [];
partampls1995ROC = [];
asymptsROC = [];
baselinesROC = [];
amplsPost = [];
partampls2003Post = [];
partampls1995Post = [];
asymptsPost = [];
baselinesPost = [];
cpt = 0;
NbSamples = [];
ParsMeans = [];
ParsModes = [];
ParsMeds = [];
Pars2p5 = [];
Pars97p5 = [];

for i = 1:100

try   
   die = 0;
       try
       if Constr
	 if strcmp(ModelType,'')
          load([SavePath 'ResGenBRBert_BRdetermConstr_' num2str(i) '.mat'])
	elseif strcmp(ModelType,'10')
          load([SavePath 'ResGenBRBert10_BRdetermConstr_' num2str(i) '.mat']) 
	elseif strcmp(ModelType,'mInf')
          load([SavePath 'ResGenBRBertmInf_BRdetermConstr_' num2str(i) '.mat'])
	else
	   load([SavePath 'ResGenBR' ModelType  '_BRdetermConstr_' num2str(i) '.mat'])
	 end
       else
	 if strcmp(ModelType,'')
          load([SavePath 'ResGenBRBert_BRdeterm_' num2str(i) '.mat'])
        elseif strcmp(ModelType,'10')
          load([SavePath 'ResGenBRBert10_BRdeterm_' num2str(i) '.mat'])
	elseif strcmp(ModelType,'mInf')
          load([SavePath 'ResGenBRBertmInf_BRdeterm_' num2str(i) '.mat'])
	else
           load([SavePath 'ResGenBR' ModelType  '_BRdeterm_' num2str(i) '.mat'])
         end	
       end
	catch
	  disp('Det crash')
	  die = 1;
	end
% load([SavePath 'ResGenBRBert10_BRdeterm_' num2str(i) '.mat'])
%          load([SavePath 'ResGenBRBert_BRdeterm_' num2str(i) '.mat'])
% load([SavePath 'ResGenBR' ModelType  '_BRdeterm_' num2str(i) '.mat'])
%load([SavePath 'ResGenBRBert10_BRdeterm_' num2str(i) '.mat'])
	ResDet = Res;

    try
        if strcmp(ModelType,'')
            load([SavePath 'ResGenBRBert_Sigmdeterm_' num2str(i) '.mat'])
        elseif strcmp(ModelType,'10')
            load([SavePath 'ResGenBRBert10_BRdeterm_' num2str(i) '.mat'])
        elseif strcmp(ModelType,'mInf')
            load([SavePath 'ResGenBRBertmInf_BRdeterm_' num2str(i) '.mat'])
        else
            load([SavePath 'ResGenBR' ModelType  '_BRdeterm_' num2str(i) '.mat'])
        end
    catch
        disp('Det crash')
        die = 1;
	end
% load([SavePath 'ResGenBRBert10_BRdeterm_' num2str(i) '.mat'])
%          load([SavePath 'ResGenBRBert_BRdeterm_' num2str(i) '.mat'])
% load([SavePath 'ResGenBR' ModelType  '_BRdeterm_' num2str(i) '.mat'])
%load([SavePath 'ResGenBRBert10_BRdeterm_' num2str(i) '.mat'])
	ResSigm = Res;


        try
	  if Constr
	    load([SavePath 'VIH_PlayingWithSigm' ModelType '_ResultsAddConstr__1_' num2str(i) '.mat'])
	  else
	    load([SavePath 'VIH_PlayingWithSigm' ModelType '_ResultsAdd__1_' num2str(i) '.mat'])
          end
	catch
	  die = 1;
          disp('Add crash')
        end
        ResAdd = Ress{1};

	try
	  if Constr
	    load([SavePath 'VIH_PlayingWithSigm' ModelType '_ResultsBertallanfyConstr__1_' num2str(i) '.mat'])
          else
	    load([SavePath 'VIH_PlayingWithSigm' ModelType '_ResultsBertallanfy__1_' num2str(i) '.mat'])
	  end
        catch
	  die = 1;
          disp('BR crash')
        end
        ResBer = Ress{1};
        
        
        try
	  if Constr
	    load([SavePath 'VIH_PlayingWithSigm' ModelType '_ResultsSigmoidConstr__1_' num2str(i) '.mat'])
          else
	    load([SavePath 'VIH_PlayingWithSigm' ModelType '_ResultsSigmoid__1_' num2str(i) '.mat'])
	  end
        catch
	  die = 1;
          disp('BR crash')
        end
        ResSigmSto = Ress{1};
        


	if or(mean(not(ResDet.Data.BuiltTraj(1:585,9) == ResSigm.Data.BuiltTraj(1:585,9))),mean(not(ResDet.Data.BuiltTraj(1:585,9) == ResBer.Data.BuiltTraj(1:585,9))))
	  disp('Careful Different ResGens..')
	  creve
	end

	vect = [];
	Names = ResDet.Parameters.Names.Estimated;
	for j = 1:length(Names)
          v = ResDet.Thetas(ResDet.Parameters.(Names{j}).Index,:);
          tmp = AutoCorrelation(v,100);
          vect(j) = length(v)/(1+2*sum(tmp(2:end)));
	end
	NbSamples(cpt+1,1) = round(size(ResDet,2)/5);;
	vect = [];
        Names = ResBer.Parameters.Names.Estimated;
        for j = 1:length(Names)
          v = ResBer.Thetas(ResBer.Parameters.(Names{j}).Index,:);
          tmp = AutoCorrelation(v,100);
          vect(j) = length(v)/(1+2*sum(tmp(2:end)));
        end
	NbSamples(cpt+1,3) = round(size(ResBer,2)/5);;
	vect = [];
        Names = ResAdd.Parameters.Names.Estimated;
        for j = 1:length(Names)
          v = ResAdd.Thetas(ResAdd.Parameters.(Names{j}).Index,:);
          tmp = AutoCorrelation(v,100);
          vect(j) = length(v)/(1+2*sum(tmp(2:end)));
        end
	NbSamples(cpt+1,2) = round(size(ResAdd,2)/5);
    vect = [];
        Names = ResSigm.Parameters.Names.Estimated;
        for j = 1:length(Names)
          v = ResSigm.Thetas(ResSigm.Parameters.(Names{j}).Index,:);
          tmp = AutoCorrelation(v,100);
          vect(j) = length(v)/(1+2*sum(tmp(2:end)));
        end
	NbSamples(cpt+1,5) = round(size(ResSigm,2)/5);
    vect = [];
        Names = ResSigmSto.Parameters.Names.Estimated;
        for j = 1:length(Names)
          v = ResSigm.StoThetas(ResSigmSto.Parameters.(Names{j}).Index,:);
          tmp = AutoCorrelation(v,100);
          vect(j) = length(v)/(1+2*sum(tmp(2:end)));
        end
	NbSamples(cpt+1,6) = round(size(ResSigm,2)/5);
    
    

	if or(ResAdd.AccRate<0.10,ResAdd.AccRate>0.50)
	  die = 1;
	end
	if or(ResBer.AccRate<0.10,ResBer.AccRate>0.50)
          die =1;
        end

	if die
	  'final'
	  creve
	end




    cpt = cpt+1;
    [cpt i]
    if Med
      FtDet = median(squeeze(ResDet.Paths(:,3,1:585)));
    else
      FtDet = mean(squeeze(ResDet.Paths(:,3,1:585)));
    end
    FtsDet = (squeeze(ResDet.Paths(:,3,1:585)));
    
    if Med
      FtSigm = median(squeeze(ResSigm.Paths(:,3,1:585)));
    else
      FtSigm = mean(squeeze(ResSigm.Paths(:,3,1:585)));
    end
    FtsSigm = (squeeze(ResSigm.Paths(:,3,1:585)));

    tmp = squeeze(ResAdd.Paths(:,3,1:585));
    if Med
      FtAdd = median(exp(tmp)./(1+exp(tmp)));
    else
      FtAdd = mean(exp(tmp)./(1+exp(tmp)));
    end
    FtsAdd = (exp(tmp)./(1+exp(tmp)));

    tmp = squeeze(ResBer.Paths(:,3,1:585));
    beta0s = squeeze(ResBer.Thetas(ResBer.Parameters.BRbase.Index,:));
    mus = squeeze(ResBer.Thetas(ResBer.Parameters.BRmu.Index,:));
    ms = squeeze(ResBer.Thetas(ResBer.Parameters.BRmm1.Index,:))+1;
    ts = squeeze(ResBer.Thetas(ResBer.Parameters.BRtinfl.Index,:));
    Bs = 1-(beta0s./mus).^(1-ms);
    ks = 1./ts.*log(Bs./(1-ms));
    if Med
      FtBer = median((repmat(1-ms',1,585).*tmp+repmat(mus.^(1-ms)',1,585)).^(repmat(1./(1-ms)',1,585)));
    else
      FtBer = mean((repmat(1-ms',1,585).*tmp+repmat(mus.^(1-ms)',1,585)).^(repmat(1./(1-ms)',1,585)));
    end
    FtsBer = ((repmat(1-ms',1,585).*tmp+repmat(mus.^(1-ms)',1,585)).^(repmat(1./(1-ms)',1,585)));
    
    tmp = squeeze(ResSigmSto.Paths(:,3,1:indend));
    rate = squeeze(ResSigmSto.Thetas(ResSigmSto.Parameters.Sigmrate.Index,:));
    base = squeeze(ResSigmSto.Thetas(ResSigmSto.Parameters.Sigmbase.Index,:));
    mu = squeeze(ResSigmSto.Thetas(ResSigmSto.Parameters.Sigmmu.Index,:));
    tinfl = squeeze(ResSigmSto.Thetas(ResSigmSto.Parameters.Sigmtinfl.Index,:));;
    c = 1./(1+exp(tinfl./rate));
    b = (mu-base).*c./(1-c);
    a = base - b;
    a = repmat(a',1,indend);
    b = repmat(b',1,indend);
    c = repmat(c',1,indend);
    FtsSigmSto = (a + b./(c.*(1+tmp)));
    if Med
      FtSigmSto = median(FtsSigmSto);
    else
      FtSigmSto = mean(FtsSigmSto);
    end
    
    
    FtTrue = squeeze(ResBer.Data.BuiltTraj(1:585,9));
    ampls(1,cpt) = max(FtDet)-min(FtDet);
    ampls(2,cpt) = max(FtAdd)-min(FtAdd);
    ampls(3,cpt) = max(FtBer)-min(FtBer);
    ampls(5,cpt) = max(FtSigm)-min(FtSigm);
    ampls(6,cpt) = max(FtSigmSto)-min(FtSigmSto);
    ampls(4,cpt) = max(FtTrue)-min(FtTrue);
    partampls2003(1,cpt) = FtDet(584)-FtDet(end-168);
    partampls2003(2,cpt) = FtAdd(584)-FtAdd(end-168);
    partampls2003(3,cpt) = FtBer(584)-FtBer(end-168);
    partampls2003(5,cpt) = FtSigm(584)-FtSigm(end-168);
    partampls2003(6,cpt) = FtSigmSto(584)-FtSigmSto(end-168);
    partampls2003(4,cpt) = FtTrue(584)-FtTrue(end-168);
    partampls1995(1,cpt) = FtDet(584)-FtDet(end-360);
    partampls1995(2,cpt) = FtAdd(584)-FtAdd(end-360);
    partampls1995(3,cpt) = FtBer(584)-FtBer(end-360);
    partampls1995(5,cpt) = FtSigm(584)-FtSigm(end-360);
    partampls1995(6,cpt) = FtSigmSto(584)-FtSigmSto(end-360);
    partampls1995(4,cpt) = FtTrue(584)-FtTrue(end-360);
    baselines(1,cpt) = FtDet(1);
    baselines(2,cpt) = FtAdd(1); 
    baselines(3,cpt) = FtBer(1);
    baselines(5,cpt) = FtSigm(1);
    baselines(6,cpt) = FtSigmSto(1);
    baselines(4,cpt) = FtTrue(1);   
    asympts(1,cpt) = FtDet(584);
    asympts(2,cpt) = FtAdd(584);
    asympts(3,cpt) = FtBer(584);
    asympts(5,cpt) = FtSigm(584);
    asympts(6,cpt) = FtSigmSto(584);
    asympts(4,cpt) = FtTrue(584);
   

    sigms(1,cpt,:) = FtDet;
    sigms(2,cpt,:) = FtAdd;
    sigms(3,cpt,:) = FtBer;
    sigms(5,cpt,:) = FtSigm;
    sigms(6,cpt,:) = FtSigmSto;
    sigms(4,cpt,:) = FtTrue;



    NPartsDet = size(FtsDet,1);
    NPartsAdd = size(FtsAdd,1);  
    NPartsBer = size(FtsBer,1);
    NPartsSigm = size(FtsSigm,1);
    NPartsSigmSto = size(FtsSigmSto,1);
    NbSamples(cpt,:) = [NPartsDet NPartsAdd NPartsBer 0 NPartsSigm NPartsSigmSto];
    amplsALL{1}(cpt,1:NPartsDet) = max(FtsDet')-min(FtsDet');
    amplsALL{2}(cpt,1:NPartsAdd) = max(FtsAdd')-min(FtsAdd');
    amplsALL{3}(cpt,1:NPartsBer) = max(FtsBer')-min(FtsBer');
    amplsALL{5}(cpt,1:NPartsSigm) = max(FtsSigm')-min(FtsSigm');
    amplsALL{6}(cpt,1:NPartsSigmSto) = max(FtsSigmSto')-min(FtsSigmSto');
    partampls2003ALL{1}(cpt,1:NPartsDet) = FtsDet(:,585)-FtsDet(:,585-168);
    partampls2003ALL{2}(cpt,1:NPartsAdd) = FtsAdd(:,585)-FtsAdd(:,585-168);
    partampls2003ALL{3}(cpt,1:NPartsBer) = FtsBer(:,585)-FtsBer(:,585-168);
    partampls2003ALL{5}(cpt,1:NPartsSigm) = FtsSigm(:,585)-FtsSigm(:,585-168);
    partampls2003ALL{6}(cpt,1:NPartsSigmSto) = FtsSigmSto(:,585)-FtsSigmSto(:,585-168);
    partampls1995ALL{1}(cpt,1:NPartsDet) = FtsDet(:,585)-FtsDet(:,585-360);
    partampls1995ALL{2}(cpt,1:NPartsAdd) = FtsAdd(:,585)-FtsAdd(:,585-360);
    partampls1995ALL{3}(cpt,1:NPartsBer) = FtsBer(:,585)-FtsBer(:,585-360);
    partampls1995ALL{5}(cpt,1:NPartsSigm) = FtsSigm(:,585)-FtsSigm(:,585-360);
    partampls1995ALL{6}(cpt,1:NPartsSigmSto) = FtsSigmSto(:,585)-FtsSigmSto(:,585-360);
    baselinesALL{1}(cpt,1:NPartsDet) = FtsDet(:,1);
    baselinesALL{2}(cpt,1:NPartsAdd) = FtAdd(:,1); 
    baselinesALL{3}(cpt,1:NPartsBer) = FtBer(:,1);  
    baselinesALL{5}(cpt,1:NPartsSigm) = FtSigm(:,1);  
    baselinesALL{6}(cpt,1:NPartsSigmSto) = FtSigmSto(:,1);  
    asymptsALL{1}(cpt,1:NPartsDet) = FtsDet(:,585);
    asymptsALL{2}(cpt,1:NPartsAdd) = FtsAdd(:,585);
    asymptsALL{3}(cpt,1:NPartsBer) = FtsBer(:,585);
    asymptsALL{5}(cpt,1:NPartsSigm) = FtsSigm(:,585);
    asymptsALL{6}(cpt,1:NPartsSigmSto) = FtsSigmSto(:,585);

    qs = 0.01:0.01:1; 

    for j = 1:length(qs)    
      amplsROC(1,j,cpt) = quantile(amplsALL{1}(cpt,1:NPartsDet),qs(j));
      amplsROC(2,j,cpt) = quantile(amplsALL{2}(cpt,1:NPartsAdd),qs(j));
      amplsROC(3,j,cpt) = quantile(amplsALL{3}(cpt,1:NPartsBer),qs(j));
      amplsROC(5,j,cpt) = quantile(amplsALL{5}(cpt,1:NPartsSigm),qs(j));
      amplsROC(6,j,cpt) = quantile(amplsALL{6}(cpt,1:NPartsSigmSto),qs(j));
      amplsROC(4,j,cpt) = max(FtTrue)-min(FtTrue);
      partampls2003ROC(1,j,cpt) = quantile(partampls2003ALL{1}(cpt,1:NPartsDet),qs(j));
      partampls2003ROC(2,j,cpt) = quantile(partampls2003ALL{2}(cpt,1:NPartsAdd),qs(j));
      partampls2003ROC(3,j,cpt) = quantile(partampls2003ALL{3}(cpt,1:NPartsBer),qs(j));
      partampls2003ROC(5,j,cpt) = quantile(partampls2003ALL{5}(cpt,1:NPartsSigm),qs(j));
      partampls2003ROC(6,j,cpt) = quantile(partampls2003ALL{6}(cpt,1:NPartsSigmSto),qs(j));
      partampls2003ROC(4,j,cpt) = FtTrue(end)-FtTrue(end-168);
      partampls1995ROC(1,j,cpt) = quantile(partampls1995ALL{1}(cpt,1:NPartsDet),qs(j));
      partampls1995ROC(2,j,cpt) = quantile(partampls1995ALL{2}(cpt,1:NPartsAdd),qs(j));
      partampls1995ROC(3,j,cpt) = quantile(partampls1995ALL{3}(cpt,1:NPartsBer),qs(j));
      partampls1995ROC(5,j,cpt) = quantile(partampls1995ALL{5}(cpt,1:NPartsSigm),qs(j));
      partampls1995ROC(6,j,cpt) = quantile(partampls1995ALL{6}(cpt,1:NPartsSigmSto),qs(j));
      partampls1995ROC(4,j,cpt) = FtTrue(end)-FtTrue(end-360);
      baselinesROC(1,j,cpt) = quantile(baselinesALL{1}(cpt,1:NPartsDet),qs(j));
      baselinesROC(2,j,cpt) = quantile(baselinesALL{2}(cpt,1:NPartsAdd),qs(j));
      baselinesROC(3,j,cpt) = quantile(baselinesALL{3}(cpt,1:NPartsBer),qs(j));
      baselinesROC(5,j,cpt) = quantile(baselinesALL{5}(cpt,1:NPartsSigm),qs(j));
      baselinesROC(6,j,cpt) = quantile(baselinesALL{6}(cpt,1:NPartsSigmSto),qs(j));
      baselinesROC(4,j,cpt) = FtTrue(1);
      asymptsROC(1,j,cpt) = quantile(asymptsALL{1}(cpt,1:NPartsDet),qs(j));
      asymptsROC(2,j,cpt) = quantile(asymptsALL{2}(cpt,1:NPartsAdd),qs(j));
      asymptsROC(3,j,cpt) = quantile(asymptsALL{3}(cpt,1:NPartsBer),qs(j));
      asymptsROC(5,j,cpt) = quantile(asymptsALL{5}(cpt,1:NPartsSigm),qs(j));
      asymptsROC(6,j,cpt) = quantile(asymptsALL{6}(cpt,1:NPartsSigmSto),qs(j));
      asymptsROC(4,j,cpt) = FtTrue(end);
    end
    
    Parameters = ResDet.Parameters;
    Names = Parameters.Names.Estimated;
    for j = 1:length(Names)
        inddet = ResDet.Parameters.(Names{j}).Index;
        indadd = ResAdd.Parameters.(Names{j}).Index;
        indber = ResBer.Parameters.(Names{j}).Index;
        indSigm = ResSigm.Parameters.(Names{j}).Index;
        indSigmSto = ResSigmSto.Parameters.(Names{j}).Index;
% 	disp(Names{j})
%         disp(inddet)
% 	disp(indber)

        ParsMeans(1,inddet,cpt) = mean(ResDet.Thetas(inddet,:));
        ParsMeds(1,inddet,cpt) = median(ResDet.Thetas(inddet,:));
        Pars2p5(1,inddet,cpt) = quantile(ResDet.Thetas(inddet,:),0.025);
        Pars97p5(1,inddet,cpt) = quantile(ResDet.Thetas(inddet,:),0.975);
        try
            ParsMeans(2,inddet,cpt) = mean(ResAdd.Thetas(indadd,:));
            ParsMeds(2,inddet,cpt) = median(ResAdd.Thetas(indadd,:));
            Pars2p5(2,inddet,cpt) = quantile(ResAdd.Thetas(indadd,:),0.025);
            Pars97p5(2,inddet,cpt) = quantile(ResAdd.Thetas(indadd,:),0.975);
       catch
            ParsMeans(2,inddet,cpt) = NaN;
            ParsMeds(2,inddet,cpt) = NaN;
            Pars2p5(2,inddet,cpt) = NaN;
            Pars97p5(2,inddet,cpt) = NaN;
        end
        try
            ParsMeans(5,inddet,cpt) = mean(ResSigm.Thetas(indadd,:));
            ParsMeds(5,inddet,cpt) = median(ResSigm.Thetas(indadd,:));
            Pars2p5(5,inddet,cpt) = quantile(ResSigm.Thetas(indadd,:),0.025);
            Pars97p5(5,inddet,cpt) = quantile(ResSigm.Thetas(indadd,:),0.975);
        catch
            ParsMeans(5,inddet,cpt) = NaN;
            ParsMeds(5,inddet,cpt) = NaN;
            Pars2p5(5,inddet,cpt) = NaN;
            Pars97p5(5,inddet,cpt) = NaN;
        end
        try
            ParsMeans(6,inddet,cpt) = mean(ResSigmSto.Thetas(indadd,:));
            ParsMeds(6,inddet,cpt) = median(ResSigmSto.Thetas(indadd,:));
            Pars2p5(6,inddet,cpt) = quantile(ResSigmSto.Thetas(indadd,:),0.025);
            Pars97p5(6,inddet,cpt) = quantile(ResSigmSto.Thetas(indadd,:),0.975);
        catch
            ParsMeans(6,inddet,cpt) = NaN;
            ParsMeds(6,inddet,cpt) = NaN;
            Pars2p5(6,inddet,cpt) = NaN;
            Pars97p5(6,inddet,cpt) = NaN;
        end
        ParsMeans(3,inddet,cpt) = mean(ResBer.Thetas(indber,:));
        ParsMeds(3,inddet,cpt) = median(ResBer.Thetas(indber,:));
        Pars2p5(3,inddet,cpt) = quantile(ResBer.Thetas(indber,:),0.025);
        Pars97p5(3,inddet,cpt) = quantile(ResBer.Thetas(indber,:),0.975);
        ParsMeans(4,inddet,cpt) = ResGens{cpt}.Parameters.(Names{j}).Value;
        ParsMeds(4,inddet,cpt) = ResGens{cpt}.Parameters.(Names{j}).Value;
        Pars2p5(4,inddet,cpt) = ResGens{cpt}.Parameters.(Names{j}).Value;
        Pars97p5(4,inddet,cpt) = ResGens{cpt}.Parameters.(Names{j}).Value;
    end
%     Pars
end
end

Res = struct();
Res.qs = qs;
Res.Parameters = Parameters;
Res.ampls = ampls;
Res.sigms = sigms;
Res.baselines = baselines;
Res.asympts = asympts;
Res.partampls1995 = partampls1995;
Res.partampls2003 = partampls2003;
%Res.amplsALL = amplsALL;
%Res.baselinesALL = baselinesALL;
%Res.asymptsALL = asymptsALL;
%Res.partampls1995ALL = partampls1995ALL;
%Res.partampls2003ALL = partampls2003ALL;
Res.amplsROC = amplsROC;
Res.baselinesROC = baselinesROC;
Res.asymptsROC = asymptsROC;
Res.partampls1995ROC = partampls1995ROC;
Res.partampls2003ROC = partampls2003ROC;
Res.amplsPost = amplsPost;
Res.baselinesPost = baselinesPost;
Res.asymptsPost = asymptsPost;
Res.partampls1995Post = partampls1995Post;
Res.partampls2003Post = partampls2003Post;
Res.NbSamples = NbSamples;
Res.ParsMeans = ParsMeans;
Res.ParsMeds = ParsMeds;
Res.Pars2p5 = Pars2p5;
Res.Pars97p5 = Pars97p5;
if Constr
  if Med
    save([SavePath 'PostThreeMethodsConstr' ModelType '_Med.mat'],'Res')
  else
    save([SavePath 'PostThreeMethodsConstr' ModelType '.mat'],'Res')
  end
else
  if Med
    save([SavePath 'PostThreeMethods' ModelType '_Med.mat'],'Res')  
  else
    save([SavePath 'PostThreeMethods' ModelType '.mat'],'Res')  
  end
  [SavePath 'PostThreeMethods' ModelType '.mat']
  
end  