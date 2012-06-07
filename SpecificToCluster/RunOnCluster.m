function [] = RunOnCluster(BigInd,SmallInd,IndMode,IndDiff)

ModelType = 'Bert';
% ModelType = 'Step10';
%ModelType = 'Affine10';
% ModelType = 'Logist';
%ModelType = 'Sigm';


s = RandStream('mcg16807','Seed',100*(BigInd+SmallInd+3*SmallInd^2));
RandStream.setDefaultStream(s)

SmallInd = SmallInd + 1;
% SmallInd = SmallInd - floor(SmallInd/10)*10 + 1;

cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])

Modes = {'','_Biased','_Restricted','_AddingObs'};
Diffs = {'Add','Bertallanfy','Sigmoid','BertallanfyConstr'};

SavePath = '/users/ecologie/dureau/src/AllData/Avahan';
%SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan';


if strcmp(ModelType,'Bert')
    load([SavePath '/ResGenBR_v2.mat'])    
elseif strcmp(ModelType,'Bert10')
    load([SavePath '/ResGenBR10.mat'])
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
% load([SavePath '/VIH_PlayingWithSigm_ResGens' Modes{IndMode} '_'  num2str(BigInd) '.mat'])
%load(['VIH_PlayingWithSigm_ResGens' num2str(BigInd) '_Restricted.mat'])

inds = 1*(SmallInd-1)+1:1*(SmallInd);
Ress = {};
for i = inds(1):inds(end)
    
       
    if or(IndMode == 1, IndMode == 4)
        Parameters = ResGens{i}.Parameters;
    elseif IndMode == 2
        Parameters = ResGens{i}.BiasedParameters;
    elseif IndMode == 3
        Parameters = ResGens{i}.RestrictedParameters;
    end




    if IndDiff == 1
        Parameters.DiffusionType = 'Add';
    elseif IndDiff == 2
        Parameters.DiffusionType = 'Bertallanfy';
    elseif IndDiff == 3
        Parameters.DiffusionType = 'Sigmoid';
    elseif IndDiff == 4
        Parameters.DiffusionType = 'BertallanfyConstr';
    else
        'unknown diffusion'
        die
    end
    
    Parameters.DiffusionType
    Parameters.TempName = ['Temp_' ModelType '_' Parameters.DiffusionType '_' num2str(i) '.mat'];
%     Temp = struct();
%     save([SavePath '/' Parameters.TempName],'Temp')
%     if strcmp(ModelType,'Bert')
%          Parameters.NameToSave = [SavePath '/VIH_PlayingWithSigm_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'];
%     elseif strcmp(ModelType,'Logist')
%          Parameters.NameToSave = [SavePath '/VIH_PlayingWithSigmLogist_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'];
%     elseif strcmp(ModelType,'Sigm')
%          Parameters.NameToSave = [SavePath '/VIH_PlayingWithSigmSigm_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'];
%       end
    Parameters.RealData = 0 ;   
    Parameters.NameToSave = [SavePath '/VIH_PlayingWithSigm_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'];
    
    Res = HIVapplyInference(ResGens{i}.Data,Parameters);
    
    Ress{end+1} = Res;
    if strcmp(ModelType,'Bert')
        save([SavePath '/VIH_PlayingWithSigm_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'],'Ress')   
    elseif strcmp(ModelType,'Bert10')
        save([SavePath '/VIH_PlayingWithSigm10_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'],'Ress')
    elseif strcmp(ModelType,'Logist')
        save([SavePath '/VIH_PlayingWithSigmLogist_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'],'Ress')   
    elseif strcmp(ModelType,'Sigm')
        save([SavePath '/VIH_PlayingWithSigmSigm_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'],'Ress')   
    elseif strcmp(ModelType,'Step')
        save([SavePath '/VIH_PlayingWithSigmStep_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'],'Ress')
    elseif strcmp(ModelType,'Step10')
        save([SavePath '/VIH_PlayingWithSigmStep10_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'],'Ress')
    elseif strcmp(ModelType,'Affine')
      save([SavePath '/VIH_PlayingWithSigmAffine_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'],'Ress')
    elseif strcmp(ModelType,'Affine10')
      save([SavePath '/VIH_PlayingWithSigmAffine10_Results' Diffs{IndDiff} '_' Modes{IndMode} '_' num2str(BigInd) '_' num2str(SmallInd) '.mat'],'Ress')
    end
end   

