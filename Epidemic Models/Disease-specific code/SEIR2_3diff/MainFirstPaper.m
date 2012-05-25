%% Main Paper 1
%%%%%%%%%%%%%%%%%%%%%%%%%

% This should have the script for all figures of paper 1

%% 
cd('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Filtering')
addpath('H:\My Documents\PhD Work\Matlab Scripts\General Tools')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes\Resampling\pf_resampling')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Joint Sampling')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\MIF')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Model Selection')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Optimization Approach')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\Models')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\SEIR')


%% Reverse-engineered tseries
% They already have been generated by MainSims
% load them (DatasGens):
SavePath = 'S:\Results\';
Name = [SavePath 'SEIR_simsDatas.mat'];

load(Name)

for i = 1:99
    plot(DatasGens{i}.RealBetaTraj)
    pause()
end

%% Apply All Methods to to Sims

Ress = {};
save([SavePath 'FirstPaperSimsRess.mat'],'Ress')

for i = 1:1
   TempName = [SavePath 'TempSimsRes.mat'];
   FullSEIRinference(DatasGens{i},'Add','Fixed',TempName);
   load(TempName);
   load([SavePath 'FirstPaperSimsRess.mat']);
   Ress{i} = Res3;
end