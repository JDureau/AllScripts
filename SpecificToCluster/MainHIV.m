function [] =MainHIV(ind,IndDiff,NbRounds,TakeClients)

ind = ind+1;

% Main HIV clean

%% load paths


mode = '';

%  / _Mysore / _MysoreMu / _MysoreCf / _CU20 / _CU20_t2003

cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])

%cd('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\')
%addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\General Tools')
%addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Disease-specific code')
%addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Disease-specific code\HIV')
%addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Toolboxes\Resampling\pf_resampling')
%addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Generic PMCMC tools')
%addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Parameter Estimation')
%addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Toolboxes')


%% Go


SavePath = '/users/ecologie/dureau/src/AllData/Avahan/';

temp = load([SavePath 'ParametersMysore' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersMysore = temp.Parameters;


temp = load([SavePath 'ParametersBelgaum' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersBelgaum = temp.Parameters;

temp = load([SavePath 'ParametersBellary' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersBellary = temp.Parameters;

temp = load([SavePath 'ParametersEastGodavry' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersEastGodavry = temp.Parameters;

temp = load([SavePath 'ParametersGuntur' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersGuntur = temp.Parameters;

temp = load([SavePath 'ParametersHyderabad' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersHyderabad = temp.Parameters;

temp = load([SavePath 'ParametersVisag' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersVisag = temp.Parameters;

temp = load([SavePath 'ParametersWarangal' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersWarangal = temp.Parameters;

temp = load([SavePath 'ParametersYevatmal' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersYevatmal = temp.Parameters;

temp = load([SavePath 'ParametersShimoga' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersShimoga = temp.Parameters;

temp = load([SavePath 'ParametersChennai' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersChennai = temp.Parameters;

temp = load([SavePath 'ParametersMadurai' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersMadurai = temp.Parameters;

temp = load([SavePath 'ParametersSalem' mode '.mat']);
if IndDiff == 1
    temp.Parameters.DiffusionType = 'Add';
elseif IndDiff == 2
    temp.Parameters.DiffusionType = 'Bertallanfy';
elseif IndDiff == 3
    temp.Parameters.DiffusionType = 'Sigmoid';
end
ParametersSalem = temp.Parameters;





ObsYearsMysore3rds = [2004.667	2006.917    2008.834	2009.26  2012.26];
CIMysore3rds = [];
ObsMysore3rds = {[0.2611 ],	[0.2424 ], [0.054], [0.111], [0.046]};
ObsVarsMysore3rds = {[7], [7], [8], [7], [8]};
NbSamples = [429 425 425 425 425];
% ObsYearsMysore2rds = [2004.667	2006.917    2008.834	];
% ObsMysore2rdsMin = [0.2193 0.1911 0.032663 ];
% ObsMysore2rdsMax = [0.3028 0.2945 0.07557 ];
% ObsMysore2rds = [0.2611	0.2424 0.054 ];
% ObsVarsMysore2rds = [7 7 8 ];

if ind == 1

  Name = ['HIV_Mysore_3rounds' mode];
  ParametersMysore.RealData = 1;
  ParametersMysore.ObsYears = ObsYearsMysore3rds;
  ParametersMysore.ObsVars = ObsVarsMysore3rds;
  ParametersMysore.Obs = ObsMysore3rds;
  ParametersMysore.NameToSave = Name;
  ParametersMysore.NbSamples = NbSamples;
  
  if IndDiff == 1
    ParametersMysore.DiffusionType = 'Add';
  elseif IndDiff == 2
    ParametersMysore.DiffusionType = 'Bertallanfy';
  elseif IndDiff == 3
    ParametersMysore.DiffusionType = 'Sigmoid';
  end
  ParametersMysore.NbRounds = NbRounds;
  ParametersMysore.TakeClients = TakeClients;
  
  Res = HIVapplyInference([],ParametersMysore);
end

%SavePath = 'S:\Results';
%Res = load([SavePath Name ]);

% if ind == 2
%   Name = 'HIV_Mysore_2rounds';
%   ParametersMysore.RealData = 1;
%   ParametersMysore.ObsYears = ObsYearsMysore2rds;
%   ParametersMysore.ObsVars = ObsVarsMysore2rds;
%   ParametersMysore.Obs = ObsMysore2rds;
%   ParametersMysore.ObsMin = ObsMysore2rdsMin;
%   ParametersMysore.ObsMax = ObsMysore2rdsMax;
%   ParametersMysore.NameToSave = Name;
%   
%   
%   if strcmp(Diff,'Logistic')
%     ParametersMysore.DiffusionType = 'Logistic';
%   end
%   
%   Res = HIVapplyInference([],ParametersMysore);
% end

 
ObsYearsBelgaum3rds = [2005.834	2007.824    2008.584	2010.71 2012.00];
ObsBelgaum3rds = {[0.339 ],	[0.062], [0.273], [0.223], [0.03]};
ObsVarsBelgaum3rds = {[7], [8], [7], [7], [8]};
NbSamples = [363 408 396 423 407];

% ObsYearsBelgaum2rds = [2005.834	2007.824    2008.584 ];
% ObsBelgaum2rdsMin = [0.2762  0.0363 0.2217];
% ObsBelgaum2rdsMax = [0.4018  0.0877 0.3251];
% ObsBelgaum2rds = [0.339	0.062 0.273 ];
% ObsVarsBelgaum2rds = [7 8 7];


if ind == 2
  Name = ['HIV_Belgaum_3rounds' mode];
  ParametersBelgaum.RealData = 1;
  ParametersBelgaum.ObsYears = ObsYearsBelgaum3rds;
  ParametersBelgaum.ObsVars = ObsVarsBelgaum3rds;
  ParametersBelgaum.Obs = ObsBelgaum3rds;
  ParametersBelgaum.NameToSave = Name;
  ParametersBelgaum.NbSamples = NbSamples;

  
  
  if IndDiff == 1
    ParametersBelgaum.DiffusionType = 'Add';
  elseif IndDiff == 2
    ParametersBelgaum.DiffusionType = 'Bertallanfy';
  elseif IndDiff == 3
    ParametersBelgaum.DiffusionType = 'Sigmoid';
  end
  ParametersBelgaum.NbRounds = NbRounds;
  ParametersBelgaum.TakeClients = TakeClients;

    
  Res = HIVapplyInference([],ParametersBelgaum);
end


% if ind == 4
%   Name = 'HIV_Belgaum_2rounds';
%   ParametersBelgaum.RealData = 1;
%   ParametersBelgaum.ObsYears = ObsYearsBelgaum2rds;
%   ParametersBelgaum.ObsVars = ObsVarsBelgaum2rds;
%   ParametersBelgaum.Obs = ObsBelgaum2rds;
%   ParametersBelgaum.ObsMin = ObsBelgaum2rdsMin;
%   ParametersBelgaum.ObsMax = ObsBelgaum2rdsMax;
%   ParametersBelgaum.NameToSave = Name;
%   
%   
%   
%   if strcmp(Diff,'Logistic')
%     ParametersBelgaum.DiffusionType = 'Logistic';
%   end
%   
%   Res = HIVapplyInference([],ParametersBelgaum);
% end


ObsYearsBellary = [2005.908	2007.873  2008.642	2010.865 2012.0];
ObsBellary = {[0.156]	[0.060] [0.142] [0.0634] [0.063]};
ObsVarsBellary = {[7] [8] [7] [7] [8]};
NbSamples = [422 424 410 400 398];


if ind == 3
  Name = ['HIV_Bellary' mode];
  ParametersBellary.RealData = 1;
  ParametersBellary.ObsYears = ObsYearsBellary;
  ParametersBellary.ObsVars = ObsVarsBellary;
  ParametersBellary.Obs = ObsBellary;
  ParametersBellary.NameToSave = Name;
  ParametersBellary.NbSamples = NbSamples;

  
  if IndDiff == 1
    ParametersBellary.DiffusionType = 'Add';
  elseif IndDiff == 2
    ParametersBellary.DiffusionType = 'Bertallanfy';
  elseif IndDiff == 3
    ParametersBellary.DiffusionType = 'Sigmoid';
  end
  ParametersBellary.NbRounds = NbRounds;
  ParametersBellary.TakeClients = TakeClients;

  
  Res = HIVapplyInference([],ParametersBellary);
end



ObsYearsEastGodavry = [2006.25  2006.84 2009.25 2009.35];
ObsEastGodavry = {[0.263],	[0.083], [0.233], [0.096]};
ObsVarsEastGodavry = {[7], [8], [7], [8]};
NbSamples = [422 422 422 422];

if ind == 4
  Name = ['HIV_EastGodavry' mode];
  ParametersEastGodavry.RealData = 1;
  ParametersEastGodavry.ObsYears = ObsYearsEastGodavry;
  ParametersEastGodavry.ObsVars = ObsVarsEastGodavry;
  ParametersEastGodavry.Obs = ObsEastGodavry;
  ParametersEastGodavry.NameToSave = Name;
  ParametersEastGodavry.NbSamples = NbSamples;

  
  if IndDiff == 1
    ParametersEastGodavry.DiffusionType = 'Add';
  elseif IndDiff == 2
    ParametersEastGodavry.DiffusionType = 'Bertallanfy';
  elseif IndDiff == 3
    ParametersEastGodavry.DiffusionType = 'Sigmoid';
  end
  ParametersEastGodavry.NbRounds = NbRounds;
  ParametersEastGodavry.TakeClients = TakeClients;

  
  Res = HIVapplyInference([],ParametersEastGodavry);
end



ObsYearsGuntur = [2006.38	2006.905   2009.53	2009.59];
ObsGuntur = {[0.213],	[0.066], [0.0839], [0.071]};
ObsVarsGuntur = {[7], [8], [7], [8]};
NbSamples = [405 405 405 405];


if ind == 5
  Name = ['HIV_Guntur' mode];
  ParametersGuntur.RealData = 1;
  ParametersGuntur.ObsYears = ObsYearsGuntur;
  ParametersGuntur.ObsVars = ObsVarsGuntur;
  ParametersGuntur.Obs = ObsGuntur;
  ParametersGuntur.NameToSave = Name;
  ParametersGuntur.NbSamples = NbSamples;

  if IndDiff == 1
    ParametersGuntur.DiffusionType = 'Add';
  elseif IndDiff == 2
    ParametersGuntur.DiffusionType = 'Bertallanfy';
  elseif IndDiff == 3
    ParametersGuntur.DiffusionType = 'Sigmoid';
  end
  ParametersGuntur.NbRounds = NbRounds;
  ParametersGuntur.TakeClients = TakeClients;

  
  Res = HIVapplyInference([],ParametersGuntur);
end




ObsYearsHyderabad = [2006.16	2006.96   2009.47	2009.52];
ObsHyderabad = {[0.143],	[0.024], [0.096], [0.037]};
ObsVarsHyderabad = {[7], [8], [7], [8]};
NbSamples = [399 399 399 399];


if ind == 6
  Name = ['HIV_Hyderabad' mode];
  ParametersHyderabad.RealData = 1;
  ParametersHyderabad.ObsYears = ObsYearsHyderabad;
  ParametersHyderabad.ObsVars = ObsVarsHyderabad;
  ParametersHyderabad.Obs = ObsHyderabad;
  ParametersHyderabad.NameToSave = Name;
  ParametersHyderabad.NbSamples = NbSamples;

  if IndDiff == 1
    ParametersHyderabad.DiffusionType = 'Add';
  elseif IndDiff == 2
    ParametersHyderabad.DiffusionType = 'Bertallanfy';
  elseif IndDiff == 3
    ParametersHyderabad.DiffusionType = 'Sigmoid';
  end
  ParametersHyderabad.NbRounds = NbRounds;
  ParametersHyderabad.TakeClients = TakeClients;

    
  Res = HIVapplyInference([],ParametersHyderabad);
end




% ObsYearsVisag = [2006.37	2006.85  2009.24	2009.34];
% % ObsVisag = [0.142	0.08 0.182 0.05];
% ObsVisagMin = [];
% ObsVisagMax = [];
% ObsVarsVisag = [7 8 7 8];
% 
% Name = '\HIV_Visag';
% HIVFullPMCMC(ParametersVisag,ObsYearsVisag,ObsVarsVisag,ObsVisag,Name)


% ObsYearsWarangal = [ 2006.152 2006.754	2009.58	2009.678];
% ObsWarangal = [0.108	0.067 0.149 0.028];
% ObsVarsWarangal = [7 8 7 8];
% 
% Name = '\HIV_Warangal';
% HIVFullPMCMC(ParametersWarangal,ObsYearsWarangal,ObsVarsWarangal,ObsWarangal,Name)


ObsYearsYevatmal = [2006.337	2006.902    2009.728	2009.9];
ObsYevatmal = {[0.373],	[0.109], [0.267], [0.117]};
ObsVarsYevatmal = {[7], [8], [7], [8]};
NbSamples = [153 153 153 153];


if ind == 7
  Name = ['HIV_Yevatmal' mode];
  ParametersYevatmal.RealData = 1;
  ParametersYevatmal.ObsYears = ObsYearsYevatmal;
  ParametersYevatmal.ObsVars = ObsVarsYevatmal;
  ParametersYevatmal.Obs = ObsYevatmal;
  ParametersYevatmal.NameToSave = Name;
  ParametersYevatmal.NbSamples = NbSamples;

  if IndDiff == 1
    ParametersYevatmal.DiffusionType = 'Add';
  elseif IndDiff == 2
    ParametersYevatmal.DiffusionType = 'Bertallanfy';
  elseif IndDiff == 3
    ParametersYevatmal.DiffusionType = 'Sigmoid';
  end
  ParametersYevatmal.NbRounds = NbRounds;
  ParametersYevatmal.TakeClients = TakeClients;

  
  Res = HIVapplyInference([],ParametersYevatmal);
end



ObsYearsShimoga = [2005.688	2007.943  2008.718 2011.5833 2011.9167 ];
ObsShimoga = {[0.0968],	 [0.023], [0.0896 ], [0.07], [0.019]};
ObsVarsShimoga = {[7], [8], [7], [7], [8]};
NbSamples = [389 426 406 396 393 377];


if ind == 8
  Name = ['HIV_Shimoga' mode];
  ParametersShimoga.RealData = 1;
  ParametersShimoga.ObsYears = ObsYearsShimoga;
  ParametersShimoga.ObsVars = ObsVarsShimoga;
  ParametersShimoga.Obs = ObsShimoga;
  ParametersShimoga.NameToSave = Name;
  ParametersShimoga.NbSamples = NbSamples;

  if IndDiff == 1
    ParametersShimoga.DiffusionType = 'Add';
  elseif IndDiff == 2
    ParametersShimoga.DiffusionType = 'Bertallanfy';
  elseif IndDiff == 3
    ParametersShimoga.DiffusionType = 'Sigmoid';
  end
  ParametersShimoga.NbRounds = NbRounds;
  ParametersShimoga.TakeClients = TakeClients;

  
  Res = HIVapplyInference([],ParametersShimoga);
end





ObsYearsChennai = [2006.588	2006.9 2009.5];
ObsChennai = {[0.0317], [0.022], [0.02]};
ObsVarsChennai = {[7], [8], [7]};
NbSamples = [410 405 397];

if ind == 9
  Name = 'HIV_Chennai';
  ParametersChennai.RealData = 1;
  ParametersChennai.ObsYears = ObsYearsChennai;
  ParametersChennai.ObsVars = ObsVarsChennai;
  ParametersChennai.Obs = ObsChennai;
  ParametersChennai.NameToSave = Name;
  ParametersChennai.NbSamples = NbSamples;

  if IndDiff == 1
    ParametersChennai.DiffusionType = 'Add';
  elseif IndDiff == 2
    ParametersChennai.DiffusionType = 'Bertallanfy';
  elseif IndDiff == 3
    ParametersChennai.DiffusionType = 'Sigmoid';
  end
  ParametersChennai.NbRounds = NbRounds;
  ParametersChennai.TakeClients = TakeClients;

  
  Res = HIVapplyInference([],ParametersChennai);
end

ObsYearsMadurai = [2006.3 2006.9	2009.3  2009.67];
ObsMadurai = {[0.057], [0.0225], [0.088], [0.037 ]};
ObsVarsMadurai = {[7], [8], [7], [8]};
NbSamples = [402 400 396 401];

if ind == 10
  Name = 'HIV_Madurai';
  ParametersMadurai.RealData = 1;
  ParametersMadurai.ObsYears = ObsYearsMadurai;
  ParametersMadurai.ObsVars = ObsVarsMadurai;
  ParametersMadurai.Obs = ObsMadurai;
  ParametersMadurai.NameToSave = Name;
  ParametersMadurai.NbSamples = NbSamples;

  if IndDiff == 1
    ParametersMadurai.DiffusionType = 'Add';
  elseif IndDiff == 2
    ParametersMadurai.DiffusionType = 'Bertallanfy';
  elseif IndDiff == 3
    ParametersMadurai.DiffusionType = 'Sigmoid';
  end
  ParametersMadurai.NbRounds = NbRounds;
  ParametersMadurai.TakeClients = TakeClients;

  
  Res = HIVapplyInference([],ParametersMadurai);
end

ObsYearsSalem = [2006.3	2006.8 2009.3  2009.47];
ObsSalem = {[0.1294], [0.035], [0.113], [0.019]};
ObsVarsSalem = {[7], [8], [7], [8]};
NbSamples = [392 396 407 407];

if ind == 11
  Name = 'HIV_Salem';
  ParametersSalem.RealData = 1;
  ParametersSalem.ObsYears = ObsYearsSalem;
  ParametersSalem.ObsVars = ObsVarsSalem;
  ParametersSalem.Obs = ObsSalem;
  ParametersSalem.NameToSave = Name;
  ParametersSalem.NbSamples = NbSamples;

  if IndDiff == 1
    ParametersSalem.DiffusionType = 'Add';
  elseif IndDiff == 2
    ParametersSalem.DiffusionType = 'Bertallanfy';
  elseif IndDiff == 3
    ParametersSalem.DiffusionType = 'Sigmoid';
  end
  ParametersSalem.NbRounds = NbRounds;
  ParametersSalem.TakeClients = TakeClients;

  
  Res = HIVapplyInference([],ParametersSalem);
end



