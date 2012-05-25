function [] =MainHIV(ind,IndDiff)

ind = ind+1;

% Main HIV clean

%% load paths

if IndDiff == 1
    Diffusion = 'Logistic';
end

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

temp = load([SavePath 'ParametersMysore.mat']);
ParametersMysore = temp.Parameters;

temp = load([SavePath 'ParametersBelgaum.mat']);
ParametersBelgaum = temp.Parameters;

temp = load([SavePath 'ParametersBellary.mat']);
ParametersBellary = temp.Parameters;

temp = load([SavePath 'ParametersEastGodavry.mat']);
ParametersEastGodavry = temp.Parameters;

temp = load([SavePath 'ParametersGuntur.mat']);
ParametersGuntur = temp.Parameters;

temp = load([SavePath 'ParametersHyderabad.mat']);
ParametersHyderabad = temp.Parameters;

temp = load([SavePath 'ParametersVisag.mat']);
ParametersVisag = temp.Parameters;

temp = load([SavePath 'ParametersWarangal.mat']);
ParametersWarangal = temp.Parameters;

temp = load([SavePath 'ParametersYevatmal.mat']);
ParametersYevatmal = temp.Parameters;

temp = load([SavePath 'ParametersShimoga.mat']);
ParametersShimoga = temp.Parameters;



ObsYearsMysore3rds = [2004.667	2006.917    2008.834	2009.26];
CIMysore3rds = [];
ObsMysore3rds = [0.2611	0.2424 0.054 0.111];
ObsMysore3rdsMin = [0.2193 0.1911 0.032663 0.06975];
ObsMysore3rdsMax = [0.3028 0.2945 0.07557 0.14818];
ObsVarsMysore3rds = [7 7 8 7];
ObsYearsMysore2rds = [2004.667	2006.917    2008.834	];
ObsMysore2rdsMin = [0.2193 0.1911 0.032663 ];
ObsMysore2rdsMax = [0.3028 0.2945 0.07557 ];
ObsMysore2rds = [0.2611	0.2424 0.054 ];
ObsVarsMysore2rds = [7 7 8 ];

if ind == 1

  Name = 'HIV_Mysore_3rounds';
  ParametersMysore.RealData = 1;
  ParametersMysore.ObsYears = ObsYearsMysore3rds;
  ParametersMysore.ObsVars = ObsVarsMysore3rds;
  ParametersMysore.Obs = ObsMysore3rds;
  ParametersMysore.ObsMin = ObsMysore3rdsMin;
  ParametersMysore.ObsMax = ObsMysore3rdsMax;
  ParametersMysore.NameToSave = Name;
  
  if strcmp(Diff,'Logistic')
    ParametersMysore.DiffusionType = 'Logistic';
  end
  
  Res = HIVapplyInference([],ParametersMysore);
end

%SavePath = 'S:\Results';
%Res = load([SavePath Name ]);

if ind == 2
  Name = 'HIV_Mysore_2rounds';
  ParametersMysore.RealData = 1;
  ParametersMysore.ObsYears = ObsYearsMysore2rds;
  ParametersMysore.ObsVars = ObsVarsMysore2rds;
  ParametersMysore.Obs = ObsMysore2rds;
  ParametersMysore.ObsMin = ObsMysore2rdsMin;
  ParametersMysore.ObsMax = ObsMysore2rdsMax;
  ParametersMysore.NameToSave = Name;
  
  
  if strcmp(Diff,'Logistic')
    ParametersMysore.DiffusionType = 'Logistic';
  end
  
  Res = HIVapplyInference([],ParametersMysore);
end

 
ObsYearsBelgaum3rds = [2005.834	2007.824    2008.584	2010.71];
ObsBelgaum3rds = [0.339	0.062 0.273 0.223];
ObsBelgaum3rdsMin = [0.2762 0.0363 0.2217  0.17552];
ObsBelgaum3rdsMax = [0.4018 0.0877 0.3251  0.2695];
ObsVarsBelgaum3rds = [7 8 7 7];
ObsYearsBelgaum2rds = [2005.834	2007.824    2008.584 ];
ObsBelgaum2rdsMin = [0.2762  0.0363 0.2217];
ObsBelgaum2rdsMax = [0.4018  0.0877 0.3251];
ObsBelgaum2rds = [0.339	0.062 0.273 ];
ObsVarsBelgaum2rds = [7 8 7];


if ind == 3
  Name = 'HIV_Belgaum_3rounds';
  ParametersBelgaum.RealData = 1;
  ParametersBelgaum.ObsYears = ObsYearsBelgaum3rds;
  ParametersBelgaum.ObsVars = ObsVarsBelgaum3rds;
  ParametersBelgaum.Obs = ObsBelgaum3rds;
  ParametersBelgaum.ObsMin = ObsBelgaum3rdsMin;
  ParametersBelgaum.ObsMax = ObsBelgaum3rdsMax;
  ParametersBelgaum.NameToSave = Name;
  
  
  
  if strcmp(Diff,'Logistic')
    ParametersBelgaum.DiffusionType = 'Logistic';
  end
  
  Res = HIVapplyInference([],ParametersBelgaum);
end


if ind == 4
  Name = 'HIV_Belgaum_2rounds';
  ParametersBelgaum.RealData = 1;
  ParametersBelgaum.ObsYears = ObsYearsBelgaum2rds;
  ParametersBelgaum.ObsVars = ObsVarsBelgaum2rds;
  ParametersBelgaum.Obs = ObsBelgaum2rds;
  ParametersBelgaum.ObsMin = ObsBelgaum2rdsMin;
  ParametersBelgaum.ObsMax = ObsBelgaum2rdsMax;
  ParametersBelgaum.NameToSave = Name;
  
  
  
  if strcmp(Diff,'Logistic')
    ParametersBelgaum.DiffusionType = 'Logistic';
  end
  
  Res = HIVapplyInference([],ParametersBelgaum);
end


ObsYearsBellary = [2005.908	2007.873  2008.642	2010.865];
ObsBellary = [0.156	0.060 0.142 0.0634];
ObsBellaryMin = [0.1106 0.0258 0.1048 0.0377];
ObsBellaryMax = [0.2003 0.0946 0.1776 0.0892];
ObsVarsBellary = [7 8 7 7];


if ind == 5
  Name = 'HIV_Bellary';
  ParametersBellary.RealData = 1;
  ParametersBellary.ObsYears = ObsYearsBellary;
  ParametersBellary.ObsVars = ObsVarsBellary;
  ParametersBellary.Obs = ObsBellary;
  ParametersBellary.ObsMin = ObsBellaryMin;
  ParametersBellary.ObsMax = ObsBellaryMax;
  ParametersBellary.NameToSave = Name;
  
  
  if strcmp(Diff,'Logistic')
    ParametersBellary.DiffusionType = 'Logistic';
  end
  
  Res = HIVapplyInference([],ParametersBellary);
end



ObsYearsEastGodavry = [2006.25  2006.84 2009.25 2009.35];
ObsEastGodavry = [0.263	0.083 0.233 0.096];
ObsEastGodavryMin = [0.2004 0.0483 0.1514 0.041215];
ObsEastGodavryMax = [0.3247 0.1184 0.3139 0.150556];
ObsVarsEastGodavry = [7 8 7 8];

if ind == 6
  Name = 'HIV_EastGodavry';
  ParametersEastGodavry.RealData = 1;
  ParametersEastGodavry.ObsYears = ObsYearsEastGodavry;
  ParametersEastGodavry.ObsVars = ObsVarsEastGodavry;
  ParametersEastGodavry.Obs = ObsEastGodavry;
  ParametersEastGodavry.ObsMin = ObsEastGodavryMin;
  ParametersEastGodavry.ObsMax = ObsEastGodavryMax;
  ParametersEastGodavry.NameToSave = Name;
  
  if strcmp(Diff,'Logistic')
    ParametersEastGodavry.DiffusionType = 'Logistic';
  end
  
  Res = HIVapplyInference([],ParametersEastGodavry);
end



ObsYearsGuntur = [2006.38	2006.905   2009.53	2009.59];
ObsGuntur = [0.213	0.066 0.0839 0.071];
ObsGunturMin = [0.1639 0.0369 0.0429 0.021891];
ObsGunturMax = [0.262 0.0956 0.125 0.1206];
ObsVarsGuntur = [7 8 7 8];


if ind == 7
  Name = 'HIV_Guntur';
  ParametersGuntur.RealData = 1;
  ParametersGuntur.ObsYears = ObsYearsGuntur;
  ParametersGuntur.ObsVars = ObsVarsGuntur;
  ParametersGuntur.Obs = ObsGuntur;
  ParametersGuntur.ObsMin = ObsGunturMin;
  ParametersGuntur.ObsMax = ObsGunturMax;
  ParametersGuntur.NameToSave = Name;
  if strcmp(Diff,'Logistic')
    ParametersEastGuntur.DiffusionType = 'Logistic';
  end
  
  Res = HIVapplyInference([],ParametersGuntur);
end




ObsYearsHyderabad = [2006.16	2006.96   2009.47	2009.52];
ObsHyderabad = [0.143	0.024 0.096 0.037];
ObsHyderabadMin = [0.0906 0.007 0.045 0];
ObsHyderabadMax = [0.1954 0.0405 0.147 0.0841];
ObsVarsHyderabad = [7 8 7 8];


if ind == 8
  Name = 'HIV_Hyderabad';
  ParametersHyderabad.RealData = 1;
  ParametersHyderabad.ObsYears = ObsYearsHyderabad;
  ParametersHyderabad.ObsVars = ObsVarsHyderabad;
  ParametersHyderabad.Obs = ObsHyderabad;
  ParametersHyderabad.ObsMin = ObsHyderabadMin;
  ParametersHyderabad.ObsMax = ObsHyderabadMax;
  ParametersHyderabad.NameToSave = Name;
  if strcmp(Diff,'Logistic')
    ParametersHyderabad.DiffusionType = 'Logistic';
  end
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
ObsYevatmal = [0.373	0.109 0.267 0.117];
ObsYevatmalMin = [0.2391 0.075527 0.1805 0.065032];
ObsYevatmalMax = [0.506 0.142355 0.3546 0.168657];
ObsVarsYevatmal = [7 8 7 8];


if ind == 9
  Name = 'HIV_Yevatmal';
  ParametersYevatmal.RealData = 1;
  ParametersYevatmal.ObsYears = ObsYearsYevatmal;
  ParametersYevatmal.ObsVars = ObsVarsYevatmal;
  ParametersYevatmal.Obs = ObsYevatmal;
  ParametersYevatmal.ObsMin = ObsYevatmalMin;
  ParametersYevatmal.ObsMax = ObsYevatmalMax;
  ParametersYevatmal.NameToSave = Name;
  if strcmp(Diff,'Logistic')
    ParametersYevatmal.DiffusionType = 'Logistic';
  end
  Res = HIVapplyInference([],ParametersYevatmal);
end



ObsYearsShimoga = [2005.688	2007.943  2008.718];
ObsShimoga = [0.0968	 0.023 0.0896];
ObsShimogaMin = [0.0632 0.0085 0.0566];
ObsShimogaMax = [0.1305 0.0514 0.1226];
ObsVarsShimoga = [7 8 7];


if ind == 10
  Name = 'HIV_Shimoga';
  ParametersShimoga.RealData = 1;
  ParametersShimoga.ObsYears = ObsYearsShimoga;
  ParametersShimoga.ObsVars = ObsVarsShimoga;
  ParametersShimoga.Obs = ObsShimoga;
  ParametersShimoga.ObsMin = ObsShimogaMin;
  ParametersShimoga.ObsMax = ObsShimogaMax;
  ParametersShimoga.NameToSave = Name;
  if strcmp(Diff,'Logistic')
    ParametersShimoga.DiffusionType = 'Logistic';
  end
  Res = HIVapplyInference([],ParametersShimoga);
end

