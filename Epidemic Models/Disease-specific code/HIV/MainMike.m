%% Mike's estimates



%% load paths

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
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\HIV')

%% Mysore prior


SavePath = 'S:\Results';
temp = load([SavePath '\ParametersMysore.mat']);
ParametersMysore = temp.Parameters;

temp = load([SavePath '\ParametersBelgaum.mat']);
ParametersBelgaum = temp.Parameters;



ObsYears1 = [2004.667	2006.917    2008.834	2009.26];
ObsVars1 = [7 7 8 7];

Obs1 = [0.284579	0.216409    0.038612	0.140494];

Obs2 = [0.298951	0.273208    0.045019	0.24908];


ObsYears2 = [2005.834	2007.824    2008.584	2010.71];
ObsVars2 = [7 8 7 7];

Obs3 = [0.335081	0.038625    0.294871   	0.271126];

Obs4 = [0.381032	0.042048    0.247365   	0.177238];



Name = '\HIV_Mike_Data1_MysorePriors';
HIVFullPMCMC(ParametersMysore,ObsYears1,ObsVars1,Obs1,Name)

Name = '\HIV_Mike_Data1_BelgaumPriors';
HIVFullPMCMC(ParametersBelgaum,ObsYears1,ObsVars1,Obs1,Name)

Name = '\HIV_Mike_Data2_MysorePriors';
HIVFullPMCMC(ParametersMysore,ObsYears1,ObsVars1,Obs2,Name)

Name = '\HIV_Mike_Data2_BelgaumPriors';
HIVFullPMCMC(ParametersBelgaum,ObsYears1,ObsVars1,Obs2,Name)

Name = '\HIV_Mike_Data3_MysorePriors';
HIVFullPMCMC(ParametersMysore,ObsYears2,ObsVars2,Obs3,Name)

Name = '\HIV_Mike_Data3_BelgaumPriors';
HIVFullPMCMC(ParametersBelgaum,ObsYears2,ObsVars2,Obs3,Name)

Name = '\HIV_Mike_Data4_MysorePriors';
HIVFullPMCMC(ParametersMysore,ObsYears2,ObsVars2,Obs4,Name)

Name = '\HIV_Mike_Data4_BelgaumPriors';
HIVFullPMCMC(ParametersBelgaum,ObsYears2,ObsVars2,Obs4,Name)


SavePath = 'S:\Results';
load([SavePath Name '.mat'])
PlotResHIV(ResRW,ResRW.Parameters)

Parameters = ResRW.Parameters;
Names = ResRW.Parameters.Names.Estimated;
k = length(Names);
for i = 1:k
    subplot(ceil(sqrt(k)),ceil(sqrt(k)),i)
    ampl = Parameters.(Names{i}).MaxLim-Parameters.(Names{i}).MinLim;
    xis = Parameters.(Names{i}).MinLim:ampl/1000:Parameters.(Names{i}).MaxLim;
    plot(xis,1/ampl*ones(size(xis)),'g')
    hold on
    [fis,xis] = ksdensity(ResRW.Thetas(i,:));
    plot(xis,fis)
    hold off
    title(Names{i})
end
legend('Prior','Posterior');






