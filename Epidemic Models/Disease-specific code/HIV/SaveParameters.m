%% Saving priors

%% Mysore

SavePath = '/Users/dureaujoseph/Documents/PhD_Data/Avahan/';
% SavePath = '/users/ecologie/dureau/src/AllData/Avahan';
% % SavePath = 'S:\Results';
% % 
% % load([SavePath '/ForCovForHIVTransfTetas.mat'])
% % Cov = cov(ResRW.TransfThetas');
% 
% cd('/users/ecologie/dureau/src/AllScripts')
% addpath([pwd '/General Tools'])
% addpath([pwd '/Toolboxes'])
% addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
% addpath([pwd '/Epidemic Models/Disease-specific code'])
% addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])




Parameters = struct();

% Initializing all model parameters
% Parameters.CovInit = Cov;
Parameters.TotalFSW.Value = 1943;%2144;
Parameters.TotalFSW.Min = 104;
Parameters.TotalFSW.Max = 30752;
Parameters.TotalFSW.Estimates = 0;
Parameters.TotalFSW.TransfType = 'Log';
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.DerInitialFt.TransfType = 'Logit';
Parameters.DerInitialFt.Init = 1;
Parameters.InitialFt.Value = 0.01;
Parameters.InitialFt.Min = -10^14;
Parameters.InitialFt.Max = 10^14;
Parameters.InitialFt.MaxLim = 0.95;
Parameters.InitialFt.MinLim = 0.0001;
Parameters.InitialFt.Estimated = 1;
Parameters.InitialFt.TransfType = 'Logit';
Parameters.InitialFt.Init = 1;
Parameters.Rho.Value = 1.2;
Parameters.Rho.Min = -10^14;
Parameters.Rho.Max = 10^14;
Parameters.Rho.MaxLim = 1.25;
Parameters.Rho.MinLim = 1;
Parameters.Rho.Estimated = 1;
Parameters.Rho.TransfType = 'Logit';
Parameters.Rho.Init = 1;
Parameters.InitialIPropF.Value = 0.002;
% Parameters.InitialIPropF.Min = 0;
% Parameters.InitialIPropF.Max = 0.04;
Parameters.InitialIPropF.Min = -10^14;
Parameters.InitialIPropF.Max =  10^14;
Parameters.InitialIPropF.MaxLim = 0.05;
Parameters.InitialIPropF.MinLim = 0;
Parameters.InitialIPropF.Estimated = 1;
Parameters.InitialIPropF.TransfType = 'Logit';
Parameters.InitialIPropF.Init = 1;
Parameters.InitialIPropM.Value = 0.002;
% Parameters.InitialIPropM.Min = 0;
% Parameters.InitialIPropM.Max = 0.02;
Parameters.InitialIPropM.Min = -10^14;
Parameters.InitialIPropM.Max =  10^14;
Parameters.InitialIPropM.MaxLim = 0.05;
Parameters.InitialIPropM.MinLim = 0;
Parameters.InitialIPropM.Estimated = 1;
Parameters.InitialIPropM.TransfType = 'Logit';
Parameters.InitialIPropM.Init = 1;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotMFactor.Value = 10; % arbitrary value, does not have a role anyway;
Parameters.TotMFactor.Min = -10^14;
Parameters.TotMFactor.Max = 10^14;
Parameters.TotMFactor.MinLim = 7;
Parameters.TotMFactor.MaxLim = 19;
Parameters.TotMFactor.Estimated = 0;
Parameters.TotMFactor.TransfType = 'Logit';
Parameters.TotMFactor.Init = 1;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.Alpham1.Value = 93;
Parameters.Alpham1.Min = -10^14;
Parameters.Alpham1.Max = 10^14;
Parameters.Alpham1.MinLim = 87; % was 92
Parameters.Alpham1.MaxLim = 138.5;
Parameters.Alpham1.Estimated = 1;
Parameters.Alpham1.TransfType = 'Logit';
Parameters.Alpham1.PlotInv = 1;
Parameters.MuFm1.Value = 50;
Parameters.MuFm1.Min = -10^14;
Parameters.MuFm1.Max =  10^14;
Parameters.MuFm1.MinLim = 45.4588; 
Parameters.MuFm1.MaxLim = 53.67;
Parameters.MuFm1.Estimated = 1;
Parameters.MuFm1.TransfType = 'Logit';
Parameters.MuFm1.PlotInv = 1;
Parameters.MuMm1.Value = 170;
Parameters.MuMm1.Min = -10^14;
Parameters.MuMm1.Max = 10^14;
Parameters.MuMm1.MinLim = 154.3;
Parameters.MuMm1.MaxLim = 191.32;
Parameters.MuMm1.Estimated = 1;
Parameters.MuMm1.PlotInv = 1;
Parameters.MuMm1.TransfType = 'Logit';
Parameters.BetaMFPerAct.Value = 0.0032;
Parameters.BetaMFPerAct.Min = -10^14;
Parameters.BetaMFPerAct.Max =  10^14;
Parameters.BetaMFPerAct.MinLim = 0.0006*1;
Parameters.BetaMFPerAct.MaxLim = 0.0011*5;
Parameters.BetaMFPerAct.Estimated = 1;
Parameters.BetaMFPerAct.TransfType = 'Logit';
Parameters.BetaFMPerAct.Value = 0.0048;
Parameters.BetaFMPerAct.Min = -10^14;
Parameters.BetaFMPerAct.Max =  10^14;
Parameters.BetaFMPerAct.MinLim = 0.0001*1;
Parameters.BetaFMPerAct.MaxLim = 0.0014*5;
Parameters.BetaFMPerAct.Estimated = 1;
Parameters.BetaFMPerAct.TransfType = 'Logit';
Parameters.NumberActsPerClient.Value = 1.7;
Parameters.NumberActsPerClient.Min = -10^14;
Parameters.NumberActsPerClient.Max =  10^14;
Parameters.NumberActsPerClient.MinLim = 1;
Parameters.NumberActsPerClient.MaxLim = 2;
Parameters.NumberActsPerClient.Estimated = 1;
Parameters.NumberActsPerClient.TransfType = 'Logit';
Parameters.BetaFM.Value = 1-(1-Parameters.BetaFMPerAct.Value)^Parameters.NumberActsPerClient.Value;
Parameters.BetaMF.Value = 1-(1-Parameters.BetaMFPerAct.Value)^Parameters.NumberActsPerClient.Value;
Parameters.eHIV.Value = 0.86;
Parameters.eHIV.Min = -10^14;
Parameters.eHIV.Max =  10^14;
Parameters.eHIV.MaxLim = 0.95;
Parameters.eHIV.MinLim = 0.8;
Parameters.eHIV.Estimated = 1;
Parameters.eHIV.TransfType = 'Logit';
Parameters.CF1.Value = 21;
Parameters.CF1.Min = -10^14;
Parameters.CF1.Max =  10^14;
Parameters.CF1.MinLim = 19.9;
Parameters.CF1.MaxLim = 23.7;
Parameters.CF1.Estimated = 1;
Parameters.CF1.TransfType = 'Logit';
Parameters.CF2.Value = 51;
Parameters.CF2.Min = -10^14;
Parameters.CF2.Max = 10^14;
Parameters.CF2.MinLim = 46.3;
Parameters.CF2.MaxLim =  54.0;
Parameters.CF2.Estimated = 1;
Parameters.CF2.TransfType = 'Logit';
Parameters.CM.Value = 9.47/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 8.64/6;
Parameters.CM.MaxLim =  10.3/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';
Parameters.SigmaRW.Value = 0.1;
Parameters.SigmaRW.Min = -10^14;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.MinLim = 0;
Parameters.SigmaRW.MaxLim = 0.5;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Logit';
Parameters.SigmaOU.Value = 0.1;
Parameters.SigmaOU.Min = -10^14;
Parameters.SigmaOU.Max = 10^14;
Parameters.SigmaOU.MinLim = 0;
Parameters.SigmaOU.MaxLim = 50;
Parameters.SigmaOU.Estimated = 0;
Parameters.SigmaOU.TransfType = 'Log';
Parameters.KappaOU.Value = 0.1;
Parameters.KappaOU.Min = -10^14;
Parameters.KappaOU.Max = 10^14;
Parameters.KappaOU.MinLim = 0;
Parameters.KappaOU.MaxLim = 50;
Parameters.KappaOU.Estimated = 0;
Parameters.KappaOU.TransfType = 'Log';
Parameters.MuOU.Value = 0;
Parameters.MuOU.Min = -10^14;
Parameters.MuOU.Max = 10^14;
Parameters.MuOU.MinLim = 0;
Parameters.MuOU.MaxLim = 1;
Parameters.MuOU.Estimated = 0;
Parameters.MuOU.TransfType = 'Logit';
Parameters.BRmm1.Value = 6;
Parameters.BRmm1.Min = -10^14;
Parameters.BRmm1.Max = 10^14;
Parameters.BRmm1.MinLim = 0;
Parameters.BRmm1.MaxLim = 60;
% Parameters.BRmm1.MaxLim = 10;
Parameters.BRmm1.Estimated = 1;
Parameters.BRmm1.TransfType = 'Logit';
Parameters.BRmm1.Init = 1;
Parameters.BRmu.Value = 0.8;
Parameters.BRmu.Min = -10^14;
Parameters.BRmu.Max = 10^14;
Parameters.BRmu.MaxLim = 0.99;
Parameters.BRmu.MinLim = 0.01;
Parameters.BRmu.Estimated = 1;
Parameters.BRmu.TransfType = 'Logit';
Parameters.BRmu.Init = 1;
Parameters.BRbase.Value = 0.05;
Parameters.BRbase.Min = -10^14;
Parameters.BRbase.Max = 10^14;
Parameters.BRbase.MaxLim = 0.9;
% Parameters.BRbase.MaxLim = 0.1;
Parameters.BRbase.MinLim = 0;
Parameters.BRbase.Estimated = 1;
Parameters.BRbase.TransfType = 'Logit';
Parameters.BRbase.Init = 1;
Parameters.BRtinfl.Value = 100;
Parameters.BRtinfl.Min = -10^14;
Parameters.BRtinfl.Max = 10^14;
Parameters.BRtinfl.MinLim = 1;
Parameters.BRtinfl.MaxLim = 300;
Parameters.BRtinfl.Estimated = 1;
Parameters.BRtinfl.TransfType = 'Logit';
Parameters.BRtinfl.Init = 1;
Parameters.BRsigma.Value = 0.1;
Parameters.BRsigma.Min = -10^14;
Parameters.BRsigma.Max = 10^14;
Parameters.BRsigma.MinLim = 0;
Parameters.BRsigma.MaxLim = 2;
Parameters.BRsigma.Estimated = 1;
Parameters.BRsigma.TransfType = 'Logit';
Parameters.BRsigma.Init = 0;
Parameters.Sigmrate.Value = 6;
Parameters.Sigmrate.Min = -10^14;
Parameters.Sigmrate.Max = 10^14;
Parameters.Sigmrate.MinLim = 0;
Parameters.Sigmrate.MaxLim = 1000;
% Parameters.Sigmmm1.MaxLim = 10;
Parameters.Sigmrate.Estimated = 1;
Parameters.Sigmrate.TransfType = 'Logit';
Parameters.Sigmrate.Init = 1;
Parameters.Sigmmu.Value = 0.8;
Parameters.Sigmmu.Min = -10^14;
Parameters.Sigmmu.Max = 10^14;
Parameters.Sigmmu.MaxLim = 0.99;
Parameters.Sigmmu.MinLim = 0.01;
Parameters.Sigmmu.Estimated = 1;
Parameters.Sigmmu.TransfType = 'Logit';
Parameters.Sigmmu.Init = 1;
Parameters.Sigmbase.Value = 0.05;
Parameters.Sigmbase.Min = -10^14;
Parameters.Sigmbase.Max = 10^14;
Parameters.Sigmbase.MaxLim = 0.9;
% Parameters.Sigmbase.MaxLim = 0.1;
Parameters.Sigmbase.MinLim = 0;
Parameters.Sigmbase.Estimated = 1;
Parameters.Sigmbase.TransfType = 'Logit';
Parameters.Sigmbase.Init = 1;
Parameters.Sigmtinfl.Value = 250;
Parameters.Sigmtinfl.Min = -10^14;
Parameters.Sigmtinfl.Max = 10^14;
Parameters.Sigmtinfl.MinLim = 0;
Parameters.Sigmtinfl.MaxLim = 300;
Parameters.Sigmtinfl.Estimated = 1;
Parameters.Sigmtinfl.TransfType = 'Logit';
Parameters.Sigmtinfl.Init = 1;
Parameters.Sigmsigma.Value = 0.1;
Parameters.Sigmsigma.Min = -10^14;
Parameters.Sigmsigma.Max = 10^14;
Parameters.Sigmsigma.MinLim = 0;
Parameters.Sigmsigma.MaxLim = 2;
Parameters.Sigmsigma.Estimated = 1;
Parameters.Sigmsigma.TransfType = 'Logit';
Parameters.Sigmsigma.Init = 0;
Parameters.InitialDeriv = 0;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.StableCUseConstraint = 0;
ParametersMysore = Parameters;
TellParsValues(Parameters)

% SavePath = 'S:\Results';
save([SavePath '/ParametersMysore.mat'],'Parameters') 

% save([SavePath '/ParametersMysoreBRmInf.mat'],'Parameters') 



%% Belgaum


Parameters = ParametersMysore;

% Initializing all model parameters
% Parameters.CovInit = Cov;
Parameters.TotalFSW.Value = 2000;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.MuFm1.Value = 140;
Parameters.MuFm1.MinLim = 128.35;  %%%% FOR TESTS
Parameters.MuFm1.MaxLim = 157.5175;
Parameters.MuFm1.Estimated = 1;
Parameters.MuMm1.Value = 110;
Parameters.MuMm1.MinLim = 102.7177;
Parameters.MuMm1.MaxLim = 124.1259;
Parameters.CF1.Value = 23;
Parameters.CF1.MinLim = 22.12;
Parameters.CF1.MaxLim = 25.82;
Parameters.CF2.Value = 95;
Parameters.CF2.MinLim = 84.82;
Parameters.CF2.MaxLim =  103.8;
Parameters.CM.Value = 3.67/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 3.3/6;
Parameters.CM.MaxLim =  4.0/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.StableCUseConstraint = 0;
TellParsValues(Parameters)
ParametersBelgaum = Parameters;



% SavePath = 'S:\Results';
save([SavePath '/ParametersBelgaum.mat'],'Parameters') 





%% Bellary

Parameters = ParametersBelgaum;

Parameters.TotalFSW.Value = 3852;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 105;
Parameters.MuFm1.MinLim = 98.7; 
Parameters.MuFm1.MaxLim = 118.35;
Parameters.MuMm1.Value = 110;
Parameters.MuMm1.MinLim = 105.602;
Parameters.MuMm1.MaxLim = 127.77;
Parameters.CF1.Value = 18;
Parameters.CF1.MinLim = 16.52;
Parameters.CF1.MaxLim = 20.08;
Parameters.CF2.Value = 110;
Parameters.CF2.MinLim = 81.18;
Parameters.CF2.MaxLim =  122.08;
Parameters.CM.Value = 4.2/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 3.8/6;
Parameters.CM.MaxLim =  4.6/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';
Parameters.InitialFt.Value = 0.4;



% SavePath = 'S:\Results';
save([SavePath '/ParametersBellary.mat'],'Parameters') 


%% EastGodavry

Parameters = ParametersBelgaum;


Parameters.TotalFSW.Value = 1695;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 100;
Parameters.MuFm1.MinLim = 88.86; 
Parameters.MuFm1.MaxLim = 106.14;
Parameters.MuMm1.Value = 150;
Parameters.MuMm1.MinLim = 138.3;
Parameters.MuMm1.MaxLim = 170.3916;
Parameters.CF1.Value = 23;
Parameters.CF1.MinLim = 21.91;
Parameters.CF1.MaxLim = 24.99;
Parameters.CF2.Value = 85;
Parameters.CF2.MinLim = 77.35;
Parameters.CF2.MaxLim =  91.17;
Parameters.CM.Value = 9.8/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 8.8/6;
Parameters.CM.MaxLim =  10.8/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';



% SavePath = 'S:\Results';
save([SavePath '/ParametersEastGodavry.mat'],'Parameters') 


%% Guntur

Parameters = ParametersBelgaum;


Parameters.TotalFSW.Value = 6055;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 72;
Parameters.MuFm1.MinLim = 67.026; 
Parameters.MuFm1.MaxLim = 79.54;
Parameters.MuMm1.Value = 230;
Parameters.MuMm1.MinLim = 208.57;
Parameters.MuMm1.MaxLim = 269.844;
Parameters.CF1.Value = 35;
Parameters.CF1.MinLim = 33;
Parameters.CF1.MaxLim = 37.1;
Parameters.CF2.Value = 110;
Parameters.CF2.MinLim = 99.9;
Parameters.CF2.MaxLim =  120;
Parameters.CM.Value = 7.4/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 6.7/6;
Parameters.CM.MaxLim =  8.2/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';



% SavePath = 'S:\Results';
save([SavePath '/ParametersGuntur.mat'],'Parameters') 



%% Hyderabad

Parameters = ParametersBelgaum;


Parameters.TotalFSW.Value = 885;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 55;
Parameters.MuFm1.MinLim = 49.907; 
Parameters.MuFm1.MaxLim = 59.23147;
Parameters.MuMm1.Value = 165;
Parameters.MuMm1.MinLim = 147.1255;
Parameters.MuMm1.MaxLim = 181.9588;
Parameters.CF1.Value = 13;
Parameters.CF1.MinLim = 12.4;
Parameters.CF1.MaxLim = 13.8;
Parameters.CF2.Value = 50;
Parameters.CF2.MinLim = 44.47;
Parameters.CF2.MaxLim =  52.94;
Parameters.CM.Value = 8.4/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 7.4/6;
Parameters.CM.MaxLim =  9.4/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';



% SavePath = 'S:\Results';
save([SavePath '/ParametersHyderabad.mat'],'Parameters') 


%% Visag

Parameters = ParametersBelgaum;



Parameters.TotalFSW.Value = 1312;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 70;
Parameters.MuFm1.MinLim = 64.2; 
Parameters.MuFm1.MaxLim = 74.3;
Parameters.MuMm1.Value = 92;
Parameters.MuMm1.MinLim = 87.8;
Parameters.MuMm1.MaxLim = 99.7;
Parameters.CF1.Value = 35;
Parameters.CF1.MinLim = 33.14;
Parameters.CF1.MaxLim = 37.1;
Parameters.CF2.Value = 92;
Parameters.CF2.MinLim = 88.3;
Parameters.CF2.MaxLim =  96.7;
Parameters.CM.Value = 5.0/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 4.7/6;
Parameters.CM.MaxLim =  5.3/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';



% SavePath = 'S:\Results';
save([SavePath '/ParametersVisag.mat'],'Parameters') 


%% Warangal

Parameters = ParametersBelgaum;

Parameters.TotalFSW.Value = 4042;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 9*12;
Parameters.MuFm1.MinLim = 8.4*12; 
Parameters.MuFm1.MaxLim = 9.9*12;
Parameters.MuMm1.Value = 13*12;
Parameters.MuMm1.MinLim = 12.5*12;
Parameters.MuMm1.MaxLim = 15.2*12;
Parameters.CF1.Value = 21;
Parameters.CF1.MinLim = 20.56;
Parameters.CF1.MaxLim = 22.83;
Parameters.CF2.Value = 75;
Parameters.CF2.MinLim = 64.15;
Parameters.CF2.MaxLim =  90.15;
Parameters.CM.Value = 10.4/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 9.4/6;
Parameters.CM.MaxLim =  11.6/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';



% SavePath = 'S:\Results';
save([SavePath '/ParametersWarangal.mat'],'Parameters') 


%% Yevatmal

Parameters = ParametersBelgaum;



Parameters.TotalFSW.Value = 969;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 45;
Parameters.MuFm1.MinLim = 41.24; 
Parameters.MuFm1.MaxLim = 54.40;
Parameters.MuMm1.Value = 100;
Parameters.MuMm1.MinLim = 91.47;
Parameters.MuMm1.MaxLim = 110.116;
Parameters.CF1.Value = 45;
Parameters.CF1.MinLim = 39.84;
Parameters.CF1.MaxLim = 51.37;
Parameters.CF2.Value = 150;
Parameters.CF2.MinLim = 125.73;
Parameters.CF2.MaxLim =  217.88;
Parameters.CM.Value = 4.8/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 4.4/6;
Parameters.CM.MaxLim =  5.2/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';
% 



% SavePath = 'S:\Results';
save([SavePath '/ParametersYevatmal.mat'],'Parameters') 


%% Shimoga

Parameters = ParametersBelgaum;



Parameters.TotalFSW.Value = 1588;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 60;
Parameters.MuFm1.MinLim = 57.13; 
Parameters.MuFm1.MaxLim = 67.89;
Parameters.MuMm1.Value = 140;
Parameters.MuMm1.MinLim = 127.69;
Parameters.MuMm1.MaxLim = 155.9;
Parameters.CF1.Value = 9.7;
Parameters.CF1.MinLim = 9.3;
Parameters.CF1.MaxLim = 10.91;
Parameters.CF2.Value = 50;
Parameters.CF2.MinLim = 44.8;
Parameters.CF2.MaxLim =  54.01;
Parameters.CM.Value = 2.9/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 2.7/6;
Parameters.CM.MaxLim =  3.2/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';



% SavePath = 'S:\Results';
save([SavePath '/ParametersShimoga.mat'],'Parameters') 

%% Bangalore

Parameters = ParametersBelgaum;



Parameters.TotalFSW.Value = 1588;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 60;
Parameters.MuFm1.MinLim = 57.13; 
Parameters.MuFm1.MaxLim = 67.89;
Parameters.MuMm1.Value = 140;
Parameters.MuMm1.MinLim = 127.69;
Parameters.MuMm1.MaxLim = 155.9;
Parameters.CF1.Value = 9.7;
Parameters.CF1.MinLim = 9.3;
Parameters.CF1.MaxLim = 10.91;
Parameters.CF2.Value = 50;
Parameters.CF2.MinLim = 44.8;
Parameters.CF2.MaxLim =  54.01;
Parameters.CM.Value = 9.47/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 8.64/6;
Parameters.CM.MaxLim =  10.3/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';




% SavePath = 'S:\Results';
save([SavePath '/ParametersBangalore.mat'],'Parameters') 


%% Chennai

Parameters = ParametersBelgaum;



Parameters.TotalFSW.Value = 1588;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 5.4*12;
Parameters.MuFm1.MinLim = 5.1*12; 
Parameters.MuFm1.MaxLim = 5.9*12;
Parameters.MuMm1.Value = 15*12;
Parameters.MuMm1.MinLim = 14.6*12;
Parameters.MuMm1.MaxLim = 18.1*12;
Parameters.CF1.Value  = 4*4.3;
Parameters.CF1.MinLim = 0*4.3;
Parameters.CF1.MaxLim = 6*4.3;
Parameters.CF2.Value = 10*4.3;
Parameters.CF2.MinLim =   7*4.3;
Parameters.CF2.MaxLim =  20*4.3;
Parameters.CM.Value = 12/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 11.7/6;
Parameters.CM.MaxLim =  13.1/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';




% SavePath = 'S:\Results';
save([SavePath '/ParametersChennai.mat'],'Parameters') 

%% Madurai

Parameters = ParametersBelgaum;



Parameters.TotalFSW.Value = 1588;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 7*12;
Parameters.MuFm1.MinLim = 6.72*12; 
Parameters.MuFm1.MaxLim = 7.8*12;
Parameters.MuMm1.Value = 8.5*12;
Parameters.MuMm1.MinLim = 7.8*12;
Parameters.MuMm1.MaxLim = 9.2*12;
Parameters.CF1.Value = 4*4.3;
Parameters.CF1.MinLim = 3.38*4.3;
Parameters.CF1.MaxLim = 4.12*4.3;
Parameters.CF2.Value = 11*4.3;
Parameters.CF2.MinLim = 10.8*4.3;
Parameters.CF2.MaxLim =  13.2*4.3;
Parameters.CM.Value = 4/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 3.7/6;
Parameters.CM.MaxLim =  4.5/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';




% SavePath = 'S:\Results';
save([SavePath '/ParametersMadurai.mat'],'Parameters') 

%% Mumbai

Parameters = ParametersBelgaum;



Parameters.TotalFSW.Value = 1588;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 60;
Parameters.MuFm1.MinLim = 57.13; 
Parameters.MuFm1.MaxLim = 67.89;
Parameters.MuMm1.Value = 140;
Parameters.MuMm1.MinLim = 127.69;
Parameters.MuMm1.MaxLim = 155.9;
Parameters.CF1.Value = 9.7;
Parameters.CF1.MinLim = 9.3;
Parameters.CF1.MaxLim = 10.91;
Parameters.CF2.Value = 50;
Parameters.CF2.MinLim = 44.8;
Parameters.CF2.MaxLim =  54.01;
Parameters.CM.Value = 9.47/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 8.64/6;
Parameters.CM.MaxLim =  10.3/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';




% SavePath = 'S:\Results';
save([SavePath '/ParametersMumbai.mat'],'Parameters') 


%% Pune

Parameters = ParametersBelgaum;



Parameters.TotalFSW.Value = 1588;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 60;
Parameters.MuFm1.MinLim = 57.13; 
Parameters.MuFm1.MaxLim = 67.89;
Parameters.MuMm1.Value = 140;
Parameters.MuMm1.MinLim = 127.69;
Parameters.MuMm1.MaxLim = 155.9;
Parameters.CF1.Value = 9.7;
Parameters.CF1.MinLim = 9.3;
Parameters.CF1.MaxLim = 10.91;
Parameters.CF2.Value = 50;
Parameters.CF2.MinLim = 44.8;
Parameters.CF2.MaxLim =  54.01;
Parameters.CM.Value = 9.47/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 8.64/6;
Parameters.CM.MaxLim =  10.3/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';




% SavePath = 'S:\Results';
save([SavePath '/ParametersPune.mat'],'Parameters') 

%% Salem

Parameters = ParametersBelgaum;



Parameters.TotalFSW.Value = 1588;
Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.MuFm1.Value = 5*12;
Parameters.MuFm1.MinLim = 4.8*12; 
Parameters.MuFm1.MaxLim = 5.6*12;
Parameters.MuMm1.Value = 12*12;
Parameters.MuMm1.MinLim = 11.5*12;
Parameters.MuMm1.MaxLim = 13.9*12;
Parameters.CF1.Value = 4*4.3;
Parameters.CF1.MinLim = 3.5*4.3;
Parameters.CF1.MaxLim = 4.07*4.3;
Parameters.CF2.Value = 10.2*4.3;
Parameters.CF2.MinLim = 9.7*4.3;
Parameters.CF2.MaxLim =  10.9*4.3;
Parameters.CM.Value = 9.47/6;
Parameters.CM.Min = -10^14;
Parameters.CM.Max = 10^14;
Parameters.CM.MinLim = 8.9/6;
Parameters.CM.MaxLim =  11/6;
Parameters.CM.Estimated = 1;
Parameters.CM.TransfType = 'Logit';




% SavePath = 'S:\Results';
save([SavePath '/ParametersSalem.mat'],'Parameters') 


