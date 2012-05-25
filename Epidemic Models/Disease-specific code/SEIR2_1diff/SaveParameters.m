% Saving priors SEIR


% 
% Parameters = struct();
% 
% Parameters.Message = 'For every data-set, calibrate accordingly EInit, RInit etc..';
% Parameters.Problem = 'MarcFlu';
% Parameters.NbVariables = 7;
% Parameters.SigmaObs = 0.1;
% Parameters.DiffusionType = 'Add';
% Parameters.ObservationLength = 7*35;
% Parameters.ComputationTStep = 0.1;
% Parameters.TotalPopulation = 100000;
% Parameters.km1.Value = 1.5647;
% Parameters.km1.Min = -10^14;
% Parameters.km1.Max = 10^14;
% Parameters.km1.MinLim = 1.55;
% Parameters.km1.MaxLim = 1.63;
% Parameters.km1.Estimated = 1;
% Parameters.km1.TransfType = 'Logit';
% Parameters.km1.PlotInv = 1;
% Parameters.km1.Sample = 1;
% Parameters.gammam1.Value = 0.94963;
% Parameters.gammam1.Min = -10^14;
% Parameters.gammam1.Max = 10^14;
% Parameters.gammam1.MinLim = 0.93;
% Parameters.gammam1.MaxLim = 1.23;
% Parameters.gammam1.Estimated = 1;
% Parameters.gammam1.TransfType = 'Logit';
% Parameters.gammam1.PlotInv = 1;
% Parameters.gammam1.Sample = 1;
% Parameters.betainit.Value = 1.5;
% Parameters.betainit.Min = -10^14;
% Parameters.betainit.Max = 10^14;
% Parameters.betainit.Estimated = 1;
% Parameters.betainit.TransfType = 'Log';
% Parameters.betainit.Init = 1;
% Parameters.betainit.Sample = 0;
% Parameters.betaderinit.Value = 0;
% Parameters.betaderinit.Min = -10^14;
% Parameters.betaderinit.Max = 10^14;
% Parameters.betaderinit.MinLim = -100;
% Parameters.betaderinit.MaxLim =  100;
% Parameters.betaderinit.Estimated = 0;
% Parameters.betaderinit.TransfType = 'Logit';
% Parameters.betaderinit.Init = 1;
% Parameters.betaderinit.Sample = 0;
% Parameters.EInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% Parameters.EInitProp.Min = -10^14;%0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% Parameters.EInitProp.Max = 10^14;%10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% Parameters.EInitProp.MinLim = 0;
% Parameters.EInitProp.MaxLim = 1;
% Parameters.EInitProp.Estimated = 1;
% Parameters.EInitProp.TransfType = 'Logit';
% Parameters.EInitProp.Init = 1;
% Parameters.EInitProp.Sample = 0;
% Parameters.IInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% Parameters.IInitProp.Min = -10^14;%0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% Parameters.IInitProp.Max = -10^14;%10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
% Parameters.IInitProp.MinLim = 0;
% Parameters.IInitProp.MaxLim = 1;
% Parameters.IInitProp.Estimated = 1;
% Parameters.IInitProp.TransfType = 'Logit';
% Parameters.IInitProp.Init = 1;
% Parameters.IInitProp.Sample = 0;
% Parameters.RInitProp.Value = 0.2;
% Parameters.RInitProp.Min = 0;
% Parameters.RInitProp.Max = 0.30;
% Parameters.RInitProp.MinLim = 0;
% Parameters.RInitProp.MaxLim = 0.6;
% Parameters.RInitProp.Estimated = 1;
% Parameters.RInitProp.TransfType = 'Logit';
% Parameters.RInitProp.Init = 1;
% Parameters.RInitProp.Sample = 1;
% Parameters.SigmaRW.Value = exp(-0.6);
% Parameters.SigmaRW.Min = -10^14;
% Parameters.SigmaRW.Max = 10^14;
% Parameters.SigmaRW.Estimated = 1;
% Parameters.SigmaRW.Sample = 0;
% Parameters.SigmaRW.TransfType = 'Log';
% Parameters.InitialCovFact.Value = 0.40;
% Parameters.InitialCovFact.Min = -10^14;
% Parameters.InitialCovFact.Max =  10^14;
% Parameters.InitialCovFact.Estimated =  0;
% Parameters.InitialCovFact.TransfType = 'Log';
% Parameters.InitialCovFact.Init = 1;
% Parameters.InitialCovFact.Sample =  0;
% Parameters = DefineEstimatedParametersIndexes(Parameters);
% Parameters = DefineTransfFunctions(Parameters);
% Parameters = DefinePriors(Parameters);
% Parameters = UpdateParsNoTransfToTransf(Parameters);
% TellParsValues(Parameters)
% 
% 
% SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc';
% 
% 
% save([SavePath '/ParametersSEIR.mat'])


%% Parameters



Parameters = struct();

Parameters.Message = 'For every data-set, calibrate accordingly EInit, RInit etc..';
Parameters.Problem = 'MarcFlu';
Parameters.NbVariables = 12;
Parameters.SigmaObs = 0.1;
Parameters.DiffusionType = 'Add';
Parameters.ObservationLength = 7*35;
Parameters.ComputationTStep = 0.1;
Parameters.TotalPopulation1 = 31113;
Parameters.TotalPopulation2 = 68887;
Parameters.km1.Value = 1.5647;
% Parameters.km1.Min = -10^14;
% Parameters.km1.Max = 10^14;
Parameters.km1.Min = 1.55;
Parameters.km1.Max = 1.63;
Parameters.km1.Estimated = 1;
Parameters.km1.TransfType = 'Log';
Parameters.km1.PlotInv = 1;
Parameters.km1.Sample = 1;
Parameters.gammam1.Value = 0.94963;
% Parameters.gammam1.Min = -10^14;
% Parameters.gammam1.Max = 10^14;
Parameters.gammam1.Min = 0.93;
Parameters.gammam1.Max = 1.23;
Parameters.gammam1.Estimated = 1;
Parameters.gammam1.TransfType = 'Log';
Parameters.gammam1.PlotInv = 1;
Parameters.gammam1.Sample = 1;
Parameters.beta11init.Value = 0.5;
Parameters.beta11init.Min = -10^14;
Parameters.beta11init.Max = 10^14;
Parameters.beta11init.Estimated = 1;
Parameters.beta11init.TransfType = 'Log';
Parameters.beta11init.Init = 1;
Parameters.beta11init.Sample = 0;
Parameters.beta12init.Value = 0.5;
Parameters.beta12init.Min = -10^14;
Parameters.beta12init.Max = 10^14;
Parameters.beta12init.Estimated = 1;
Parameters.beta12init.TransfType = 'Log';
Parameters.beta12init.Init = 1;
Parameters.beta12init.Sample = 0;
Parameters.beta22init.Value = 0.5;
Parameters.beta22init.Min = -10^14;
Parameters.beta22init.Max = 10^14;
Parameters.beta22init.Estimated = 1;
Parameters.beta22init.TransfType = 'Log';
Parameters.beta22init.Init = 1;
Parameters.beta22init.Sample = 0;
Parameters.EInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation1;
Parameters.EInitProp.Min = -10^14;%0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Max = 10^14;%10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.MinLim = 0;
Parameters.EInitProp.MaxLim = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.EInitProp.TransfType = 'Logit';
Parameters.EInitProp.Init = 1;
Parameters.EInitProp.Sample = 0;
Parameters.IInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation1;
Parameters.IInitProp.Min = -10^14;%0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Max = 10^14;%10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.MinLim = 0;
Parameters.IInitProp.MaxLim = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.IInitProp.TransfType = 'Logit';
Parameters.IInitProp.Init = 1;
Parameters.IInitProp.Sample = 0;
Parameters.RInitProp.Value = 0.2;
Parameters.RInitProp.Min = 0;
Parameters.RInitProp.Max = 0.30;
Parameters.RInitProp.MinLim = 0;
Parameters.RInitProp.MaxLim = 0.6;
Parameters.RInitProp.Estimated = 1;
Parameters.RInitProp.TransfType = 'Logit';
Parameters.RInitProp.Init = 1;
Parameters.RInitProp.Sample = 1;
Parameters.SigmaRW11.Value = exp(-0.6);
Parameters.SigmaRW11.Min = -10^14;
Parameters.SigmaRW11.Max = 10^14;
Parameters.SigmaRW11.Estimated = 1;
Parameters.SigmaRW11.Sample = 0;
Parameters.SigmaRW11.TransfType = 'Log';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)
% 
% SavePath = 'S:\Results\';
% save([SavePath 'ParametersSEIR.mat'],'Parameters')


SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc';


save([SavePath '/ParametersSEIR2_cst.mat'])


