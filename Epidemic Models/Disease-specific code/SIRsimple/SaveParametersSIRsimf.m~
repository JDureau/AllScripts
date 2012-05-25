% Parameters to validate Simforence EKF

Parameters = struct();

Parameters.Problem = 'Simf';
Parameters.NbVariables = 9;
Parameters.DiffusionType = 'Add';
Parameters.p_t = 1000000;
Parameters.r0_init.Value = 20;
Parameters.r0_init.MinLim = 17;
Parameters.r0_init.MaxLim = 25;
Parameters.r0_init.Min = -10^14;
Parameters.r0_init.Max = 10^14;
Parameters.r0_init.Estimated = 1;
Parameters.r0_init.TransfType = 'Logit';
Parameters.vm1.Value = 11;
Parameters.vm1.MinLim = 7;
Parameters.vm1.MaxLim = 16;
Parameters.vm1.Min = -10^14;
Parameters.vm1.Max = 10^14;
Parameters.vm1.Estimated = 1;
Parameters.vm1.TransfType = 'Logit';
Parameters.iota.Value = 2.17215;
Parameters.iota.MinLim = 0.5;
Parameters.iota.MaxLim = 10;
Parameters.iota.Min = -10^14;
Parameters.iota.Max = 10^14;
Parameters.iota.Estimated = 1;
Parameters.iota.TransfType = 'Logit';
Parameters.e.Value = 0.02;
Parameters.e.MinLim = 0.01;
Parameters.e.MaxLim = 0.1;
Parameters.e.Min = -10^14;
Parameters.e.Max = 10^14;
Parameters.e.Estimated = 1;
Parameters.e.TransfType = 'Logit';
Parameters.d.Value = 0.917029;
Parameters.d.MinLim = 0.8;
Parameters.d.MaxLim = 0.99;
Parameters.d.Min = -10^14;
Parameters.d.Max = 10^14;
Parameters.d.Estimated = 1;
Parameters.d.TransfType = 'Logit';
Parameters.vol__r0.Value = 0.6;
Parameters.vol__r0.MinLim = 1e-6;
Parameters.vol__r0.MaxLim = 1;
Parameters.vol__r0.Min = -10^14;
Parameters.vol__r0.Max = 10^14;
Parameters.vol__r0.Estimated = 1;
Parameters.vol__r0.TransfType = 'Logit';
Parameters.Sinit.Value = 0.07;
Parameters.Sinit.MinLim = 0.07;
Parameters.Sinit.MaxLim = 0.07;
Parameters.Sinit.Min = -10^14;
Parameters.Sinit.Max = 10^14;
Parameters.Sinit.Estimated = 0;
Parameters.Sinit.TransfType = 'Logit';
Parameters.Iinit.Value = 1e-05;
Parameters.Iinit.MinLim = 7e-6;
Parameters.Iinit.MaxLim = 1e-4;
Parameters.Iinit.Min = -10^14;
Parameters.Iinit.Max = 10^14;
Parameters.Iinit.Estimated = 1;
Parameters.Iinit.TransfType = 'Logit';
Parameters.rep2.Value = 0.6;
Parameters.rep2.MinLim = 0.5;
Parameters.rep2.MaxLim = 0.7;
Parameters.rep2.Min = -10^14;
Parameters.rep2.Max = 10^14;
Parameters.rep2.Estimated = 1;
Parameters.rep2.TransfType = 'Logit';
Parameters.phi.Value = 0.1;
Parameters.phi.MinLim = 0.07;
Parameters.phi.MaxLim = 0.2;
Parameters.phi.Min = -10^14;
Parameters.phi.Max = 10^14;
Parameters.phi.Estimated = 1;
Parameters.phi.TransfType = 'Logit';
Parameters.mu_b = 0.000273972602739726;
Parameters.mu_d = 0.000273972602739726;
Parameters.rep1 = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)
% 
% SavePath = 'S:\Results\';
% save([SavePath 'ParametersSEIR.mat'],'Parameters')


SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/SimfValidation';

save([SavePath '/Parameters.mat'])


