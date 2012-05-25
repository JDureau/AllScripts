
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Dengue';
% SavePath = 'S:\Results';
% 
% load([SavePath '/ForCovForHIVTransfTetas.mat'])
% Cov = cov(ResRW.TransfThetas');


Parameters = struct();

% Initializing all model parameters
% Parameters.CovInit = Cov;
Parameters.SInit.Value = 0.25;
Parameters.SInit.Min = -10^14;
Parameters.SInit.Max =  10^14;
Parameters.SInit.MinLim = 0.1;
Parameters.SInit.MaxLim = 0.4;
Parameters.SInit.Estimated = 1;
Parameters.SInit.TransfType = 'Logit';
Parameters.R1Init.Value = 0.11;
Parameters.R1Init.Min = -10^14;
Parameters.R1Init.Max =  10^14;
Parameters.R1Init.MinLim = 0.02;
Parameters.R1Init.MaxLim = 0.3;
Parameters.R1Init.Estimated = 1;
Parameters.R1Init.TransfType = 'Logit';
Parameters.R2Init.Value = 0.13;
Parameters.R2Init.Min = -10^14;
Parameters.R2Init.Max =  10^14;
Parameters.R2Init.MinLim = 0.02;
Parameters.R2Init.MaxLim = 0.3;
Parameters.R2Init.Estimated = 1;
Parameters.R2Init.TransfType = 'Logit';
Parameters.I1Init.Value = 10^-6;
Parameters.I1Init.Min = -10^14;
Parameters.I1Init.Max =  10^14;
Parameters.I1Init.MinLim = 10^-9;
Parameters.I1Init.MaxLim = 10^-4;
Parameters.I1Init.Estimated = 0;
Parameters.I1Init.TransfType = 'Logit';
Parameters.I2Init.Value = 10^-6;
Parameters.I2Init.Min = -10^14;
Parameters.I2Init.Max =  10^14;
Parameters.I2Init.MinLim = 10^-9;
Parameters.I2Init.MaxLim = 10^-4;
Parameters.I2Init.Estimated = 0;
Parameters.I2Init.TransfType = 'Logit';
Parameters.I12Init.Value = 10^-6;
Parameters.I12Init.Min = -10^14;
Parameters.I12Init.Max =  10^14;
Parameters.I12Init.MinLim = 10^-9;
Parameters.I12Init.MaxLim = 10^-4;
Parameters.I12Init.Estimated = 0;
Parameters.I12Init.TransfType = 'Logit';
Parameters.I21Init.Value = 10^-6;
Parameters.I21Init.Min = -10^14;
Parameters.I21Init.Max =  10^14;
Parameters.I21Init.MinLim = 10^-9;
Parameters.I21Init.MaxLim = 10^-4;
Parameters.I21Init.Estimated = 0;
Parameters.I21Init.TransfType = 'Logit';
Parameters.Q1Init.Value = 0.0005;
Parameters.Q1Init.Min = -10^14;
Parameters.Q1Init.Max =  10^14;
Parameters.Q1Init.MinLim = 10^-5;
Parameters.Q1Init.MaxLim = 10^-3;
Parameters.Q1Init.Estimated = 0;
Parameters.Q1Init.TransfType = 'Logit';
Parameters.Q2Init.Value = 0.0008;
Parameters.Q2Init.Min = -10^14;
Parameters.Q2Init.Max =  10^14;
Parameters.Q2Init.MinLim = 10^-5;
Parameters.Q2Init.MaxLim = 10^-3;
Parameters.Q2Init.Estimated = 0;
Parameters.Q2Init.TransfType = 'Logit';
Parameters.r0Init.Value = 2.5;
Parameters.r0Init.Min = -10^14;
Parameters.r0Init.Max =  10^14;
Parameters.r0Init.MinLim = 1.1;
Parameters.r0Init.MaxLim = 4;
Parameters.r0Init.Estimated = 1;
Parameters.r0Init.TransfType = 'Logit';
Parameters.vm1.Value = 3.7;
Parameters.vm1.Min = -10^14;
Parameters.vm1.Max =  10^14;
Parameters.vm1.MinLim = 2;
Parameters.vm1.MaxLim = 12;
Parameters.vm1.Estimated = 1;
Parameters.vm1.TransfType = 'Logit';
Parameters.qm1.Value = 6.1;
Parameters.qm1.Min = -10^14;
Parameters.qm1.Max =  10^14;
Parameters.qm1.MinLim = 3;
Parameters.qm1.MaxLim = 20;
Parameters.qm1.Estimated = 1;
Parameters.qm1.TransfType = 'Logit';
Parameters.eta.Value = 0.25;
Parameters.eta.Min = -10^14;
Parameters.eta.Max =  10^14;
Parameters.eta.MinLim = 0.1;
Parameters.eta.MaxLim = 10;
Parameters.eta.Estimated = 1;
Parameters.eta.TransfType = 'Logit';
Parameters.e.Value = 0.17;
Parameters.e.Min = -10^14;
Parameters.e.Max =  10^14;
Parameters.e.MinLim = 0.07;
Parameters.e.MaxLim = 0.3;
Parameters.e.Estimated = 1;
Parameters.e.TransfType = 'Logit';
Parameters.d.Value = 0.75;
Parameters.d.Min = -10^14;
Parameters.d.Max =  10^14;
Parameters.d.MinLim = 0.2;
Parameters.d.MaxLim = 0.9;
Parameters.d.Estimated = 1;
Parameters.d.TransfType = 'Logit';
Parameters.ade.Value = 0.5;
Parameters.ade.Min = -10^14;
Parameters.ade.Max =  10^14;
Parameters.ade.MinLim = 0.2;
Parameters.ade.MaxLim = 0.95;
Parameters.ade.Estimated = 1;
Parameters.ade.TransfType = 'Logit';
Parameters.rep2.Value = 0.09;
Parameters.rep2.Min = -10^14;
Parameters.rep2.Max =  10^14;
Parameters.rep2.MinLim = 0.03;
Parameters.rep2.MaxLim = 0.2;
Parameters.rep2.Estimated = 1;
Parameters.rep2.TransfType = 'Logit';
Parameters.phi.Value = 0.3;
Parameters.phi.Min = -10^14;
Parameters.phi.Max =  10^14;
Parameters.phi.MinLim = 0.1;
Parameters.phi.MaxLim = 0.5;
Parameters.phi.Estimated = 1;
Parameters.phi.TransfType = 'Logit';
Parameters.SigmaRW.Value = 0.001;
Parameters.SigmaRW.Min = -10^14;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.MinLim = 0;
Parameters.SigmaRW.MaxLim = 0.4;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Logit';
Parameters.NbVariables = 11;
Parameters.ObservationLength = 277;
Parameters.PopulationSize.Value = 1649457;
Parameters.rep1.Value = 1;
Parameters.MuB.Value = 0.00128205128205128;
Parameters.MuD.Value = 0.00128205128205128;
Parameters.Problem = 'Dengue';
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.StableCUseConstraint = 0;
TellParsValues(Parameters)

% SavePath = 'S:\Results';
save([SavePath '/ParametersDengue.mat'],'Parameters') 
