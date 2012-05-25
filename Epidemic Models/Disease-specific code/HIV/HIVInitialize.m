function Parameters = HIVInitialize(Parameters)

Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;
Parameters.TotM.Value = Parameters.TotalFSW.Value*Parameters.TotMFactor.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
Parameters.BetaFM.Value = 1-(1-Parameters.BetaFMPerAct.Value)^Parameters.NumberActsPerClient.Value;
Parameters.BetaMF.Value = 1-(1-Parameters.BetaMFPerAct.Value)^Parameters.NumberActsPerClient.Value;

Parameters = UpdateParsNoTransfToTransf(Parameters);

InitialState = [];
InitialState(1) = Parameters.InitialSF1;
InitialState(2) = Parameters.InitialHIVF1;
InitialState(3) = Parameters.InitialSF2;
InitialState(4) = Parameters.InitialHIVF2;
InitialState(5) = Parameters.InitialSM;
InitialState(6) = Parameters.InitialHIVM;
InitialState(7) = (InitialState(2)+InitialState(4))/sum(InitialState(1:4))*100;
InitialState(8) = (InitialState(6))/sum(InitialState(5:6))*100;
InitialState(9) = Parameters.InitialFt.Value;
InitialState(10) = Parameters.InitialDeriv;

Parameters.InitialState = InitialState';