function Parameters = HIV_Initialize(Parameters)

Parameters.TotF1.Value = Parameters.TotalFSW.Value/2;
Parameters.TotF2.Value = Parameters.TotalFSW.Value/2;
Parameters.InitialSF1 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF1.Value;
Parameters.InitialHIVF1 = Parameters.InitialIPropF.Value*Parameters.TotF1.Value;
Parameters.InitialSF2 = (1-Parameters.InitialIPropF.Value)*Parameters.TotF2.Value;
Parameters.InitialHIVF2 = Parameters.InitialIPropF.Value*Parameters.TotF2.Value;



InitialState(1) = Parameters.InitialSF1;
InitialState(2) = Parameters.TotF1.Value - Parameters.InitialSF1;
InitialState(3) = Parameters.InitialSF2;
InitialState(4) = Parameters.TotF2.Value - Parameters.InitialSF2;
Parameters.TotM.Value = Parameters.TotalFSW.Value*(Parameters.CF1.Value+Parameters.CF2.Value)/2/Parameters.CM.Value;
Parameters.InitialSM = (1-Parameters.InitialIPropM.Value)*Parameters.TotM.Value;
Parameters.InitialHIVM = Parameters.InitialIPropM.Value*Parameters.TotM.Value;
% Parameters.NbContactsForMen = (Parameters.CF1.Value + Parameters.CF2.Value)/2*Parameters.TotalFSW.Value/Parameters.TotM.Value;
InitialState(5) = Parameters.InitialSM;
InitialState(6) = Parameters.TotM.Value - Parameters.InitialSM;
InitialState(7) = (InitialState(2) + InitialState(4))/Parameters.TotalFSW.Value*100;
InitialState(8) = (InitialState(6))/Parameters.TotM.Value*100;



try
    if  or(strcmp(Parameters.DiffusionType,'Bertallanfy'),strcmp(Parameters.DiffusionType,'BertallanfyConstr'))
        m = max(eps,Parameters.BRmm1.Value) + 1;
        mu = Parameters.BRmu.Value;
        tinfl = Parameters.BRtinfl.Value;
%         Parameters.BRbase.MaxLim = 0.99*m^(1/(1-m))*mu;
%         Parameters.BRbase.MinLim = 0;
%         Parameters = UpdateParsTransfToNoTransf(Parameters);
        baseline = Parameters.BRbase.Value;%,Parameters.BRbase.Value^(1/100)*0.99*m^(1/(1-m))*mu);
%         baseline = Parameters.BRbase.Value;
        InitialState(9) = (baseline^(1-m)-mu^(1-m))/(1-m);
        B  = (1 - (baseline/mu)^(1-m));
        Parameters.k = 1/tinfl*log(B/(1-m)); 
        
    elseif strcmp(Parameters.DiffusionType,'Sigmoid')
        InitialState(9) = exp(Parameters.Sigmtinfl.Value/Parameters.Sigmrate.Value);
    else
        InitialState(9) = log(Parameters.InitialFt.Value/(1-Parameters.InitialFt.Value));
    end
catch
    InitialState(9) = log(Parameters.InitialFt.Value/(1-Parameters.InitialFt.Value));

end
% InitialState(10) = 0;
Parameters.BetaFM.Value = 1-(1-Parameters.BetaFMPerAct.Value)^Parameters.NumberActsPerClient.Value;
Parameters.BetaMF.Value = 1-(1-Parameters.BetaMFPerAct.Value)^Parameters.NumberActsPerClient.Value;

Parameters.InitialState = InitialState';