function res = TrajLogLik(Fts,Parameters,Model,Data)
    
    % creates data from a Beta trajectory

    Parameters = Model.InitializeParameters(Parameters);
    Variables = Parameters.InitialState;


    TempVariables = Variables;

    

    TStep = Parameters.ComputationTStep;
    NbTSteps = Parameters.NbTSteps;
    
    IndBeta = 1;
    Record = [];
    
    LogLik = 0;
    for IndTime = 2:length(Data.Instants)    
        for IndDiscr = 1:Data.NbComputingSteps(IndTime)
            % Variables
            Variables(9) = Fts(IndBeta);
            TempVariables(9) = Fts(IndBeta);
            TotF1 = Parameters.TotF1.Value;
            TotF2 = Parameters.TotF2.Value;
            TotM  = Parameters.TotM.Value;
            TempVariables(1) = TempVariables(1) + ( Parameters.MuFm1.Value^-1*(TotF1) + Parameters.Alpham1.Value^-1*Variables(2) - Parameters.BetaMF.Value*Parameters.CF1.Value*(1-Parameters.eHIV.Value*Variables(9)).*Variables(1).*Variables(6)./TotM-Parameters.MuFm1.Value^-1*Variables(1))*TStep;
            TempVariables(2) = TempVariables(2) + ( Parameters.BetaMF.Value*Parameters.CF1.Value*(1-Parameters.eHIV.Value*Variables(9)).*Variables(1).*Variables(6)./TotM - (Parameters.MuFm1.Value^-1 + Parameters.Alpham1.Value^-1)*Variables(2))*TStep;
            TempVariables(3) = TempVariables(3) + ( Parameters.MuFm1.Value^-1*(TotF2) + Parameters.Alpham1.Value^-1*Variables(4) - Parameters.BetaMF.Value*Parameters.CF2.Value*(1-Parameters.eHIV.Value*Variables(9)).*Variables(3).*Variables(6)./TotM-Parameters.MuFm1.Value^-1*Variables(3))*TStep;
            TempVariables(4) = TempVariables(4) + ( Parameters.BetaMF.Value*Parameters.CF2.Value*(1-Parameters.eHIV.Value*Variables(9)).*Variables(3).*Variables(6)./TotM - (Parameters.MuFm1.Value^-1 + Parameters.Alpham1.Value^-1)*Variables(4))*TStep;
            TempVariables(5) = TempVariables(5) + ( Parameters.MuMm1.Value^-1*(TotM) + Parameters.Alpham1.Value^-1*Variables(6) - Parameters.BetaFM.Value*(1-Parameters.eHIV.Value*Variables(9)).*Variables(5).*(Parameters.CF1.Value*TotF1./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(2)./TotF1+Parameters.CF2.Value*TotF2./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(4)./TotF2)-Parameters.MuMm1.Value^-1*Variables(5))*TStep;
            TempVariables(6) = TempVariables(6) + ( Parameters.BetaFM.Value*(1-Parameters.eHIV.Value*Variables(9)).*Variables(5).*(Parameters.CF1.Value*TotF1./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(2)./TotF1+Parameters.CF2.Value*TotF2./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(4)./TotF2) - (Parameters.MuMm1.Value^-1+Parameters.Alpham1.Value^-1)*Variables(6))*TStep;
            TempVariables(7) = (TempVariables(2) + TempVariables(4))/(TotF1+TotF2)*100; 
            TempVariables(8) = (TempVariables(6))/(TotM)*100;         
            TempVariables(2) = min(TotF1,TempVariables(2));
            TempVariables(4) = min(TotF2,TempVariables(4));
            TempVariables(6) = min(TotM,TempVariables(6));
            Variables = TempVariables;
            Record(IndBeta,:) = Variables;
            IndBeta = IndBeta+1;
        end    
        LogLik = LogLik + log(normpdf(Variables(Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),Parameters.MultNoise*Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))');
    end
    res.LogLik = LogLik;
    res.Record = Record;
   