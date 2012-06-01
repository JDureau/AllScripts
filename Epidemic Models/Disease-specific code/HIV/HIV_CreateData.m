function Data = HIV_CreateData(Fts,Parameters,Model,Data)
    
    % creates data from a Beta trajectory

    Parameters = Model.InitializeParameters(Parameters);
    Variables = Parameters.InitialState;


    TempVariables = Variables;


    TStep = Parameters.ComputationTStep;
    NbTSteps = Parameters.NbTSteps;
    
    IndBeta = 1;
    Record = [];
    crash = 0;
    if mean(not(isreal(Fts)))
        crash = 1;
        Fts = Fts(1)*ones(size(Fts));
    end
            
    for IndDiscr = 1:NbTSteps
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
        TempVariables(5) = TempVariables(5) + ( Parameters.MuMm1.Value^-1*(TotM) + Parameters.Alpham1.Value^-1*Variables(6) - Parameters.BetaFM.Value*Parameters.NbContactsForMen*(1-Parameters.eHIV.Value*Variables(9)).*Variables(5).*(Parameters.CF1.Value*TotF1./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(2)./TotF1+Parameters.CF2.Value*TotF2./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(4)./TotF2)-Parameters.MuMm1.Value^-1*Variables(5))*TStep;
        TempVariables(6) = TempVariables(6) + ( Parameters.BetaFM.Value*Parameters.NbContactsForMen*(1-Parameters.eHIV.Value*Variables(9)).*Variables(5).*(Parameters.CF1.Value*TotF1./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(2)./TotF1+Parameters.CF2.Value*TotF2./(Parameters.CF2.Value*TotF2+Parameters.CF1.Value*TotF1).*Variables(4)./TotF2) - (Parameters.MuMm1.Value^-1+Parameters.Alpham1.Value^-1)*Variables(6))*TStep;
        TempVariables(7) = (TempVariables(2) + TempVariables(4))/(TotF1+TotF2)*100; 
        TempVariables(8) = (TempVariables(6))/(TotM)*100;         
        TempVariables(2) = min(TotF1,TempVariables(2));
        TempVariables(4) = min(TotF2,TempVariables(4));
        TempVariables(6) = min(TotM,TempVariables(6));
        Variables = TempVariables;
        Record(IndBeta,:) = Variables;
        IndBeta = IndBeta+1;
    end    

    
   
%     TellParsValues(Parameters)
    
    Data.Parameters = Parameters;
    Data.BuiltTraj = Record;
    Data.Observations = [];
    for i = 2:length(Data.Instants)
        inds = Data.ObservedVariables(i);
        if Parameters.MultNoise
            Data.Observations(inds,i) = binornd(425,Record(sum(Data.NbComputingSteps(1:i)),inds)/100)/4.25;  
        else
            Data.Observations(inds,i) = Record(sum(Data.NbComputingSteps(1:i)),inds)*(1 + randn(1,1)*Parameters.MultNoise);
        end
    end
    
    baseline = Parameters.BRbase.Value;
    m = Parameters.BRmm1.Value+1;
    mu = Parameters.BRmu.Value;

    
    subplot(3,1,1)
    xis = TStep:Parameters.ComputationTStep:NbTSteps*TStep;
    plot(xis,Record(:,7))
    hold on
    try
        for i = 2:length(Data.Instants)
            if Data.ObservedVariables(i)==7
                plot(sum(Data.NbComputingSteps(1:i))*Parameters.ComputationTStep,Data.Observations(7,i),'or')
            end
        end
    end
    hold off
    title('FSWs')
    
    subplot(3,1,2)
    plot(xis,Record(:,8))
    hold on
    try
        for i = 2:length(Data.Instants)
            if Data.ObservedVariables(i)==8
                plot(sum(Data.NbComputingSteps(1:i))*Parameters.ComputationTStep,Data.Observations(8,i),'or')
            end
        end
    end
    hold off
    title('Clients')
    
    subplot(3,1,3)
    plot(xis,Fts(1:length(xis)))
    title('Condom Use')
    ylim([0 1])
    
    
        
%     if mean(not(isreal(Fts))
%         crash = 1;
%         Fts = Fts(1)*ones(size(Fts));
%     end
%     
    Data.Crash = crash;
    Data.Fts = Fts;
    
    
    
    