function Result = MIFHIV(Data,Parameters,Model)

Result = struct();
Result.SamplesX = [];
Result.SamplesTheta = [];

Observations = Data.Observations;
ObservationInstants = Data.Instants;
ComputationTStep = Parameters.ComputationTStep;
ObservationTStep = mean(diff(ObservationInstants));
NbTSteps = ceil(ObservationTStep/ComputationTStep);
TStep = ComputationTStep;
NbIterations = Parameters.MIFNbIterations;
NbParticules = Parameters.MIFNbParticules;
a = Parameters.MIFCoolingParameters;
b = Parameters.MIFb;

ThetasRecord = [];

% Initialize Theta
NbPars = length(Parameters.Names.Estimated); 
NbParsInit = length(Parameters.Names.EstimatedInit);
NbParsNotInit = length(Parameters.Names.EstimatedNotInit); 
Names = Parameters.Names.Estimated;
NamesInit = Parameters.Names.EstimatedInit;
NamesNotInit = Parameters.Names.EstimatedNotInit;
ThetasInit = ones(NbParticules,NbParsInit);
ThetasNotInit = ones(NbParticules,NbParsNotInit);
InitIndices = [];
NotInitIndices = [];
for i = 1:NbParsNotInit
    ind = Parameters.(NamesNotInit{i}).Index;
    NotInitIndices(i) = ind;
    ThetasNotInit(:,i) = Parameters.(NamesNotInit{i}).TransfValue*ThetasNotInit(:,i);
end
for i = 1:NbParsInit
    ind = Parameters.(NamesInit{i}).Index;
    InitIndices(i) = ind;
    ThetasInit(:,i) = Parameters.(NamesInit{i}).TransfValue*ThetasInit(:,i);
end
ThetasRecord(InitIndices,1) = mean(ThetasInit)';
ThetasRecord(NotInitIndices,1) = mean(ThetasNotInit)';
ThetasRecordTemp = zeros(length(ObservationInstants),NbParsNotInit);
ThetasVarRecordTemp = zeros(length(ObservationInstants),NbParsNotInit);
LogLiks = [];


for IndIt = 2:NbIterations

    disp(['Iteration ' num2str(IndIt)])
    
    ThetasInit = repmat(ThetasRecord(InitIndices,IndIt-1),1,NbParticules)';
    ThetasNotInit = repmat(ThetasRecord(NotInitIndices,IndIt-1),1,NbParticules)';
    
    
    ThetasInit = ThetasInit + a^(IndIt-1)*b*((chol(Parameters.MIFCov(InitIndices,InitIndices)))*(randn(NbParticules,NbParsInit)'))';
    ThetasNotInit = ThetasNotInit + a^(IndIt-1)*b*((chol(Parameters.MIFCov(NotInitIndices,NotInitIndices)))*(randn(NbParticules,NbParsNotInit)'))';
    ThetasRecordTemp(1,:) = mean(ThetasNotInit) ; 
    ThetasVarRecordTemp(1,:) = var(ThetasNotInit) ; 
    RecordStates = zeros(length(ObservationInstants),NbParticules,Parameters.NbVariables-1);
    RecordThetas = zeros(length(ObservationInstants),NbParticules,NbParsNotInit);
    LogLiksTemp = zeros(1,length(ObservationInstants));
    NoTransfThetasInit = UpdateParsTransfToNoTransf(Parameters,ThetasInit,NamesInit);
    NoTransfThetasNotInit = UpdateParsTransfToNoTransf(Parameters,ThetasNotInit,NamesNotInit);
    
    
    InitialSF1 = (1-NoTransfThetasInit(:,Parameters.InitialIPropF.IndexInit))*Parameters.TotF1.Value;
    InitialHIVF1 = NoTransfThetasInit(:,Parameters.InitialIPropF.IndexInit)*Parameters.TotF1.Value;
    InitialSF2 = (1-NoTransfThetasInit(:,Parameters.InitialIPropF.IndexInit))*Parameters.TotF2.Value;
    InitialHIVF2 = NoTransfThetasInit(:,Parameters.InitialIPropF.IndexInit)*Parameters.TotF2.Value;
    TotF1 =  InitialSF1 + InitialHIVF1;
    TotF2 =  InitialSF2 + InitialHIVF2;
    TotM = Parameters.TotalFSW.Value*NoTransfThetasInit(:,Parameters.TotMFactor.IndexInit);
    InitialSM = (1-NoTransfThetasInit(:,Parameters.InitialIPropM.IndexInit)).*TotM;
    InitialHIVM = NoTransfThetasInit(:,Parameters.InitialIPropM.IndexInit).*TotM;
   
    Xf(:,1) = InitialSF1;
    Xf(:,2) = InitialHIVF1;
    Xf(:,3) = InitialSF2;
    Xf(:,4) = InitialHIVF2;
    Xf(:,5) = InitialSM;
	Xf(:,6) = InitialHIVM;
    Xf(:,7) = (Xf(:,2)+Xf(:,4))./(sum(Xf(:,1:4)')')*100;
    Xf(:,8) = (Xf(:,6))./(sum(Xf(:,5:6)')')*100;
    Xf(:,9) = NoTransfThetasInit(:,Parameters.InitialFt.IndexInit);
    

    for IndTime = 2:length(ObservationInstants)
        
        BetaFM = 1-(1-NoTransfThetasNotInit(:,Parameters.BetaFMPerAct.IndexNotInit)).^NoTransfThetasNotInit(:,Parameters.NumberActsPerClient.IndexNotInit);
        BetaMF = 1-(1-NoTransfThetasNotInit(:,Parameters.BetaMFPerAct.IndexNotInit)).^NoTransfThetasNotInit(:,Parameters.NumberActsPerClient.IndexNotInit);

        
        TempXf = Xf;
        rands = randn(NbParticules,Data.NbComputingSteps(IndTime));
        
        MuF = NoTransfThetasNotInit(:,Parameters.MuFm1.IndexNotInit).^-1;
        MuM = NoTransfThetasNotInit(:,Parameters.MuMm1.IndexNotInit).^-1;
        Alpha = NoTransfThetasNotInit(:,Parameters.Alpham1.IndexNotInit).^-1;
        sigmasRW = NoTransfThetasNotInit(:,Parameters.SigmaRW.IndexNotInit);
        CF1 = NoTransfThetasNotInit(:,Parameters.CF1.IndexNotInit);
        CF2 = NoTransfThetasNotInit(:,Parameters.CF2.IndexNotInit);
        eHIV = NoTransfThetasNotInit(:,Parameters.eHIV.IndexNotInit);
        

        for IndDiscr = 1:Data.NbComputingSteps(IndTime)
            
            
            TempXf(:,1) = TempXf(:,1) + ( MuF.*(TotF1) + Alpha.*Xf(:,2) - BetaMF.*CF1.*(1-eHIV.*Xf(:,9)).*Xf(:,1).*Xf(:,6)./TotM-MuF.*Xf(:,1))*TStep;
            TempXf(:,2) = TempXf(:,2) + ( BetaMF.*CF1.*(1-eHIV.*Xf(:,9)).*Xf(:,1).*Xf(:,6)./TotM - (MuF + Alpha).*Xf(:,2))*TStep;
            TempXf(:,3) = TempXf(:,3) + ( MuF.*(TotF2) + Alpha.*Xf(:,4) - BetaMF.*CF2.*(1-eHIV.*Xf(:,9)).*Xf(:,3).*Xf(:,6)./TotM-MuF.*Xf(:,3))*TStep;
            TempXf(:,4) = TempXf(:,4) + ( BetaMF.*CF2.*(1-eHIV.*Xf(:,9)).*Xf(:,3).*Xf(:,6)./TotM - (MuF + Alpha).*Xf(:,4))*TStep;
            TempXf(:,5) = TempXf(:,5) + ( MuM.*(TotM) + Alpha.*Xf(:,6) - BetaFM.*(1-eHIV.*Xf(:,9)).*Xf(:,5).*(CF1.*TotF1./(CF2.*TotF2+CF1.*TotF1).*Xf(:,2)./TotF1+CF2.*TotF2./(CF2.*TotF2+CF1.*TotF1).*Xf(:,4)./TotF2)-MuM.*Xf(:,5))*TStep;
            TempXf(:,6) = TempXf(:,6) + ( BetaFM.*(1-eHIV.*Xf(:,9)).*Xf(:,5).*(CF1.*TotF1./(CF2.*TotF2+CF1.*TotF1).*Xf(:,2)./TotF1+CF2.*TotF2./(CF2.*TotF2+CF1.*TotF1).*Xf(:,4)./TotF2) - (MuM+Alpha).*Xf(:,6))*TStep;
            TempXf(:,7) = (TempXf(:,2) + TempXf(:,4))./(TotF1+TotF2)*100; 
            TempXf(:,8) = (TempXf(:,6))./(TotM)*100; 
            TempXf(:,9) = min(1,max(0,TempXf(:,9) + sqrt(TStep)*sigmasRW.*rands(:,IndDiscr)));

            Xf = TempXf;
            
        end
        
        Weigths = ones(NbParticules,1);
        for i = 1:length(Names)
            if Parameters.(Names{i}).Init
                temptransf = ThetasInit(:,Parameters.(Names{i}).IndexInit);
                tempnotransf = NoTransfThetasInit(:,Parameters.(Names{i}).IndexInit);
            else
                temptransf = ThetasNotInit(:,Parameters.(Names{i}).IndexNotInit);
                tempnotransf = NoTransfThetasNotInit(:,Parameters.(Names{i}).IndexNotInit);
            end
            tmp = Parameters.(Names{i}).Prior(Names{i},Parameters,temptransf,tempnotransf);
            Weigths =  Weigths.*(tmp);
        end
        
        LogLikWeigths = normpdf(Xf(:,Data.ObservedVariables(:,IndTime)),Observations(Data.ObservedVariables(:,IndTime), IndTime)*ones(NbParticules,1),Parameters.SigmaObs*Observations(Data.ObservedVariables(:,IndTime),IndTime)*ones(NbParticules,1));
        Weigths = Weigths.*LogLikWeigths;
        ResampledIndexes = Resample(Weigths);
        if length(unique(ResampledIndexes))<NbParticules/10
            disp(['Nb Kept Particules : ' num2str(length(unique(ResampledIndexes)))]) 
        end
        Xf(:,1) = Xf(ResampledIndexes,1);
        Xf(:,2) = Xf(ResampledIndexes,2);
        Xf(:,3) = Xf(ResampledIndexes,3);
        Xf(:,4) = Xf(ResampledIndexes,4);
        Xf(:,5) = Xf(ResampledIndexes,5);
        Xf(:,6) = Xf(ResampledIndexes,6);
        Xf(:,7) = Xf(ResampledIndexes,7);
        Xf(:,8) = Xf(ResampledIndexes,8);
        Xf(:,9) = Xf(ResampledIndexes,9);
        
        ThetasNotInit = ThetasNotInit(ResampledIndexes,:) + a^(IndIt-1)*((chol(Parameters.MIFCov(NotInitIndices,NotInitIndices)))*(randn(NbParticules,NbParsNotInit)'))';
        ThetasInit = ThetasInit(ResampledIndexes,:);
        NoTransfThetasInit = UpdateParsTransfToNoTransf(Parameters,ThetasInit,NamesInit);
        NoTransfThetasNotInit = UpdateParsTransfToNoTransf(Parameters,ThetasNotInit,NamesNotInit);

        NoTransfThetasRecordTemp(IndTime,:) = mean(NoTransfThetasNotInit)';
        ThetasRecordTemp(IndTime,:) = mean(ThetasNotInit)';
        ThetasVarRecordTemp(IndTime,:) = var(ThetasNotInit)';        
        RecordStates(IndTime,:,:) = Xf;
        RecordThetas(IndTime,:,:) = ThetasNotInit;
        LogLiksTemp(IndTime) = log(sum(LogLikWeigths)/NbParticules);
        if or(not(isreal(log(sum(Weigths)/NbParticules))),or(isinf(log(sum(Weigths)/NbParticules)),isnan(log(sum(Weigths)/NbParticules))))
            disp('inf or nan')
        end
    end
%     plot(mean(squeeze(RecordStates(:,:,5))'))
%     hold on
%     plot(Data.Observations(5,:),'g')
%     hold off
%     die
    temp = zeros(1,NbParsNotInit);
    for IndTime = 2:length(ObservationInstants)
        temp = temp + min(ThetasVarRecordTemp).*(ThetasVarRecordTemp(IndTime,:).^-1).*(ThetasRecordTemp(IndTime,:)-ThetasRecordTemp(IndTime-1,:));
    end
    ThetasRecord(NotInitIndices,IndIt) = ThetasRecord(NotInitIndices,IndIt-1) + temp';
    ThetasRecord(InitIndices,IndIt) = mean(ThetasInit)';
    tmp = [ThetasRecord(:,IndIt) ThetasRecord(:,IndIt)];
    NoTransfThetasRecord(:,IndIt) = UpdateParsTransfToNoTransf(Parameters,ThetasRecord(:,IndIt)',Names);
    
    LogLiks(IndIt-1) = sum(LogLiksTemp);
    TempParameters = Parameters;
    for i = 1:length(Names)
        ind  = Parameters.(Names{i}).Index;
        TempParameters.(Names{i}).TransfValue = ThetasRecord(ind,IndIt);
    end
    TempParameters = UpdateParsTransfToNoTransf(TempParameters);
%     Temp = EstimationSMCsmoothGen(Data, Model, TempParameters);
%     figure(2)
%     for i = 1:length(NamesNotInit)
%         subplot(length(NamesNotInit),1,i)
%         plot(NoTransfThetasRecordTemp(2:end,i))
%         title(NamesNotInit{i})
%         hold on
%         xis = ones(size(NoTransfThetasRecordTemp(:,i)));
%         plot(Parameters.(NamesNotInit{i}).MinLim*xis,'--r')
%         plot(Parameters.(NamesNotInit{i}).MaxLim*xis,'--r')
%         hold off
%     end
    

    subplot(3,1,1)
    plot(cumsum(Data.NbComputingSteps),mean(squeeze(RecordStates(:,:,7))'))
    hold on
    for i = 2:length(Data.NbComputingSteps)
        if Data.ObservedVariables(i) == 7
            plot(sum(Data.NbComputingSteps(1:i)),Data.Observations(Data.ObservedVariables(i),i),'og')
        end
    end
    hold off
    title('Sex Workers')

    subplot(3,1,2)
    plot(cumsum(Data.NbComputingSteps),mean(squeeze(RecordStates(:,:,8))'))
    hold on
    for i = 2:length(Data.NbComputingSteps)
        if Data.ObservedVariables(i) == 8
            plot(sum(Data.NbComputingSteps(1:i)),Data.Observations(Data.ObservedVariables(i),i),'og')
        end
    end
    hold off
    title('Clients')

    subplot(3,1,3)
    plot(cumsum(Data.NbComputingSteps),mean((squeeze(RecordStates(1:end,:,9))')))
    title(['Condom use ' num2str(sum(LogLiksTemp))])
    
%     figure(1)
%     PlotMarc(Temp,8)
%     title(Temp.LogLik)
    pause(0.05)
end

figure(1)
subplot(4,1,1)
plot(LogLiks(2:end))
title('LogLiks')
subplot(4,1,2)
plot(cumsum(Data.NbComputingSteps),mean(squeeze(RecordStates(:,:,7))'))
hold on
for i = 2:length(Data.NbComputingSteps)
    if Data.ObservedVariables(i) == 7
        plot(sum(Data.NbComputingSteps(1:i)),Data.Observations(Data.ObservedVariables(i),i),'og')
    end
end
hold off
title('Sex Workers')

subplot(4,1,3)
plot(cumsum(Data.NbComputingSteps),mean(squeeze(RecordStates(:,:,8))'))
hold on
for i = 2:length(Data.NbComputingSteps)
    if Data.ObservedVariables(i) == 8
        plot(sum(Data.NbComputingSteps(1:i)),Data.Observations(Data.ObservedVariables(i),i),'og')
    end
end
hold off
title('Clients')

subplot(4,1,4)
plot(cumsum(Data.NbComputingSteps),mean((squeeze(RecordStates(1:end,:,9))')))
title('Condom use')

figure(2)
for i = 1:size(ThetasRecord,1)
    subplot(size(ThetasRecord,1),1,i)
    plot(NoTransfThetasRecord(i,2:end))
    title(Names{i})
end

Result.NoTransfThetasRecord = NoTransfThetasRecord;
Result.ThetasRecord = ThetasRecord;
Result.RecordStates = RecordStates;
Result.LogLiks = LogLiks;





