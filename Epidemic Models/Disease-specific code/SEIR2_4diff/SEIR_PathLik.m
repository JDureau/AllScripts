function Temp = SEIR_PathLik(Data, Model,Parameters)

    Parameters = Model.InitializeParameters(Parameters);

    Variables = Parameters.InitialState';
    gamma = Parameters.gammam1.Value^-1;
    k = Parameters.km1.Value^-1;
    TotPop = Parameters.TotalPopulation;
    TempVariables = Variables;
    TempVariables(:,5) = zeros(size(TempVariables(:,5)));
    TempPath = Parameters.CurrentPaths;
    Betas = exp(TempPath(6,:));
    RecordVariables = zeros(Parameters.NbVariables,length(Data.Instants));
    ComputationTStep = Parameters.ComputationTStep;
    


    
    LogLik = 0;
    ind = 0;
    LogLiks = [];
    for IndTime = 2:length(Data.Instants)
        TempVariables(5) = 0;
        for IndDiscr = 1:Data.NbComputingSteps(IndTime)
            ind = ind+1;
            % Variables
            beta = Betas(ind);
            
            TempVariables(1) = TempVariables(1) + (-beta.*Variables(1).*Variables(3)/TotPop)*ComputationTStep ;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
            TempVariables(2) = TempVariables(2) + ( beta.*Variables(1).*Variables(3)/TotPop-k*Variables(2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands(:,IndDiscr,2);
            TempVariables(3) = TempVariables(3) + (-gamma*Variables(3) + k*Variables(2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(IndDiscr,3);
            TempVariables(4) = TempVariables(4) + ( gamma*Variables(3))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(IndDiscr,3);
            TempVariables(5) = TempVariables(5) + ( k*Variables(2))*ComputationTStep ;
            TempVariables(1) = max(TempVariables(1),0);
            TempVariables(2) = max(TempVariables(2),0);
            TempVariables(3) = max(TempVariables(3),0);
            TempVariables(6) = log(beta);
            Variables = TempVariables;
            
        end
        RecordVariables(:,IndTime) = Variables;
        try
            coeff = Parameters.MultCoeff.Value/10;
        catch
            coeff = 1;
        end
        temp =  log(max(0,eval(Model.LikFunction)));
        LogLiks(IndTime) = temp;
        LogLik = LogLik + temp;
    end
    
    subplot(2,1,1)
    plot(Data.Observations(5,:),'g')
    hold on
    plot(RecordVariables(5,:))
    hold off
    subplot(2,1,2)
    plot(RecordVariables(6,2:end))
    hold on
    plot(log(Parameters.Betas(Data.Instants(2:end)+1)),'g')
    hold off
    pause(0.01)
%     disp(LogLik)
%     plot(Betas)
%     ylim([0.5 2.5])
%     pause(0.01)
    Temp.CompletePaths = RecordVariables;
    Temp.LogLik = LogLik;
    Temp.RecordVariables = RecordVariables;
%     Temp.TimeVarPath = Path;
%     Temp.RtPath = gamma^-1*exp(Path).*RecordVariables(1,:)/Parameters.TotalPopulation;
%     temp = cumsum(Data.NbComputingSteps(2:end));
%     Temp.IncEst = RecordVariables(5,temp);
