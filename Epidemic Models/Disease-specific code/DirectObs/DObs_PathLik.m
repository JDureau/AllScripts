function Temp = DObs_PathLik(Data, Model,Parameters)

    Parameters = Model.InitializeParameters(Parameters);

    Variables = Parameters.InitialState';
   
    TempVariables = Variables;

    TempPath = Parameters.CurrentPaths;
    Betas = TempPath(1,:);
    RecordVariables = zeros(Parameters.NbVariables,length(Data.Instants));
    ComputationTStep = Parameters.ComputationTStep;
    


    
    
    LogLik = 0;
    ind = 1;
    LogLiks = [];
    for IndTime = 2:length(Data.Instants)

        for IndDiscr = 1:Data.NbComputingSteps(IndTime)
            ind = ind+1;
            % Variables
            beta = Betas(ind);
            
            TempVariables(1) = beta;
            Variables = TempVariables;
            
        end
        RecordVariables(:,IndTime) = Variables;
        temp =  log(max(eps,eval(Model.LikFunction)));
        LogLiks(IndTime) = temp;
        LogLik = LogLik + temp;
    end
   
%     disp(LogLik)
%     plot(Betas)
%     ylim([0.5 2.5])
%     pause(0.01)

    subplot(2,1,1)
    plot(Data.Observations(1,:),'g')
    hold on
    plot(RecordVariables(1,:))
    hold off
    subplot(2,1,2)
    plot(RecordVariables(1,2:end))
    hold on
    plot((Parameters.Betas(Data.Instants(2:end)+1)),'g')
    hold off
    pause(0.01)

    Temp.CompletePaths = RecordVariables;
    Temp.LogLik = LogLik;
    Temp.RecordVariables = RecordVariables;
%     Temp.TimeVarPath = Path;
%     Temp.RtPath = gamma^-1*exp(Path).*RecordVariables(1,:)/Parameters.TotalPopulation;
%     temp = cumsum(Data.NbComputingSteps(2:end));
%     Temp.IncEst = RecordVariables(5,temp);
